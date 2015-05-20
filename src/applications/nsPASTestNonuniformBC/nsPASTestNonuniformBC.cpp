#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "fieldMath/f.h"
#include "fieldMath/postCalculation.h"
#include "fieldMath/weakGradV.h"

#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "utilities/VectorUtils.h"
#include "fields/DiscontinousField.h"
#include "fieldMath/WeakFormDefinition.h"
#include "fieldMath/DivergenceOperator.h"
#include "fieldMath/GradientOperator.h"
#include "fieldMath/weakDivW.h"
#include "fieldMath/BoundaryIntegralNRotRotU.h"
#include "fieldMath/ConvectiveDerivative.h"
#include "fieldMath/NonlinearOperatorSplitting.h"
#include "fieldMath/BIntegral_ddtUdotN.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using SEM::weak::phi;
    using SEM::weak::ddt;
    using SEM::field::grad;
    using SEM::weak::ddn;
    using SEM::weak::grad;
    using SEM::weak::a;
    using SEM::weak::bInt_rotUdotNxGradPhi;
    using SEM::weak::bInt_UdotNPhi;
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
//     //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
    VectorField Usol("Usol",mesh, NO_READ);
    VectorField Uerr("Uerr",mesh, NO_READ);
    VectorField force("force",mesh, NO_READ);
    
    ScalarField p("p",mesh, NO_READ);
    ScalarField p_real("p_real",mesh,NO_READ);
    ScalarField psol("psol",mesh,NO_READ);
    ScalarField perr("perr",mesh,NO_READ);
    ScalarField p_realErr("p_realErr",mesh,NO_READ);

    ScalarField divU("divU",mesh,NO_READ);
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    //initial conditions:
    U = Vector(0,0);
    U.registryFile()->fireWriting();
    
    p = 0.;
    p.registryFile()->fireWriting();
    
    p_real = p;
    p_real.registryFile()->fireWriting();
    
    DiscontinousField<Vector> K=cDeriv(U,U);
    DiscontinousField<Vector> Kold(K);
    DiscontinousField<Vector> Kold2(K);
    
    auto X = xCmps(mesh.spectralNodes());
    auto Y = yCmps(mesh.spectralNodes());
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar t = Case::time().time();
        
        xCmps(Usol) = -cos(X)*sin(Y)*std::sin(2.*t);
        yCmps(Usol) = sin(X)*cos(Y)*std::sin(2.*t);
        
        psol = -0.25*(cos(2.*X)+cos(2.*Y))*(std::sin(2.*t)*std::sin(2.*t));
        
        xCmps(force) = -cos(X)*sin(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        yCmps(force) = sin(X)*cos(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        DiscontinousField<Vector> F = force;
        
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(U);
        
        DiscontinousField<Vector> Residue = -(2.*K-Kold) + F;//-(3.*K - 3.*Kold + Kold2) + F;
        
        for(unsigned int i=1; i<timeCoeffs.size(); ++i)
        {
            Residue-=timeCoeffs[i]/dt*U.oldField(i-1);
        }
        las::solve
        (
            a(grad(p),grad(phi)) == a(Residue,grad(phi))  + bInt_rotUdotNxGradPhi(U,nu,true) - bInt_UdotNPhi(U,timeCoeffs[0]/dt)
        );
        
        // Velocity equation
        DiscontinousField<Vector> oldTimeTerms(mesh);
        for(unsigned int i=1; i<timeCoeffs.size();++i)
        {
            oldTimeTerms += timeCoeffs[i]/dt*oifsConvDerivative(U,dt,i-1);
        }
        las::solve
        (
            a(timeCoeffs[0]/dt,U,phi) + a(nu*grad(U),grad(phi)) == a(p,grad(phi)) + a(F,phi) - a(oldTimeTerms,phi)
        );
        
        Kold2=Kold;
        Kold=K;
        K = cDeriv(U,U);
        
        p_real = p - nu*div(U);
        divU = div(U);
        
        // ---- fix pressure solution, because it is solved with respect to constatn
        size_t n = X.size()/2;
        p_real += (psol[n]-p_real[n]);
        p += (psol[n]-p[n]);
        
        // --- calculate error ---------------
        perr = abs(psol-p)/(1.+abs(psol));
        p_realErr = abs(psol-p_real)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        
        Case::time().fireWriting();
    }
    
    return 0;
}