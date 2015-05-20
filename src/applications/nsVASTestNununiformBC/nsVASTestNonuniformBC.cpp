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
#include "fieldMath/BIntegralAvgPressure.h"
#include "fieldMath/BIntegral_ddtUdotN.h"
#include "fields/InletOutletVelocity.h"
#include "fieldMath/NonlinearOperatorSplitting.h"
#include "fieldMath/TimeExtrapolation.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using SEM::iomanagment::NO_READ;
    using SEM::iomanagment::NO_WRITE;
    using SEM::weak::phi;
    using SEM::weak::ddt;
    using SEM::field::grad;
    using SEM::weak::ddn;
    using SEM::weak::grad;
    using SEM::weak::a;
    using SEM::weak::bInt_rotUdotNxGradPhi;
    using SEM::weak::bInt_UdotNPhi;
    using SEM::bInt_avgPressure;
    using SEM::bInt_ddtUdotN;
    
    using SEM::array::sin;
    using SEM::array::cos;
    using SEM::array::abs;
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
//     //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
    VectorField Usol("Usol",mesh,NO_READ);
    VectorField Uerr("Uerr",mesh,NO_READ);
    ScalarField p("p",mesh, READ_IF_PRESENT);
    ScalarField psol("psol",mesh,NO_READ);
    ScalarField perr("perr", mesh,NO_READ);
    
    ScalarField p_ext("p_ext",mesh,NO_READ,NO_WRITE);
    VectorField force("force",mesh, NO_READ, NO_WRITE);
    
    ScalarField divU("divU", mesh, NO_READ);
    //--------------------- GET COEFFS -----------------------------------
    
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    U.registryFile()->fireWriting();
    
    auto X = xCmps(mesh.spectralNodes());
    auto Y = yCmps(mesh.spectralNodes());
    
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar t=Case::time().time();
        
        //--setup assumed solution--------------
        xCmps(Usol) = -cos(X)*sin(Y)*std::sin(2.*t);
        yCmps(Usol) = sin(X)*cos(Y)*std::sin(2.*t);
        
        psol = -0.25*(cos(2.*X)+cos(2.*Y))*(std::sin(2.*t)*std::sin(2.*t));
        
        xCmps(force) = -cos(X)*sin(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        yCmps(force) = sin(X)*cos(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        DiscontinousField<Vector> F = force;
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(U);
        
        p_ext = extInTime(p);
        
        DiscontinousField<Vector> oldTimeTerms(mesh);
        oldTimeTerms = 0.;
        for(unsigned int i=1; i<timeCoeffs.size();++i)
        {
            oldTimeTerms += timeCoeffs[i]/dt*oifsConvDerivative(U,dt,i-1);
        }
        
        las::solve
        (
            a(timeCoeffs[0]/dt,U,phi) + a(nu*grad(U),grad(phi)) == a(p_ext,grad(phi)) + a(F,phi) - a(oldTimeTerms,phi) - bInt_avgPressure(p,U)
        );
        
        
        las::solve
        (
            a(grad(p),grad(phi)) == a(F-cDeriv(U,U),grad(phi)) + bInt_rotUdotNxGradPhi(U,nu) - bInt_ddtUdotN(U)
        );

        
        divU=div(U);
        
        size_t n = X.size()/2;
        p += (psol[n]-p[n]);
        
        perr = abs(psol-p)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        
        Case::time().fireWriting();
        
    }
    
    return 0;
}