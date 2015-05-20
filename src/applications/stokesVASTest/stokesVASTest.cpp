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
    using namespace SEM::field;
    using namespace SEM::mesh;
    using iomanagment::NO_READ;
    using iomanagment::NO_WRITE;
    using iomanagment::AUTO;
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
    
    //--------------------- CREATE FIELDS -----------------------------------
    VectorField U("U", mesh);
    VectorField Residue("Residue", mesh, NO_READ, NO_WRITE);
    VectorField Usol("Usol", mesh, NO_READ);
    VectorField Uerr("Uerr", mesh, NO_READ);
    
    ScalarField p("p", mesh,NO_READ);
    ScalarField p_real("p_real", mesh, NO_READ);
    ScalarField psol("psol", mesh, NO_READ);
    ScalarField perr("perr", mesh, NO_READ);
    ScalarField p_realErr("p_realErr", mesh, NO_READ);
    ScalarField p_ext("p_ext",mesh, NO_READ,NO_WRITE);
    
    VectorField force("force",mesh, NO_READ, NO_WRITE);
    
    //--------------------- GET COEFFS -----------------------------------
    std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( Case::solutionControl().timeDiscretization() );
    Scalar dt = Case::time().timeStep();
    
    auto X=xCmps(mesh.spectralNodes());
    auto Y=yCmps(mesh.spectralNodes());

    Scalar Pi = 4.*atan(1.);
    Scalar nu = Case::material()["nu"];
    
    using SEM::array::sin;
    using SEM::array::cos;
    using SEM::array::abs;
    using SEM::array::mag;
    using SEM::array::max;
    
    numArray<Scalar> sx = sin(Pi*X);
    numArray<Scalar> s2x = sin(2.*Pi*X);
    numArray<Scalar> sy = sin(Pi*Y);
    numArray<Scalar> s2y = sin(2.*Pi*Y);
    numArray<Scalar> cx = cos(Pi*X);
    numArray<Scalar> cy = cos(Pi*Y);
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    //initial conditions:
    xCmps(U) = Pi*s2y*sx*sx;
    yCmps(U) =-Pi*s2x*sy*sy;
    U |= "fixedValue";
    U |= Vector(0.,0.);
    U.registryFile()->fireWriting();
    
    p =-cx*sy; 
    p.registryFile()->fireWriting();
    
    
    Usol = U;
    Usol.registryFile()->fireWriting();
    
    psol = p;
    psol.registryFile()->fireWriting();
    
    DiscontinousField<Vector> force2(mesh);
    
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar t = Case::time().time();
        
        Scalar st = std::sin(t);
        Scalar ct = std::cos(t);
        
        xCmps(force) = -Pi*st*s2y*sx*sx + Pi*ct*sx*sy - 2.*Pi*Pi*Pi*ct*s2y*(cx*cx-sx*sx) + 4.*Pi*Pi*Pi*ct*s2y*sx*sx;
        yCmps(force) =  Pi*st*s2x*sy*sy - Pi*ct*cx*cy + 2.*Pi*Pi*Pi*ct*s2x*(cy*cy-sy*sy) - 4.*Pi*Pi*Pi*ct*s2x*sy*sy;
        
        xCmps(Usol) = Pi*ct*s2y*sx*sx;
        yCmps(Usol) =-Pi*ct*s2x*sy*sy;
        psol        =-ct*cx*sy; 
        
        
        p_ext = extInTime(p);
        
        force2 = force;
        
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi)) == a(force2,phi) + a(p_ext,grad(phi))
        );
        
        las::solve
        (
            a(grad(p),grad(phi)) == a(force2,grad(phi)) + bInt_rotUdotNxGradPhi(U,nu) - bInt_ddtUdotN(U)
        );
        
        p_real = p - nu*div(U);
        
        //fix pressure to one exact solution node
        size_t n = psol.size()/2;
        p += ( psol[n]-p[n] );
        p_real += (psol[n]-p_real[n]);
        
        
        p_realErr =abs(psol-p_realErr)/(1.+abs(psol));
        perr = abs(psol-p)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        
        Case::time().fireWriting();
    }
    
    return 0;
}