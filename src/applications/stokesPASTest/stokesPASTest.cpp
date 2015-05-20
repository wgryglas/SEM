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
    
    //--------------------- CREATE FIELDS -----------------------------------
    VectorField U("U",mesh);
    VectorField Residue("Residue",mesh, NO_READ,NO_WRITE);
    VectorField Usol("Usol",mesh, NO_READ);
    VectorField Uerr("Uerr",mesh, NO_READ);
    ScalarField p("p",mesh, NO_READ);
    ScalarField p_real("p_real",mesh, NO_READ);
    ScalarField psol("psol",mesh, NO_READ);
    ScalarField perr("perr",mesh, NO_READ);
    ScalarField p_realError("p_realError",mesh, NO_READ);
    
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    auto X=xCmps(mesh.spectralNodes());
    auto Y=yCmps(mesh.spectralNodes());

    numArray<Vector> force(X.size());

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
    
    //--------------------- INITIAL CONDITIONS -----------------------------------
    //initial conditions:
    xCmps(U) = Pi*s2y*sx*sx;
    yCmps(U) =-Pi*s2x*sy*sy;
    U |= "fixedValue";
    U |= Vector(0.,0.);
    
    p_real =-cx*sy; 
    
    p=p_real;
    
    Usol = U;
    
    psol = p;
    

    Case::time().fireWriting();


    //--------------------- SOLVE EQUATION -----------------------------------

    
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
        
        
        Residue = force;
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( U );
        for(int i=1; i<timeCoeffs.size(); ++i)
        {
            Residue-=timeCoeffs[i]/dt*U.oldField(i-1);
        }
        
        las::solve
        (
            a(grad(p),grad(phi)) == a(Residue,grad(phi)) + bInt_rotUdotNxGradPhi(U,nu,true) - bInt_UdotNPhi(U,timeCoeffs[0]/dt)
        );
        
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi)) == a(force,phi) - a(grad(p),phi)
        );
        
        p_real = p - nu*div(U);
        
        
        //fix pressure to one exact solution node
        p_real += ( psol[0]-p_real[0] );
        p += ( psol[0]-p[0] );
        
        perr = abs(psol-p)/(1.+abs(psol));
        p_realError = abs(psol-p_real)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        
        Case::time().fireWriting();
    }
    
    return 0;
}