#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fields/DiscontinousField.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "fieldMath/f.h"
#include "fieldMath/postCalculation.h"
#include "fieldMath/weakGradV.h"
#include "fieldMath/GradientOperator.h"
#include "fieldMath/DivergenceOperator.h"
#include "fieldMath/WeakFormDefinition.h"

#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "utilities/VectorUtils.h"
#include "fieldMath/NonlinearOperatorSplitting.h"


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
    
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
    VectorField Ureal("Ureal",mesh,NO_READ);
    VectorField Usol("Usol",mesh,NO_READ);
    VectorField UrealErr("Uerr",mesh,NO_READ);
    VectorField Uerr("UrealErr",mesh,NO_READ);
    
    ScalarField p("p",mesh,NO_READ);
    ScalarField psol("psol",mesh,NO_READ);
    ScalarField perr("perr",mesh,NO_READ);
    ScalarField p_ext("p_ext",mesh,NO_READ, NO_WRITE);
    ScalarField Phi("Phi",mesh,NO_READ);
    ScalarField Divergence("divU", mesh, NO_READ);
    VectorField force("force",mesh,NO_READ);
    
    Phi = 0.;
    Phi |= "fixedGradient";//"fixedValue"
    Phi |= 0.;
    
    
    DiscontinousField<Scalar> p_tmp(p);
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    Scalar nu = Case::material()["nu"];
    auto X = xCmps(mesh.spectralNodes());
    auto Y = yCmps(mesh.spectralNodes());
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    //initial conditions:
    Ureal = U;
    Ureal.registryFile()->fireWriting();
    
    VectorField DU("DU",mesh,NO_READ,NO_WRITE);
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar t = Case::time().time();
        //--setup assumed solution--------------
        
        xCmps(Usol) = -cos(X)*sin(Y)*std::sin(2.*t);
        yCmps(Usol) = sin(X)*cos(Y)*std::sin(2.*t);
        
        psol = -0.25*(cos(2.*X)+cos(2.*Y))*(std::sin(2.*t)*std::sin(2.*t));
        
        xCmps(force) = -cos(X)*sin(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        yCmps(force) = sin(X)*cos(Y)*(2.*std::cos(2*t) + nu*2.*std::sin(2.*t));
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( Phi );
        
        DiscontinousField<Vector> F = force;
        DiscontinousField<Scalar> p_ext = p_tmp;
        DiscontinousField<Vector> velNTerms(mesh);
        for(size_t i=1; i<timeCoeffs.size(); ++i)
        {
            p_tmp -=  (timeCoeffs[i]/timeCoeffs[0])*Phi.oldField(i-1);
            velNTerms += timeCoeffs[i]/dt*oifsConvDerivative(U,dt,i-1);
        }
        
        las::solve
        (
            a(timeCoeffs[0]/dt,U,phi) + a(nu*grad(U),grad(phi) ) == a(p_tmp, grad(phi)) + a(F,phi) - a(velNTerms,phi)
        );
        
        DiscontinousField<Scalar> divU = div(U);
        las::solve
        (
            a(grad(Phi),grad(phi)) == a( (-timeCoeffs[0]/dt)*divU, phi )
        );
        
        
        //calculate real pressure
        p_tmp = p_ext + Phi - nu*divU;
        
        
        //Transform discontinous field into continous result:
        p = p_tmp;
        Divergence = divU;
        Ureal = U -  (dt/timeCoeffs[0])*grad(Phi);
        
        //fix p to match solution in one node, because numerical solution is accurate to constant
        size_t n = p.size()/2;
        p += (psol[n] - p[n]);
        
        //calculate error fields
        perr = abs(psol-p)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        UrealErr = abs(Usol-Ureal)/(1.+abs(Usol));
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}