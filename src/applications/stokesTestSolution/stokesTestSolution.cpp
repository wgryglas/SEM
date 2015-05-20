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
    VectorField Uerr("Uerr",mesh,NO_READ);
    VectorField UrealError("UrealError",mesh,NO_READ);
    
    //Create p field
    ScalarField p("p",mesh,NO_READ);
    ScalarField psol("psol",mesh,NO_READ);
    ScalarField perr("perr",mesh,NO_READ);
    
    //Create p* - pressure extrapolation field
    ScalarField p_ext("pext",mesh,NO_READ,NO_WRITE);
    
    //Create Phi temporal field 
    ScalarField Phi("Phi",mesh,NO_READ);
    ScalarField divergence("divU",mesh,NO_READ);
    
    
    Phi = 0.;
    Phi |= "fixedGradient";//"fixedValue"
    Phi |= 0.;
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    numArrayVectorComponent<const numArrayBase<numArray<Vector> > > X=xCmps(mesh.spectralNodes());
    numArrayVectorComponent<const numArrayBase<numArray<Vector> > > Y=yCmps(mesh.spectralNodes());

    numArray<Vector> force(X.size());

    Scalar Pi = 4.*atan(1.);
    Scalar nu = 1.;
    
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
    
    Ureal = U;
    Ureal.registryFile()->fireWriting();
    
    
    p =-cx*sy; 
    p |= "fixedGradient";
    p |= 0.;
    p.registryFile()->fireWriting();
    
    Usol = U;
    Usol.registryFile()->fireWriting();
    
    psol = p;
    psol.registryFile()->fireWriting();
    
    
    DiscontinousField<Scalar> p_tmp(p);
    
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
        
        DiscontinousField<Scalar> p_ext = p_tmp;
//         p_ext = p;//2.*p.cachedField(0)-p.cachedField(1);//p
//         p_tmp = p_ext;
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( Phi );
        for(size_t i=1; i<timeCoeffs.size(); ++i)
        {
            p_tmp -= (timeCoeffs[i]/timeCoeffs[0])*Phi.oldField(i-1);
        }
        
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi) ) == a(force,phi) + a(p_tmp, grad(phi))
        );
        
        DiscontinousField<Scalar> divU = div(U);
        las::solve
        (
            a(grad(Phi),grad(phi)) == a( (-timeCoeffs[0]/dt)*div(U), phi )
        );
        
        
        //calculate real pressure
        p_tmp = p_ext + Phi - nu*divU;
        
        //Transform discontinous field into continous result:
        p = p_tmp;
        Ureal = U -  (dt/timeCoeffs[0])*grad(Phi);
        divergence = divU;
        
        
        //Calculate error
        size_t pos = p.size()/2;
        p += psol[pos] - p[pos];
        perr = abs(psol-p)/(1.+abs(psol));
        Uerr = abs(Usol-U)/(1.+abs(Usol));
        UrealError= abs(Usol-Ureal)/(1.+abs(Usol));
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}