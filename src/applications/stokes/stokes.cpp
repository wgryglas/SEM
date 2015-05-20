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
    ScalarField p("p",mesh,READ_IF_PRESENT);
    ScalarField p_ext("p_ext",mesh,NO_READ, NO_WRITE);
    ScalarField Phi("Phi",mesh,NO_READ);
    
    
    Phi = 0.;
    Phi |= "fixedGradient";//"fixedValue"
    Phi |= 0.;
    
    
    DiscontinousField<Scalar> p_tmp(p);
    
    //--------------------- GET COEFFS -----------------------------------
    std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( Case::solutionControl().timeDiscretization() );
    Scalar dt = Case::time().timeStep();
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    //initial conditions:
    Ureal = U;
    Ureal.registryFile()->fireWriting();
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        DiscontinousField<Scalar> p_ext = p_tmp;
        
        for(size_t i=1; i<timeCoeffs.size(); ++i)
        {
            p_tmp -=  (timeCoeffs[i]/timeCoeffs[0])*Phi.oldField(i-1);
        }
        
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi) ) == a(p_tmp, grad(phi))
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
        
        Ureal = U -  (dt/timeCoeffs[0])*grad(Phi);
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}