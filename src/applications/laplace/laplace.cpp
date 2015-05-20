#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
// #include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "fieldMath/postCalculation.h"

#include "solver/CompoundDiscretOperator.h"

#include "fieldMath/GradientOperator.h"
#include "fieldMath/WeakFormDefinition.h"
#include "solver/FixedOperator.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using SEM::weak::phi;
    using SEM::weak::ddt;
    using SEM::field::grad;
    using SEM::weak::a;
    using SEM::weak::bInt;
    using SEM::weak::ddn;
    

    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    
    //Create T field
    ScalarField T("T",mesh);
    
    //Create post grad T field
    VectorField gradT("gradT", mesh, NO_READ);
    
    
    //get material coeff
    Scalar coeff = Case::material()["kt"]/(Case::material()["rho"]*Case::material()["Cp"]);
    
    
    //auto equation = las::fixOperator( a(ddt(T),phi) + a(coeff*grad(T), grad(phi)) );
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        las::solve
        (
            a(ddt(T),phi) + a(coeff*grad(T),grad(phi)) -bInt(coeff,ddn(T),phi)
//             equation - bInt(coeff,ddn(T),phi)
        );
        
        gradT = -grad(T);
        
        
        
     
        Case::time().fireWriting();
    }
    
    return 0;
    
}




