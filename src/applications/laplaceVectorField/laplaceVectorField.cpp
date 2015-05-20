#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
// #include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "solver/CompoundDiscretOperator.h"

#include "fieldMath/WeakFormDefinition.h"

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
    
    //Create U field
    VectorField U("U",mesh);
    
    //get material coeff
    Scalar nu = Case::material()["nu"];
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        las::solve
        (
          a(weak::ddt(U),phi) + a(grad(U),grad(phi)) - bInt(1.,ddn(U),phi)
        );
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}




