#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/ScalarField.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"


int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    
    //Setup case
    Case::setup("/home/wojtek/Desktop/SMESH_test3");
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    ScalarField T
    (
        new RegistryFile
        (
            Case::time(),
         "T",
         Case::time().localPath(),
         iomanagment::READ_ONCE,
         iomanagment::AUTO
        ),
     mesh
    );
    
    //get material coeff
    Scalar coeff = Case::material()["kt"]/(Case::material()["rho"]*Case::material()["Cp"]);
    
    //solve time loop
    while(!Case::time().end())
    {
        las::solve(ddt(T)-laplacian(coeff,T));
        ++Case::time();
    }
    
    return 0;
    
}
