#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
// #include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/VectorUtils.h"
#include "fieldMath/f.h"
#include "fieldMath/WeakFormDefinition.h"


int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using namespace SEM::array;
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

    ScalarField T("T",mesh);
    ScalarField Sol("Sol",mesh,NO_READ);
    ScalarField Err("Err",mesh,NO_READ);

    //get list of nodes coords
    const numArray<Vector>& nodes=mesh.spectralNodes();

    auto X = xCmps(nodes);
    auto Y = yCmps(nodes);
    
    Scalar pi = 4*std::atan(1.);

    using SEM::array::abs;
    using SEM::array::sin;

    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Sol = sin(pi*X)*sin(pi*Y);
        
        las::solve
        (
            a(grad(T),grad(phi)) == a(2.*pi*pi*Sol,phi)
        );

        Err = abs(T-Sol)/(1.0+abs(Sol));
        
        Case::time().fireWriting();
    }

    return 0;

}




