#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/ddt.h"
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
    using SEM::field::grad;
    using SEM::weak::a;
    using SEM::weak::bInt;
    using SEM::weak::ddn;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    ScalarField T("T", mesh, NO_READ);
    
    T |="fixedGradient";
    T |=0.;
    
    ScalarField Sol("Sol", mesh, NO_READ);
    ScalarField Err("Err", mesh, NO_READ);
    
    
    //get list of nodes coords
    const numArray<Vector>& nodes=mesh.spectralNodes();
    
    numArray<Scalar> X( xCmps(nodes) );
    numArray<Scalar> Y( yCmps(nodes) );
    
    Scalar Pi = 4.*std::atan(1.);
    
    using SEM::array::sin;
    using SEM::array::cos;
    using SEM::array::abs;
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();

        Sol = cos(Pi*X)*cos(Pi*Y);
        
        las::solve
        (
            a(grad(T),grad(phi)) -bInt(1.,ddn(T),phi)== a(2.*Pi*Pi*cos(Pi*X)*cos(Pi*Y), phi)//, SEM_ConjugateGradient
        );
        
        //fix solution in one point due to possible diffrence in constant value (equation matrix is singular)
        T = T+(Sol[0]-T[0]);
        
        Err = abs(T-Sol)/(1.+abs(Sol));
        
        Case::time().fireWriting();
        
    }
    
    return 0;
    
}




