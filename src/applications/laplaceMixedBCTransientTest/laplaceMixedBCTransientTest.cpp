#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "solver/Solver.h"
#include "solver/CompoundDiscretOperator.h"
#include "solver/DiscretEquation.h"
#include "utilities/VectorUtils.h"
#include "fieldMath/WeakFormDefinition.h"


int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using SEM::field::ScalarField;
    using SEM::mesh::Mesh;
    using SEM::array::abs;
    using SEM::weak::phi;
    using SEM::weak::grad;
    using SEM::field::grad;
    using SEM::weak::a;
    using SEM::weak::bInt;
    using SEM::weak::ddn;
    using SEM::array::abs;
    using SEM::weak::ddt;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    
    //get list of nodes coords
    const numArray<Vector>& nodes=mesh.spectralNodes();
    numArray<Scalar> X( xCmps(nodes) );
    numArray<Scalar> Y( yCmps(nodes) );
    
    
    //create fields
    ScalarField T("T",mesh);
    ScalarField Sol("Sol", mesh, NO_READ);
    ScalarField Err("Err", mesh, NO_READ);
    
    
    numArray<Scalar> solBase = (X*X-1)*(Y*Y-1);
    
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar t = Case::time().time();
        
        Sol = solBase*t;
        las::solve
        (
            a(ddt(T),phi) + a(grad(T),grad(phi)) - bInt(1.,ddn(T),phi)  == a(solBase - 2.*( (X*X-1.) + (Y*Y-1.) )*t, phi)
        );
        
        Err = abs(Sol-T)/(1.+abs(Sol));
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}




