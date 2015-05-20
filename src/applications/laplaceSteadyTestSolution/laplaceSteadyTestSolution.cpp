#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "solver/Solver.h"
#include "solver/CompoundDiscretOperator.h"
#include "solver/DiscretEquation.h"
#include "utilities/VectorUtils.h"
#include "fieldMath/WeakFormDefinition.h"


void solveEq(SEM::field::GeometricField<SEM::Scalar> &f)
{
    using namespace SEM;
    const mesh::Mesh & m = f.mesh();
    
    las::SEMMatrix matrix(1);
    las::SEMVector vector(1);
    vector[0].resize(f.size());
    
    const numArray<Vector> nodes=m.spectralNodes();
    
    
    matrix[0].resize(m.size());
    
    for(size_t e=0; e<m.size(); ++e)
    {
        numArray2D<Scalar>::const_mappedArray massVector(m[e].massMatrix(), m[e].localNodesInMatrix());
        const std::vector<int> &mask=m[e].indexVectorMask();
        
        for(size_t n=0; n<mask.size(); ++n)
        {
            Scalar x = nodes[mask[n]].x();
            Scalar y = nodes[mask[n]].y();
            vector[0][mask[n]] = massVector[n]*(-2.0*( (6.0*y-3.0)*x*(x-1.0) + y*y*(2.0*y-3.0)));
        }
            
        size_t nodesNum =m[e].indexVectorMask().size();
        matrix[0][e].resize2D(nodesNum,nodesNum);
        matrix[0][e] = m[e].stiffMatrix(1.);
    }
    
    las::solve(f,matrix,vector);
}

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using SEM::field::ScalarField;
    using SEM::mesh::Mesh;
    
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
    
    using SEM::array::abs;
    using SEM::weak::phi;
    using SEM::weak::grad;
    using SEM::field::grad;
    using SEM::weak::a;
    using SEM::weak::bInt;
    using SEM::weak::ddn;
    
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        las::solve
        (
            a(grad(T),grad(phi)) == a(-2.0*( (6.0*Y-3.0)*X*(X-1.0) + Y*Y*(2.0*Y-3.0)), phi)
        );
        
        Sol = Y*Y*(2.0*Y-3.0)*X*(X-1.0) + X;
        Err = abs(T-Sol)/(1.0+abs(T));
        
        Case::time().fireWriting();
    }
    
    return 0;
    
}





