#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "fieldMath/postCalculation.h"
#include "utilities/VectorUtils.h"
#include "fieldMath/f.h"

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
    
    //Create T field
    ScalarField T("T", mesh);
    ScalarField Sol("Sol", mesh, NO_READ);
    ScalarField Err("Err", mesh, NO_READ);
    
    numArray<Scalar> X( xCmps( mesh.spectralNodes() ) );
    numArray<Scalar> Y( yCmps( mesh.spectralNodes() ) );
    
    Scalar Pi = 4.*std::atan(1.);
    numArray<Scalar> sx = SEM::array::sin(Pi*X);
    numArray<Scalar> sy = SEM::array::sin(Pi*Y);

    T = sx*sy;
    T|="fixedValue";
    T|=0.;
    T.registryFile()->fireWriting();
    
    Sol=T;
    Sol.registryFile()->fireWriting();
    
    Err=0;
    Err.registryFile()->fireWriting();
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        Scalar st = std::sin(Case::time().time());
        Scalar ct = std::cos(Case::time().time());
        
        Sol = ct*sx*sy;
        
        las::solve(
                    a(ddt(T),phi) + a(grad(T),grad(phi)) - bInt(1.,ddn(T),phi)  == a(sx*sy*(-st + 2.*Pi*Pi*ct), phi)
                  );
        
        Err= SEM::array::abs(Sol-T)/(1.+SEM::array::abs(Sol));
        
        Case::time().fireWriting();
        
    }
    
    return 0;
}



