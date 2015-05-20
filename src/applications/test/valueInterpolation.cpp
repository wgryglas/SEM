#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "fieldMath/postCalculation.h"
#include "utilities/VectorUtils.h"

#include "postprocessing/PostObject.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
 
    //Check if some point can be found in element
    std::cout << "elements number: "<<mesh.size()<<std::endl;
    
    Vector p(0.25,0.75);
//     Vector p(0.25,0.);
    
    for(size_t i=0; i<mesh.size(); ++i)
    {  
        std::cout<<"does element contain point=("<<p[0]<<", "<<p[1]<<")?"<<mesh[i].containsPoint(p)<<std::endl;
        if(mesh[i].containsPoint(p))        
        {
            Vector l = mesh[i].spectralElement().mapToLocal(p,mesh[i].nodes());
            std::cout<<"local coords for this point are: ("<<l[0]<<", "<<l[1]<<")"<<std::endl;
        }
    }
    
    iomanagment::Dictionary dict("valueOnLine");
    
    dict.add( new DictEntry("writeType","columns") );
    dict.add( new DictEntry("field","T") );
    dict.add( new DictEntry("fileName","lineData.dat") );
    dict.add( new DictEntry("start","(0.5 0)") );
    dict.add( new DictEntry("end","(0.5 1)") );
    dict.add( new DictEntry("resolution","50") );
    
    
//     ScalarField T("T",mesh, iomanagment::NO_READ);
//     T = SEM::array::sin(yCmps(mesh.spectralNodes()) );
//     T.registryFile()->fireWriting();
    
    Case::time().setCurrentTime(0.);
    
    PostObject* valuesOnLine = PostObject::Impl("scalarValuesOnLine")( dict, mesh, Case::time() );
    
    valuesOnLine->calculate();
    
}