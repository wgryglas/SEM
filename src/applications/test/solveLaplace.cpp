#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"

#include "utilities/foreach.h"
#include "utilities/TypeDefs.h"
#include "utilities/numArray.h"
#include "utilities/numArray2D.h"
#include "utilities/ArrayFunctions.h"


#include "components/Basic.h"

#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"
#include "iomanagment/case.h"

#include "mesh/Mesh.h"

#include "fields/GeometricField.h"

#include "solver/DiscretOperator.h"
#include "solver/EquationPart.h"
#include "solver/Solver.h"

#include "fieldMath/laplacian.h"
#include "fieldMath/ddt.h"
#include "fieldMath/f.h"

#include "materials/Property.h"

int main(int argc, char* argv[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;

    //Setup case paths
//     Case::setup("/home/wojtek/Desktop/SEM_test",//WD
//                 "/home/wojtek/Dropbox/PRACA_MAGISTERSKA/SEM_PROJECT/src");//INSTALL DIR
    Case::setup("/home/wojtek/Desktop/SMESH_test3");

    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());

    //Create T field
    GeometricField<Scalar> T(Case::time(),
                             new RegistryFile(
                                            Case::time(),
                                            "T",
                                            Case::time().localPath(),
                                            READ_ONCE,
                                            AUTO),
                            mesh
                            );


    //get material coeff
    Scalar coeff = Case::material()["DT"];
    
    //solve time loop
    while(!Case::time().end())
    {
        las::solve(ddt(T)-laplacian(coeff,T));
        ++Case::time();
    }

    return 0;
}


