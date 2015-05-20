#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"
#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"

#include "components/Basic.h"
#include "components/CmpTraits.h"

#include "mesh/Mesh.h"
#include "fields/GeometricField.h"

#include "utilities/TypeDefs.h"
#include "utilities/numArray.h"
#include "utilities/numArray2D.h"
#include "utilities/VectorUtils.h"

#include "iomanagment/case.h"


#define PRINT(EXPR) std::cout<<#EXPR<<'='<<EXPR<<std::endl

int main(int argc, char* argv[])
{
    using namespace SEM;


// //     numArray<Scalar> a1(5,1.);
//     Scalar pi = 4.*std::atan(1.);
// //     a1 = pi*a1 + pi*a1;
// //     
// //     
// //     iomanagment::write(a1,std::cout<<"2*pi*5[1]"<<std::endl)<<std::endl;
//     
//     
//     Case::setup(argc, argv);
//     
//     mesh::Mesh mesh(Case::meshPath(),Case::elementsPath());
// 
//     field::GeometricField<Scalar> f1(mesh);
//     
//     auto x = xCmps(mesh.spectralNodes());
//     auto y = yCmps(mesh.spectralNodes());
//     
//     f1 = SEM::array::sin(x)*SEM::array::sin(y);
//     
//     field::GeometricField<Scalar> f2(mesh);
//     
//     f2 = pi*pi*f1;
//     
//     numArray<Scalar> result = (f2-f1)/(f1*(pi*pi-1));
//     
//     iomanagment::write(result,std::cout<<"(f2-f1)/(pi*pi-1)="<<std::endl)<<std::endl;
//     
//     iomanagment::write(mesh.spectralNodes(),std::cout<<"coords="<<std::endl)<<std::endl;
//     
    
//     numArray<numArray<Scalar> > list2d(2);
//     list2d[0].resize(2,0.);
//     list2d[1].resize(2,0.);
//     
//     list2d=2.;
//     list2d[1]=1.;
//     
//     iomanagment::write(list2d, std::cout<<"list 2d(2x2, [0][:]=2, [1][:]=1"<<std::endl)<<std::endl;
    
    numArray2D<Scalar> a2d(2,2,0.);
    iomanagment::write(a2d,std::cout<<"array2d"<<std::endl)<<std::endl;
    
    
    std::stringstream ss("[[4]]");
    
    iomanagment::read(ss,a2d);
    
    iomanagment::write(a2d,std::cout<<"array2d"<<std::endl)<<std::endl;
    
    
    return 0;
}
