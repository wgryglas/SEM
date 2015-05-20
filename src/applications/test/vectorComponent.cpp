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

#include "fields/GenericField.h"

#include "utilities/TypeDefs.h"
#include "utilities/numArray.h"
#include "utilities/numArray2D.h"
#include "utilities/VectorUtils.h"

int main(int argc, char* argv[])
{
    using namespace SEM;
    
    numArray<Vector> vs(10);
    
    for(int i=0; i<10; ++i)
    {
        vs[i].x() = i;
        vs[i].y() = 2*i;
    }
    
    SEM::iomanagment::write(vs, std::cout<<"initial vectors list"<<std::endl)<<std::endl;
    
    numArray<Vector> vs2(vs);
    
    xCmps(vs) +=xCmps(vs) + yCmps(vs);
    
    SEM::iomanagment::write(vs, std::cout<<" x-es = x-es + y-es"<<std::endl)<<std::endl;

    numArray<double> x_part = xCmps(vs).slice(2,4);
    
    SEM::iomanagment::write(x_part, std::cout << "list slice of x-cmps(2-4) = "<<std::endl)<<std::endl;
    
    numArray<double> comp_mult = xCmps(vs) * yCmps(vs);
   
    SEM::iomanagment::write(comp_mult, std::cout<<"components multiplication of each vector in current list"<<std::endl)<<std::endl;
    
    
    SEM::iomanagment::write(vs2,std::cout<<"initial list"<<std::endl)<<std::endl;
    
    SEM::iomanagment::write(numArray<Scalar>(array::dotProd(vs2,vs2)),std::cout<<"dot product of initial list with initial list"<<std::endl)<<std::endl;
    
    SEM::iomanagment::write(numArray<Scalar>(array::crossProd(vs2,vs2)),std::cout<<"cross product of initial list with initial list"<<std::endl)<<std::endl;
    
    SEM::iomanagment::write(numArray<Scalar>(array::mag(vs2)), std::cout<<"magnintued of each vector in initial list"<<std::endl) << std::endl;
    
    numArray<Vector> vs3(vs2);
    
//     CmpTraits<Vector>::const_CmpArray cmpArray  = CmpTraits<Vector>::cmpArray(vs2,0);
//     numArray<
    CmpTraits<Vector>::cmpArray(vs3,0).slice(2,4) -= CmpTraits<Vector>::cmpArray(vs2,0).slice(2,4)*2.;
    
    SEM::iomanagment::write(vs3, std::cout<<"plot vector list with x component (2:4) of initail list multiplied by 2"<<std::endl) << std::endl;
    
    numArray<Vector> vs4 = vs3 + vs2;
    SEM::iomanagment::write(vs4, std::cout<<"vs4 = sum of vs3 and vs2"<<std::endl) << std::endl;
    
    vs4 = -1.*vs4;
    SEM::iomanagment::write(vs4, std::cout<<"vs4 =-vs4"<<std::endl) << std::endl;
    
    return 0;
}
