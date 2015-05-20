#include <iostream>

#include "iomanagment/case.h"
#include "utilities/numArray.h"
#include "utilities/ArrayFunctions.h"

#include "fields/ElementFieldBase.h"
#include "fields/GeometricField.h"

#include "utilities/RK4Integrator.h"

int main(int argc, char* argv[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    
//     Case::setup(argc,argv);
    //Create mesh
//     mesh::Mesh mesh(Case::meshPath(),Case::elementsPath());

    SEM::RK4Integrator<double> integrator
    (
        [](double time, double y)->double
        {
            return std::cos(time);
        }
    );
    
    std::vector<double> result(1,0.);
    
    double start=0.;
    double end=3.;
    size_t steps=100;
    double step = (end-start)/(steps-1);
    integrator.integrate(start,end, steps, result);
    
    
    for(size_t t=0; t<result.size(); ++t)
    {
        double var = start+t*step;
        std::cout << "diff = "<<result[t]-std::sin(var)<<",num="<<result[t]<<", exact="<<std::sin(var)<<std::endl;
    }
    
    double lastRes = integrator.integrate(start,step,steps-1,0.);
    std::cout<<std::endl<<std::endl<< "diff as sigle result="<<lastRes-std::sin(end)<<", num="<<lastRes<<", exact="<<std::sin(end)<<std::endl;
    
    
    
    
    
}