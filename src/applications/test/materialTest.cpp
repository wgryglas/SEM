#include <iostream>

#include "iomanagment/case.h"
#include "materials/Material.h"


int main(int argc, char* argv[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    
    
    //with out env setup, direct pass
    //Case::setup("/home/wojtek/Desktop/SEM_test","/home/wojtek/Dropbox/PRACA_MAGISTERSKA/SEM_PROJECT/export");
    
    //env setup (or not then runtime error)
    Case::setup("/home/wojtek/Desktop/SEM_test");
    
    materials::Material & material=Case::material();
    
    std::cout<<"value fo DT="<<material["DT"]<<std::endl;
    
    
    return 0;
}