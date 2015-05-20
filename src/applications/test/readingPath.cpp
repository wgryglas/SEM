
#include <iostream>
#include <vector>

#include "boost/filesystem.hpp"
#include "iomanagment/Dictionary2.h"

int main(int argc, char* args[])
{
    using namespace SEM::iomanagment;
//     boost::filesystem::path path="/home/wojtek/.smesh/smesh.config";
//     
//     
//     Dictionary dict;
//     path>>dict;
//     
//     std::cout<<dict<<std::endl;
    
    boost::filesystem::path pathString="/home/wojtek/.smesh/solvers/laplace.sem";
    
    Dictionary dictString;
    
    pathString>>dictString;
    
    std::cout<<dictString<<std::endl;
  
    std::vector<std::string> strings;
    dictString.subDictionary("fields").subDictionary("T").entry("boundaryConditions")>>strings;
    
    for(int i=0; i<strings.size(); ++i)
        std::cout<<strings[i]<<std::endl;
 
    return 0;
}