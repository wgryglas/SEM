
#include <iostream>
#include <vector>

#include "boost/filesystem.hpp"
#include "iomanagment/Dictionary2.h"
#include "mesh/Mesh.h"
#include "iomanagment/case.h"

int main(int argc, char* args[])
{
    using namespace SEM::iomanagment;
    using namespace SEM;
    
//     boost::filesystem::path meshFile="/home/wojtek/Desktop/SMESH_test3/mesh/mesh.sem";
//     boost::filesystem::path elementsFile="/home/wojtek/Desktop/SMESH_test3/mesh/elements.sem";
        
    boost::filesystem::path casePath="/home/wojtek/Desktop/SMESH_test3";
    
    Case::setup(casePath);
    
    mesh::Mesh mesh(Case::meshPath(),Case::elementsPath());
    
    std::cout<<mesh<<std::endl;
    
    
    return 0;
}