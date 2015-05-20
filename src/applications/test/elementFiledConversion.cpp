
#include <iostream>

#include "iomanagment/case.h"
#include "utilities/numArray.h"
#include "utilities/ArrayFunctions.h"

#include "fields/ElementFieldBase.h"

#include "fields/GeometricField.h"

#include "fields/DiscontinousField.h"

int main(int argc, char* argv[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    
    Case::setup(argc,argv);
    
    //Create mesh
    mesh::Mesh mesh(Case::meshPath(),Case::elementsPath());
    
    numArray<Scalar> contField(mesh.nodesNumber(),2.);
    
    iomanagment::write(contField,std::cout);
    
    SEM::field::GeometricField<Scalar> cf(mesh);
    
    field::DiscontinousField<Scalar> dF(mesh);
    
}