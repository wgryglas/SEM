#include <iostream>
#include <fstream>
#include <iterator>

#include "utilities/foreach.h"

#include "iomanagment/Dictionary2.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/case.h"
#include "iomanagment/RegistryFile.h"
#include "components/Basic.h"
#include "mesh/Mesh.h"
#include "fields/GeometricField.h"

int main(int argc, char* argv[])
{
    using namespace std;
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;

    Case::setup("/home/wojtek/Desktop/SEM_test");

    mesh::Mesh mesh;

    GeometricField<Scalar> T(Case::time(),
                             new RegistryFile(Case::time(),"T",Case::time().localPath(),READ_ONCE,AUTO),
                             mesh);

    iomanagment::write(T,cout);

    return 0;
}
