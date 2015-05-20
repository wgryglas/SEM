#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"


#include "utilities/foreach.h"
#include "utilities/ListWrap.h"

#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"
#include "iomanagment/RegistryFile.h"
#include "iomanagment/case.h"

#include "components/Basic.h"

#include "fields/GenericField.h"

#include "time/Time.h"

int main(int argc, char* argv[])
{
  using namespace std;
  using namespace SEM;
  using namespace SEM::iomanagment;
  using namespace SEM::field;

  Case::setup("/home/wojtek/Desktop/SEM_test");

  Time runTime;

  RegistryFile::ref file(new RegistryFile( runTime,"testField",runTime.localPath(), iomanagment::READ_ONCE, iomanagment::AUTO) );

  cout<<file->filePath()<<endl;

  GenericField<Vector> v(file);

  GenericField<Vector> v2(v);

  v+=v2;
  cout<<runTime<<endl;
  iomanagment::write(v,cout);
  ++runTime;

  v+=v2;
  cout<<runTime<<endl;
  iomanagment::write(v,cout);
  ++runTime;

  return 0;
}
