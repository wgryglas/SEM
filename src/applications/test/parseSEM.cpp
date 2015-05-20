
#include <iostream>
#include <string>
#include <fstream>

#include "boost/filesystem.hpp"

#include "iomanagment/case.h"

int main(int argc, char* argv[])
{
  using namespace std;
  using SEM::Case;

  Case::setup("/home/wojtek/Desktop/SEM_test");

    cout<<Case::path()<<endl;

    for(int i=0; i<5; ++i,++Case::time())
         cout<<Case::time()<<endl;

  return 0;
}
