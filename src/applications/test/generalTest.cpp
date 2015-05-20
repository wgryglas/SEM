#include <iostream>
#include <fstream>
#include <iterator>

#include "iomanagment/Dictionary2.h"

#include "fields/GeometricField.h"

#include "elements/SpectralElement.h"
#include "elements/SE_QLGL.h"
#include "elements/SE_TLGL.h"

#include "utilities/ImplementationFactory.h"

#include "boost/lexical_cast.hpp"

struct T1 
{
    DECLARE_IMPLEMENTATION_FACTORY(T1,unsigned int)
public:
    virtual std::string kurwa() =0;
};

DEFINE_IMPLEMENTATION_FACTORY(T1,unsigned int)

struct T2 : public T1
{
    unsigned int _a;
    T2(unsigned int a):_a(a){}
    std::string kurwa() {return "to chuj "+boost::lexical_cast<std::string>(_a);}
};

REGISTER_IMPLEMENATION(T1,T2,"T2");

struct T3 : public T1
{
    unsigned int _a;
    T3(unsigned int a):_a(a){}
    std::string kurwa() {return "to pizda szmata "+boost::lexical_cast<std::string>(_a);}
};

REGISTER_IMPLEMENATION(T1,T3,"T3");




int main(int argc, char* argv[])
{
  using namespace std;

  using namespace SEM;
  using namespace SEM::field;


//  boost::array<int,3> a;

//  stringstream ss("(1 2 3)");
//  iomanagment::read<boost::array<int,3> >(ss, a);
//  iomanagment::write<boost::array<int,3> >(a,std::cout);
//  cout<<endl;

//  stringstream ss2("( (1 2 3) (4 5 6) )");
//  typedef boost::array<boost::array<int,3>,2> mat;
//  mat b;
//  iomanagment::read<mat>(ss2, b);
//  iomanagment::write<mat>(b,std::cout);
//  cout<<endl;

//   stringstream ss3("((1 2)(3 4))");
// 
//   typedef GenericField<Vector> listVec;
//   listVec v(2);
//   iomanagment::read<listVec>(ss3,v);
// 
//   iomanagment::write<listVec>(v,cout);

    mesh::SpectralElement* se = mesh::SpectralElement::Impl("SE_QLGL")(5);
  
  
  T1* c = T1::Impl("T2")(5);
  std::cout << c->kurwa() << std::endl;
  
  return 0;
}
