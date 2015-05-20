#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"

#include "utilities/foreach.h"
#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"

#include "components/Basic.h"
#include "fields/GenericField.h"

int main(int argc, char* argv[])
{
    using namespace std;
    using namespace SEM;
    using namespace SEM::field;

    stringstream ss("((1 2)(3 4))");
    GenericField<Vector> field;
    iomanagment::read(ss,field);

    stringstream ss2("((1 2)(3 4))");
    GenericField<Vector> field2;
    iomanagment::read(ss2,field2);
    iomanagment::write(field2,cout<<"field2"<<endl)<<endl;

    field2+=field;
    GenericField<Vector> field3(field.size());
    field3=field+field2;

    iomanagment::write(field3,cout<<"field3"<<endl)<<endl;

    GenericField<Scalar> field4;
    stringstream ss3("(1 2 3)");
    iomanagment::read(ss3,field4);
    iomanagment::write(field4,cout<<"field4"<<endl)<<endl;


    iomanagment::Dictionary d("pole");
    field.write(d);
    cout<<"dictionary pole"<<endl<<d<<endl;

    //uniform field read
    stringstream unS("[(1 2)]");
    GenericField<Vector> v7;
    v7.resize(3);

    iomanagment::read(unS,v7);

    iomanagment::write(v7, cout<<"v7"<<endl)<<endl;

    cout<<endl<<"field7 elemts:"<<endl;
    for(int i=0; i<v7.size(); ++i)
        iomanagment::write(v7[i],cout<<"v7[i]=")<<endl;

  return 0;
}

/////////////////////////////////
// Dictionary read/write
/////////////////////////////////
//    SEM::iomanagment::Dictionary dict("control");
//    "/home/wojtek/Desktop/control.sem">>dict;
//    cout<<dict;

/////////////////////////////////
// String read/write
/////////////////////////////////
//    stringstream ss;
//    string in="cos tam";
//    SEM::iomanagment::write<string>(in,ss);
//    string out;
//    SEM::iomanagment::read<string>(ss,out);
//    cout<<out<<endl;

/////////////////////////////////
// List read/write
/////////////////////////////////
//    vector<int> inVec2;
//    inVec2.push_back(0);
//    inVec2.push_back(1);
//    inVec2.push_back(2);
//    inVec2.push_back(3);
//
//    vector<int> inVec3;
//    inVec3.push_back(0);
//    inVec3.push_back(1);
//    inVec3.push_back(2);
//
//    vector<vector<int> > inVec;
//    inVec.push_back(inVec2);
//    inVec.push_back(inVec3);
//
//    stringstream ss;
//    SEM::iomanagment::write<vector<vector<int> > >(inVec,ss);
//
//    std::cout<<ss.str()<<endl;
//
//    vector<vector<int> > outVec;
//    SEM::iomanagment::read<vector<vector<int> > >(ss, outVec);
//
//    cout<<"out vecotr::"<<endl;
//    foreach(vector<int> v, outVec)
//    {
//        cout<<"---------"<<endl;
//        foreach(int i,v)
//        {
//            cout<<i<<endl;
//        }
//    }

/////////////////////////////////
// List read
/////////////////////////////////
//    string s="(3 4 5 6 7)";
//    stringstream ss(s);
//
//    vector<int> v;
//    SEM::iomanagment::read<vector<int> >(ss,v);
//
//    stringstream ss2;
//    SEM::iomanagment::write<vector<int> >(v,ss2);
//
//    SEM::iomanagment::Indenter().indetStream<ostream>(ss2,cout);

///////////////////////////////////
// Vector derived list read
//////////////////////////////////
//    using namespace SEM::field;
//    using namespace SEM::iomanagment;
//
//    string s="(3 4 5 6 7)";
//    stringstream ss(s);
//
//    GenericField<int> f;
//    ss>>f;
//    cout<<f;

