
#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"

#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"

#include "components/Basic.h"

//#include "fields/ContinousField.h"

#include "utilities/TypeDefs.h"
#include "utilities/numArray.h"
#include "utilities/numArray2D.h"

#include "solver/EquationAssignment.h"

int main(int argc, char* argv[])
{
  using namespace std;

  using namespace SEM;
  using namespace SEM::iomanagment;
//   using namespace SEM::field;

   numArray<double> d;
   stringstream ss("(1 2 3)");
   iomanagment::read(ss,d);
   iomanagment::write(d,cout<<endInfo<<"d{read(\"1 2 3\",d)}="<<endl);


   numArray<double> d3(3);
   iomanagment::write(d3,cout<<endInfo<<"d3{d3(3)}="<<endl)<<endl;

   las::AddEqAssigment assign;
   assign(d + sin(d) + 10.0,d3);
   
   iomanagment::write(d3,cout<<endInfo<<"d3{d3=d+sin(d)+10.0}="<<endl)<<endl;
   
   numArray2D<Scalar> a1(3,3,0.);
   numArray2D<Scalar> a2(3,3,1.);
//    assign(a2,a1);
   a1+=a2;
   
   std::cout<<endInfo<<"a2.size()="<<a2.size()<<std::endl;
   std::cout<<endInfo<<"a2[1].size()="<<a2[1].size()<<std::endl;
   std::cout<<endInfo<<"a2[1][1]="<<a2.evalAt(5)<<std::endl;
   
   iomanagment::write(a1,cout<<endInfo<<"a1(3x3)=1."<<std::endl);
   
   
   iomanagment::write(d3,cout<<endInfo<<"d3{d3 += d + sin(d) + 10.0}="<<endl);

   d3=4*atan(1);
   iomanagment::write(d3,cout<<endInfo<<"d3{d3=4*atan(1)}="<<endl);

//    ContinousField<double> genF(d3);
//    iomanagment::write<ContinousField<double> >(genF,cout<<endInfo<<"genF{genF(d3)}="<<endl);
//    genF+=d3+cos(d3);
//    iomanagment::write<ContinousField<double> >(genF,cout<<endInfo<<"genF{genF+=d3+cos(d3)}="<<endl);
// 
//    numArray<double> v(5);
//    v = 2;
//    iomanagment::write(v,cout<<endInfo<<"v{v=2} ="<<endl);
// 
//    numArray<double> v2(v);
//    iomanagment::write(v2,cout<<endInfo<<"v{v2(v)} ="<<endl);
// 
//    v.slice(1,4)=sin( v.slice(1,4) ) + v2.slice(1,4);
//    iomanagment::write(v,cout<<endl<<"v{v.slice(1,4)=sin( v.slice(1,4) ) + v2.slice(1,4)} ="<<endl);
// 
// 
//    numArray2D<double> v2d(2,5);
//    iomanagment::write(v2d,cout<<endl<<"v2d{v2d(2,5)}=");
//    v2d[0]=1;
//    v2d[1]=2;
//    iomanagment::write(v2d,cout<<endl<<"v2d{v2d[0]=1, v2d[1]=2}=");
// 
//    numArray2D<double> v2d2(v2d);
//    v2d2[0]=3;
//    v2d2[1]=4;
//    iomanagment::write(v2d2,cout<<endl<<"v2d2{v2d2[0]=3, v2d[1]=4}=")<<endl;
// 
//    
//    v2d.slice(0,1,1,3) = -1.;
//    iomanagment::write(v2d,cout<<endl<<"v2d{v2d.slice(0,1,1,3) = -1.}=")<<endl;
//    
//    v2d.slice(0,1,1,3) = v2d.slice(0,1,1,3) + v2d2.slice(0,1,1,3);
//    iomanagment::write(v2d,cout<<endl<<"v2d{v2d.slice(0,1,1,3) = v2d + v2d2}=")<<endl;
// 
//    MatrixToMatrixMap map;
//    VectorToMatrixMap vMap(3);
//    boost::array<int,2> indexes;
// 
//    indexes[0]=0;
//    indexes[1]=1;
//    vMap[0]=indexes;
//    indexes[1]=2;
//    vMap[1]=indexes;
//    indexes[1]=3;
//    vMap[2]=indexes;
//    map.push_back(vMap);
// 
//    indexes[0]=1;
//    indexes[1]=1;
//    vMap[0]=indexes;
//    indexes[1]=2;
//    vMap[1]=indexes;
//    indexes[1]=3;
//    vMap[2]=indexes;
//    map.push_back(vMap);
// 
//    v2d[0]=1;
//    v2d[1]=2;
//    iomanagment::write(v2d,cout<<endl<<"v2d{v2d[0]=1, v2d[1]=2}=");
// 
//    v2d.slice(map) = v2d.slice(map) + v2d2.slice(map);
//    iomanagment::write(v2d,cout<<endl<<"v2d{map{0-1,1-3} = v2d + v2d2}=")<<endl;
// 
//    MatrixToVectorMap mTVMapping(2);
//    std::vector<int> ind(5);
//    ind[0]=0;ind[1]=1;ind[2]=2;ind[3]=3;ind[4]=4;
//    mTVMapping[0]=ind;
//    ind[0]=5;ind[1]=6;ind[2]=7;ind[3]=8;ind[4]=9;
//    mTVMapping[1]=ind;
// 
//    numArray<double> toMap(10);
//    toMap = 3;
//    toMap.slice(5,9) = 2;
//    iomanagment::write(toMap, cout<<"toMap{1-4=3,5-9=2}=")<<endl;
// 
//    numArray2D<double> mappedVecToArray2D(toMap.sliceArray2D(mTVMapping));
//    iomanagment::write(mappedVecToArray2D,cout<<"mappedVec{=toMap.sliceArray2D(mTVMapping[0]=1-4, mTVMapping[1]=5-9})}=")<<endl;
// 
// 
//    VectorToMatrixMap vToMMap(3);
//    boost::array<int,2> mIndex;
//    mIndex[0]=0;
//    mIndex[1]=1;
//    vToMMap[0]=mIndex;
//    mIndex[0]=0;
//    mIndex[1]=3;
//    vToMMap[1]=mIndex;
//    mIndex[0]=1;
//    mIndex[1]=2;
//    vToMMap[2]=mIndex;
//    numArray<double> mappedArray2DToVector(toMap.sliceArray2D(mTVMapping).sliceArray(vToMMap));
//    iomanagment::write(mappedArray2DToVector,cout<<"mappedArray2DToVector{(0,1)(0,3)(1,2) from 2 (5{3} 5{2} ) }=")<<endl;
// 
// //   numArraySlice<double>vs(v,1,4);
// //   numArraySlice<double>vs2(v2,1,4);
// 
//     std::vector<size_t> sizes(2);
//     sizes[0]=5;
//     sizes[1]=5;
// 
//     numArray2D<double> mapped2DBySplitting(toMap.sliceArray2D(sizes));
//     iomanagment::write(mapped2DBySplitting,cout<<"mapped2DBySplitting{split array to [5elem][5elm] from 2 (5{3} 5{2} ) }=")<<endl;


   //for(int i=0; i<vs.size();++i) cout<<vs[i]<<endl;
//   vs=vs+sin(vs2);


  return 0;
}


