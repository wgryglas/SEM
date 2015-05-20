#include <iostream>
#include <fstream>
#include <iterator>

#include "Eigen/Dense"

#include "utilities/foreach.h"
#include "utilities/TypeDefs.h"
#include "utilities/numArray.h"
#include "utilities/numArray2D.h"

#include "components/Basic.h"

#include "iomanagment/Dictionary2.h"
#include "iomanagment/FileProcessing.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"
#include "iomanagment/case.h"

#include "mesh/Mesh.h"

//#include "fields/GenericField.h"
//#include "fields/GeometricField.h"
//#include "fields/PatchField.h"

//#include "solver/EquationMatrix.h"
//#include "solver/Equation.h"
//#include "solver/Solver.h"

//#include "fieldMath/laplacian.h"
//#include "fieldMath/ddt.h"



#include "solver/Equation.h"
#include "solver/EquationMatrix.h"

int main(int argc, char* argv[])
{
    using namespace std;

    using namespace SEM;
    using namespace SEM::iomanagment;
//    using namespace SEM::field;
    using namespace SEM::mesh;

    Case::setup("/home/wojtek/Desktop/SEM_test");
    Mesh mesh;

    cout<<mesh<<endl;

//    GeometricField<Scalar> T(Case::time(),
//                             new RegistryFile(
//                                            Case::time(),
//                                            "T",
//                                            Case::time().localPath(),
//                                            MUST_READ,
//                                            AUTO),
//                            mesh
//                            );


//    LaplaceBuilder<double> eq=laplacian(1.,T);
//    las::SEMMatrix matrix;
//    numArray<Scalar> vector;
//    las::buildMatrixAndVector(eq,matrix,vector);

//    int nonZero=0;
//    for(int e=0; e<matrix.size(); ++e)
//    {
//        numArray<Scalar>::indexMapped eVec=vector.slice(mesh[e].indexVectorMask());
//        for(int row=0;row<matrix[e].size(); ++row)
//        {
//            for(int col=0; col<matrix[e][row].size();++col)
//            {
//                std::cout<<matrix[e][row][col]<<", ";
//                if(matrix[e][row][col]!=0)
//                    ++nonZero;
//            }
//            std::cout<<"\t"<<eVec[row]<<endl;
//        }
//        std::cout<<std::endl<<std::endl;
//    }

//    cout<<"nonZero elements: "<<nonZero<<std::endl;

//    las::solve(eq);


//    TimeFirstDerivativeBuilder eq1=ddt(T);
//    LaplaceBuilder<double> eq2=laplacian(1.,T);

//    las::CompoundEquationBuilder<Scalar, TimeFirstDerivativeBuilder,
//                                 LaplaceBuilder<double>,las::PlusEqAssigment> cEq(eq1,eq2);
//    las::solve( cEq );
//las::solve( eq2 );

//    iomanagment::write(T,cout)<<endl;

//    int i=0;
//    foreach (double& v, T)
//    {
//        if(v==1)
//            ++i;
//    }
//    cout<<"number of ones in rhsVecotr "<<i<<endl;


//    Dictionary d("T");
//    while(!Case::time().end())
//    {
//        T = sin(T);
//        T.relax(0.5);
//        d.clear();
//        T.write(d);
//        cout<<d<<endl;
//        ++Case::time();
//    }




//    GeometricField<double> geomF(mesh);
//    geomF=1;
//    cout <<endl;
//    iomanagment::write(geomF,cout<<"\n geomF[geomF=1]=\n")<<endl;


//    numArray2D<double> array(2,2);
//    array.row(1)=2;

//    numArray2D<double>::ColumnType col = array.column(1);

//    for(int i=0; i<col.size(); ++i)
//        cout<<"col("<<i<<") ="<<col[i]<<endl;

//    iomanagment::write(array,cout)<<endl;



    return 0;
}

