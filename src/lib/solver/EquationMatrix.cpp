
//#include"EquationMatrix.h"

//namespace SEM
//{
//namespace FIELDMATH
//{

//    EquationMatrix::EquationMatrix()
//    {
//    }

//    EquationMatrix::EquationMatrix(int size): mat(Eigen::MatrixXd(size,size))
//	{
//    }

//    double& EquationMatrix::operator()(int i, int j)
//	{
//		return mat(i,j);
//    }

//    void EquationMatrix::assign(vector<vector<double>>& subMat, int rowStart, int colStart)
//	{
//		for(int i=0; i<subMat.size(); i++)
//		{
//			for(int j=0;j<subMat[i].size(); j++)
//			{
//				mat(i+rowStart,j+colStart) = subMat[i][j];
//			}
//		}
//    }

//    void EquationMatrix::asign(vector<vector<double>>& subMat,vector< vector< array<int,2> > >& indexes)
//	{
//		for(int i=0;i<indexes.size();i++)
//		{
//			for(int j=0;j<indexes[i].size();j++)
//			{
//				mat(indexes[i][j][0],indexes[i][j][1]) = subMat[i][j];
//			}
//		}
//    }

//    void EquationMatrix:: add(vector<vector<double>>& subMat,vector< vector< array<int,2> > >& indexes)
//	{
//		for(int i=0;i<indexes.size();i++)
//		{
//			for(int j=0;j<indexes[i].size();j++)
//			{
//				mat(indexes[i][j][0],indexes[i][j][1]) += subMat[i][j];
//			}
//		}
//    }

//    void EquationMatrix::clear()
//    {
//		//mat.
//    }

//}

//}
