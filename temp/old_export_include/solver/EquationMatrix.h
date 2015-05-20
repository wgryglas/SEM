
#ifndef _EquationMatrix_H_
#define _EquationMatrix_H_

#include<Eigen/Dense>
#include<vector>
#include<array>

using namespace std;
namespace SEM
{
namespace FIELDMATH
{
	class EquationMatrix
	{
	public:
		Eigen::MatrixXd mat;

		EquationMatrix
		(
		);

		EquationMatrix
		(
			int size
		);

		double& operator()
		(
			int i,
			int j
		);

		void assign
		(
			vector<vector<double>>& subMat,
			int rowStart,
			int colStart
		);

		void asign
		(
			vector<vector<double>>& subMat,
			vector<vector<array<int,2>>>& indexes
		);

		void add
		(
			vector<vector<double>>& subMat,
			vector< vector< array<int,2> > >& indexes
		);

		void clear
		(
		);

	};
}
}




#endif _EquationMatrix_H_
