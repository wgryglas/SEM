
#ifndef _Equation_H_
#define _Equation_H_

#include <vector>
#include <Eigen/Sparse>

#include "ScalarField.h"

typedef Eigen::SparseMatrix<double> Matrix; //matirx type
typedef Eigen::VectorXd Vector; // vector type
typedef Eigen::Triplet<double> SMCoeff;//<-- to fill matrix-sparse matrix values


using namespace SEM::FIELDS;
using namespace std;
namespace SEM
{
namespace FIELDMATH
{
	class Equation
	{
	public:
		Matrix matrix;
		vector<SMCoeff> coeffs;
		Vector rhsVector;
		Vector initialSol;
		ScalarField & field;
		bool sparseBuild;

		Equation
		(
			ScalarField & field
		);

		void applyDirichletBoundaryCondition
		(
		);

		void buildSparseMatrix
		(
		);

		inline bool isBuild();

		inline void setAsUnBuild();

		Equation& operator+
		(
			Equation& eq
		);

		Equation& operator-
		(
			Equation& eq
		);


		Equation& operator=
		(
			Vector& v
		);

	};


	ostream& operator<<(ostream& o, Equation& e);

}
}






#endif _Equation_H_
