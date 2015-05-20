/*
#ifndef _ScalarField_H_
#define _ScalarField_H_

#include "Eigen/Dense"
#include <vector>

#include "Mesh.h"
#include "IOObject.h"
#include "Field.h"



using namespace SEM::MESH;
using namespace SEM::IOStream;

namespace SEM
{
namespace FIELDS
{
	class ScalarField: public Field
	{
	public:
		
		vector<double> Tacc;
		vector<double> Tprev;
		vector<bool> dirichletMask;
		map< boundaryEdge*, vector<double> > neumanValues;
		map<string,boundaryType> bType;

		ScalarField
		(
			IOFile& io,
			Mesh& mesh
		);
		


		ScalarField& operator=
		(
			Eigen::VectorXd& vector
		);

		~ScalarField
		(
		);

		ostream& displayNodal
		(
			ostream&
		);
		
		ostream& displayNodalWithCoords
		(
			ostream&
		);

		Eigen::VectorXd getAccSolution
		(
		);

		Eigen::VectorXd getPrevSolution
		(
		);


		void writeListWithNodes
		(
		);

		void writeList
		(
		);

		ScalarField operator+
		(
			ScalarField& f
		);
		
		ScalarField operator-
		(
			ScalarField& f
		);
		

	private:
		void ScalarField::applyInteriorInintialCondition
		(
			IOSubObject& io
		);

		void ScalarField::applyBoundaryCondition
		(
			IOSubObject& io
		);



	};


	ostream& operator<<(ostream& o, ScalarField field);

	
}
}



#endif _ScalarField_H_
*/
