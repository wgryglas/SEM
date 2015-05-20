
#ifndef _ddt_H_
#define _ddt_H_



#include "Equation.h"
#include "ScalarField.h"
#include "SolutionControl.h"


using namespace SEM::FIELDS;

namespace SEM
{
namespace FIELDMATH
{

	Equation& ddt
	(
		ScalarField& field
	);


	void Euler
	(
		Equation& e
	);

	void BDF2
	(
		Equation& e
	);

}
}






#endif _ddt_H_
