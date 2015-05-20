
#ifndef _f_H_
#define _f_H_

#include <Eigen\Dense>
#include "ScalarField.h"
#include "SolutionControl.h"
#include "Equation.h"


using namespace SEM::FIELDS;

namespace SEM
{
namespace FIELDMATH
{
	Eigen::VectorXd& f(ScalarField T);

}
}

#endif _f_H_