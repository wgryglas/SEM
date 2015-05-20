
#ifndef _Solver_H_
#define _Solver_H_

#include"Equation.h"


namespace SEM
{
namespace FIELDMATH
{

	Eigen::VectorXd& solve(Equation& equation);

}
}







#endif _Solver_H_
