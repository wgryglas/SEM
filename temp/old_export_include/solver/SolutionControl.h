
#ifndef _SolutionControl_H_
#define _SolutionControl_H_

#include "IOFile.h"

using namespace SEM::IOStream;
namespace SEM
{
	enum TimeDiscretization
	{
		Euler,
		BDF2
	};

	enum SolverType
	{
		ConjugateGradient,
		BiCGSTAB,
		SimplicialLLT,
		SimplicialLDLT
	};

	class SolutionControl
	{
	public:
		static SolverType solverType;
		static TimeDiscretization timeDiscretization;
		static int maxIteration;
		static double timeStep;
		static int numberOfTimeSteps;


		SolutionControl
		(
			IOFile& io
		);

		SolverType str2SolverType
		(
			string s
		);

		TimeDiscretization str2TimeDiscretization
		(
			string s
		);

	};
}







#endif _SolutionControls_H_
