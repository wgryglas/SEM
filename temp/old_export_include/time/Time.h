
#ifndef _Time_H_
#define _Time_H_

#include <iostream>
#include <string>

#include "SolutionControl.h"

using namespace std;

namespace SEM
{
	class Time
	{
	public:

		Time();

		static string timeName();
		static bool end();
		static double time();

		void operator++();
		void operator++(int t);

	private :
		static double currTime;
		static int currTimeStep;


	};

	ostream& operator<<(ostream& o,Time& t);
}

#endif _Time_H_
