
#ifndef _Utilities_H_
#define _Utilities_H_

#include<vector>

using namespace std;

namespace SEM
{
	/*class Utils
	{
	public: 
		static IOStream::InfoStream Info;

	};*/
	
	
	/*template<typename T>
	bool str2num
	(
		T& num,
		const std::string& s
	)
	{
		std::istringstream iss(s);
		return !(iss >>std::dec >>num).fail();
	};*/

	/*template<typename T>
	T str2num(std::string& s)
	{
		std::ifastream <basic_formatters, string_reader> myString(&s);
		T value;
		myString >> value;

		return value;
	}*/

	void operator+=(vector< vector<vector<vector<double>>>>& a, vector< vector<vector<vector<double>>>>& b);


}



#endif _Utilities_H_