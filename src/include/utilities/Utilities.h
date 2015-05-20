#ifndef _Utilities_H_
#define _Utilities_H_

#include <map>
#include <cmath>
#include "utilities/TypeDefs.h"
//#include <applications/poisson/createFileds.h>

namespace SEM
{

    template<typename itrType> 
    void deleteAll(itrType  begin, itrType end)
    {
        for(; begin!=end; ++begin)
        {
            delete *begin;
        }
    }
    
    template<typename Container>
    void deleteAll(const Container &c)
    {
        deleteAll(c.begin(),c.end());
    }
    

    template<typename T1, typename T2>
    void deleteAllKeys(const std::map<T1,T2> & map)
    {
        for(auto e : map)
            delete map.first;
    }
    
    template<typename T1, typename T2>
    void deleteAllValues(const std::map<T1,T2> & map)
    {
        for(auto e : map)
            delete map.second;
    }
    
//     inline
//     bool isValueValid(const Scalar &value)
//     {
//         return std::isnan(value);
//     }
    
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

//	void operator+=(vector< vector<vector<vector<double>>>>& a, vector< vector<vector<vector<double>>>>& b);

#define  SEM_UNUSED(VARIABLE) (void) VARIABLE;


}



#endif //_Utilities_H_
