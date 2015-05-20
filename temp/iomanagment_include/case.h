#ifndef _Case_H_
#define _Case_H_

#include <iostream>
#include <string>

#include "boost/filesystem.hpp"

namespace SEM
{
	class Case
	{
	public:
        static boost::filesystem::path CASE_PATH;

		Case(int size, char* args[]);

        Case(std::string workingPath);

        boost::filesystem::path operator+ (const std::string& s);
        boost::filesystem::path operator+ (const boost::filesystem::path& p);

	};

    std::ostream& operator<<(std::ostream& o, Case );
}

#endif //_Case_H_
