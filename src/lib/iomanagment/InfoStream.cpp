
#include "InfoStream.h"

namespace SEM { namespace iomanagment {
    std::ostream& INFO(std::ostream& o)
	{
        o<<"SEM::INFO::";

		return o;
	}

    std::ostream& arrow(std::ostream &o)
    {
        return o<<"------------>";
    }

    std::ostream& endInfo( std::ostream& o )
	{
        o<<std::endl<<"---------------------------------------"<<std::endl;
		return o;
	}

    std::ostream& ERROR_IN_FUNCTION(std::ostream& o, const char* function, const char* file, int line )
    {
        o<<"SEM::FATAL ERROR!"<<std::endl
         <<arrow<<" in file: "<<file<<", in function: "<<function<<", in line: "<<line<<std::endl;

        return o;
    }

    std::ostream& WARNING(std::ostream & o)
    {
        o<<"SEM::WARNING"<<std::endl
         <<arrow;

        return o;
    }

    std::ostream& endProgram(std::ostream& o)
    {
        o<<std::endl<<"---------------------------------------"<<std::endl<<
                      "PROGRAM EXITING"<<std::endl;
        o.flush();
        exit(1);

        return o;
    }
}//iomanagment
}//SEM
