#ifndef _InfoStream_H_
#define _InfoStream_H_

#include <stdlib.h>
#include <iostream>
#include <iomanip>

namespace SEM{ namespace iomanagment{
    //////////////////////////////////////////////
    ///  PROGRAM-USER COMMUNICATION UTILITIES
    /// ----------------------------------------
    /// set of manipulator and predefined function
    /// for easy to use outputing information for
    /// users
    //////////////////////////////////////////////

    /** INFO -output prepared info text
     * @param o -output stream
     */
    std::ostream& INFO(std::ostream& o);

    std::ostream& arrow(std::ostream &o);

    std::ostream& endInfo(std::ostream&);

    std::ostream& ERROR_IN_FUNCTION(std::ostream& o, const char* function, const char* file, int line );

    std::ostream& WARNING(std::ostream& o);

    std::ostream& endProgram(std::ostream& o);

}//IOMANAGMENT
}//SEM

#define OUTPUT_STREAM std::cout
#define ERROR_STREAM  std::cerr

#define Warning WARNING(OUTPUT_STREAM)

#define ErrorInFunction \
    SEM::iomanagment::ERROR_IN_FUNCTION(ERROR_STREAM, __func__, __FILE__, __LINE__)

#define Info SEM::iomanagment::INFO(OUTPUT_STREAM)

#define InfoArrow SEM::iomanagment::arrow(OUTPUT_STREAM)

#define EndInfo SEM::iomanagment::endInfo(OUTPUT_STREAM)


#endif //_InfoStream_H_
