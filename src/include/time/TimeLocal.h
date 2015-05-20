#ifndef TIMELOCAL_H
#define TIMELOCAL_H

#include "boost/filesystem/path.hpp"

#include "iomanagment/LocalPath.h"

namespace SEM {

class Time;

class TimeLocal : public iomanagment::LocalPath
{
    const Time& m_time;

    TimeLocal& operator = (const TimeLocal& other){}

public:
    TimeLocal(const Time& time);
    TimeLocal(const TimeLocal& other);

    ~TimeLocal();

    boost::filesystem::path path() const;
};

}//SEM



#endif // TIMELOCAL_H
