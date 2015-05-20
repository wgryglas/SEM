#include "TimeLocal.h"

#include "time/Time.h"
#include "iomanagment/case.h"

namespace SEM {

TimeLocal::TimeLocal(const Time &time):m_time(time)
{
}

TimeLocal::TimeLocal(const TimeLocal &other):m_time(other.m_time)
{
}

TimeLocal::~TimeLocal()
{
}

boost::filesystem::path TimeLocal::path() const
{
    return m_time.timeName();
}

}//SEM
