#include "Timmer.h"

#include "iomanagment/InfoStream.h"

namespace SEM { namespace time {

Timer::Timer(): m_startClock(std::clock())
{
}

Timer &Timer::restart()
{
    m_startClock = std::clock();
}


double Timer::elapsed() const
{
    double diff = std::clock()-m_startClock;
    return (diff)/CLOCKS_PER_SEC;
}


}//time
}//SEM

