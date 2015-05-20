#ifndef TIMER_H
#define TIMER_H

#include <ctime>

namespace SEM { namespace time {

class Timer
{
    std::clock_t m_startClock;

public:
    Timer();
    virtual ~Timer(){}

    Timer & restart();

    double elapsed() const;
};


}//time
}//SEM

#endif // TIMER_H
