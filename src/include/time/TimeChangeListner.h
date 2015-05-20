#ifndef TIMECHANGELISTNER_H
#define TIMECHANGELISTNER_H

#include "components/Basic.h"

namespace SEM { namespace time {

    /// \brief The TimeChangeListner struct
    /// base interface for time change listners
    struct TimeChangeListner
    {
        virtual void timeChanged(Scalar time)=0;
    };

}//time
}//SEM




#endif // TIMECHANGELISTNER_H
