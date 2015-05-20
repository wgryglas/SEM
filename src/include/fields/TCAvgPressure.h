
#ifndef TCAVGPRESSURE_H
#define TCAVGPRESSURE_H

#include "TimeChangingValuePatch.h"
#include "AvgPressure.h"

namespace SEM { namespace field { 

class TCAvgPressure : public AvgPressure
{
protected:
    LinearInterpolator<Scalar,Scalar> m_values;
public:
    static std::string TYPE_NAME;
    
    TCAvgPressure(const iomanagment::Register &reg, const mesh::Boundary &edges);
    
    TCAvgPressure(const TCAvgPressure & other)=delete;
    TCAvgPressure & operator=(const TCAvgPressure & other)=delete;
    
    
    ~TCAvgPressure();
    
    std::string typeName() const;
    
    void timeChanged(Scalar time);
    
    //IO
    void read(const iomanagment::Dictionary &fieldDict);
    
    void write(iomanagment::Dictionary &fieldDict);
    
};

}//field
}//SEM

#endif // TCAVGPRESSURE_H
