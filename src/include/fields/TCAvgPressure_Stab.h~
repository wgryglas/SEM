#ifndef TCAVGPRESSURE_STAB_H
#define TCAVGPRESSURE_STAB_H

#include "TCAvgPressure.h"

namespace SEM{
namespace field{

class TCAvgPressure_Stab : public TCAvgPressure
{
  
    Scalar m_delta;
    Scalar m_Uref;
    
public:
    static std::string TYPE_NAME;
    
    TCAvgPressure_Stab(const iomanagment::Register &reg, const mesh::Boundary &edges);
    
    TCAvgPressure_Stab(const TCAvgPressure_Stab & other)=delete;
    TCAvgPressure_Stab & operator=(const TCAvgPressure_Stab & other)=delete;
    
    
    virtual ~TCAvgPressure_Stab();
    
    std::string typeName() const;
    
    void updateValues();
    
    //IO
    void read(const iomanagment::Dictionary &fieldDict);
    
    void write(iomanagment::Dictionary &fieldDict);
    
    inline Scalar delta() const;
    
    inline Scalar Uref() const;
};


Scalar TCAvgPressure_Stab::delta() const
{
    return m_delta;
}

Scalar TCAvgPressure_Stab::Uref() const
{
    return m_Uref;
}


} //field
} //SEM

#endif // TCAVGPRESSURE_STAB_H
