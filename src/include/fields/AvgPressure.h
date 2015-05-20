#ifndef AVGPRESSURE_H
#define AVGPRESSURE_H

#include <string>

#include "components/Basic.h"
#include "PatchField.h"


namespace SEM { namespace field { 

template<typename T>
class GeometricField;
    

extern const std::string PRESSURE_DO_NOTHING;
    
class AvgPressure: public PatchField<Scalar>
{
protected:
    typedef PatchField<Scalar> baseType;
    Scalar m_avgValue;
    
    const GeometricField<Vector> & m_velocity;
    Scalar m_viscosity;
    
public:
    static std::string TYPE_NAME;
    
    AvgPressure(const iomanagment::Register &reg, const mesh::Boundary &edges);
    //AvgPressure(const mesh::Boundary & edges);

    AvgPressure(const AvgPressure & other)=delete;
    AvgPressure & operator=(const AvgPressure & other)=delete;
    
    
    ~AvgPressure();
    
    void updateValues();
    
    void timeChanged(Scalar time);
    
    Scalar avgPressure() const;
    
    std::string typeName() const;
    
    std::string typeFamily() const;
    
    //IO
    void read(const iomanagment::Dictionary &fieldDict);
    
    void write(iomanagment::Dictionary &fieldDict);
    
};


}
}

#endif // AVGPRESSURE_H
