
#ifndef TIMECHANGINGVALUEPATCH_H
#define TIMECHANGINGVALUEPATCH_H

#include <vector>
#include <array>
#include <sstream>

#include "fields/FixedValuePatch.h"
#include "components/CmpTraits.h"
#include "components/Basic.h"
#include <iomanagment/case.h>
#include "utilities/LinearInterpolator.h"

namespace SEM { namespace field { 

    
template<typename T>
class TimeChangingValuePatch : public FixedValuePatch<T>
{
    LinearInterpolator<Scalar,T> m_values;
    typedef FixedValuePatch<T> baseType;
public:
    static std::string TYPE_NAME;
    
    using baseType::name;
    
    TimeChangingValuePatch(const mesh::Boundary &edges)
    : baseType(edges)
    {
    }
    
    TimeChangingValuePatch(const TimeChangingValuePatch<T>& other)
    : baseType(other), m_values(other.m_values)
    {
    }
    
    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    std::string typeName() const {return TimeChangingValuePatch<T>::TYPE_NAME;}
    
    /// copy constructor from expresion
    template<typename Exp>
    TimeChangingValuePatch(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }
    
    /// assign from single value
    TimeChangingValuePatch<T> & operator = (const T& singVal)
    {
        baseType::operator =(singVal);
        return *this;
    }
    
    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, TimeChangingValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, TimeChangingValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, TimeChangingValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, TimeChangingValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, TimeChangingValuePatch<T>,baseType)
    
    void timeChanged(Scalar time)
    {
        (*this) = m_values(time);
    }
    
    virtual void read(const iomanagment::Dictionary &fieldDict)
    {
        m_values.update(fieldDict.subDictionary("boundary").subDictionary(name()).entry("value"));
        (*this) = m_values(Case::time().time() + Case::time().timeStep());
    }
    
    virtual void write(iomanagment::Dictionary &fieldDict)
    {
        iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
        
        iomanagment::DictEntry *valEntry =new iomanagment::DictEntry("value");
        m_values.writeValues(*valEntry);
        
        iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
        bDict->add(typeEntry);
        bDict->add(valEntry);
        
        if(!fieldDict.hasSubDictionary("boundary"))
            fieldDict.add(new iomanagment::Dictionary("boundary"));
        
        fieldDict.subDictionary("boundary").add(bDict);
    }
};


template<typename T>
std::string TimeChangingValuePatch<T>::TYPE_NAME = "timeChangingValue";


}//field
}//SEM

#endif // TIMECHANGINGVALUEPATCH_H
