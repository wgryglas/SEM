#ifndef TIMECHANGINGGRADIENTPATCH_H
#define TIMECHANGINGGRADIENTPATCH_H

#include "fields/TimeChangingValuePatch.h"

namespace SEM { namespace field {

    
template<typename T>    
class TimeChangingGradientPatch : public TimeChangingValuePatch<T>
{
    typedef TimeChangingValuePatch<T> baseType;
public:
    static const std::string TYPE_NAME;
    
    TimeChangingGradientPatch(const mesh::Boundary &edges)
    : baseType(edges)
    {
    }
    
    TimeChangingGradientPatch(const TimeChangingGradientPatch<T>& other)
    : baseType(other)
    {
    }
    
    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    std::string typeName() const {return TimeChangingGradientPatch<T>::TYPE_NAME;}
    std::string typeFamily() const { return NEUMANN;}
    
    template<typename Exp>
    TimeChangingGradientPatch(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }
    
    /// assign from single value
    TimeChangingGradientPatch<T> & operator = (const T& singVal)
    {
        baseType::operator =(singVal);
        return *this;
    }
    
    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, TimeChangingGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, TimeChangingGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, TimeChangingGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, TimeChangingGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, TimeChangingGradientPatch<T>,baseType)
};

template<typename T>
const std::string TimeChangingGradientPatch<T>::TYPE_NAME = "timeChangingGradient";

}//field
}//SEM

#endif // TIMECHANGINGGRADIENTPATCH_H
