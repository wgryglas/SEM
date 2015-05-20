#ifndef FIXEDVALUEPATCH_H
#define FIXEDVALUEPATCH_H

#include "PatchField.h"

namespace SEM { namespace field {

template<typename T>
class FixedValuePatch : public PatchField<T>
{
    typedef PatchField<T> baseType;
public:
    static std::string TYPE_NAME;
    
//     FixedValuePatch(iomanagment::RegistryFile::ref regFile,const mesh::Boundary &edges)
//         : baseType(regFile,edges)
//     {
//     }

    FixedValuePatch(const mesh::Boundary &edges)
        : baseType(edges)
    {
    }

    FixedValuePatch(const FixedValuePatch<T>& other)
        : baseType(other)
    {
    }

    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    
    std::string typeName() const {return FixedValuePatch<T>::TYPE_NAME;}
    std::string typeFamily() const {return DIRICHLET;}
    
    // ----------------------------------------------------------------//
    //      EXPRESION TEMPLATES HANDLING
    // ----------------------------------------------------------------//
    /// copy constructor from expresion
    template<typename Exp>
    FixedValuePatch(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }

    /// assign from single value
    FixedValuePatch<T> & operator = (const T& singVal)
    {
        baseType::operator =(singVal);
        return *this;
    }
    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, FixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, FixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, FixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, FixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, FixedValuePatch<T>,baseType)

};

template<typename T>
std::string FixedValuePatch<T>::TYPE_NAME = "fixedValue";

}//field
}//SEM





#endif // FIXEDVALUEPATCH_H
