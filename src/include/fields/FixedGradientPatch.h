#ifndef FIXEDGRADIENTPATCH_H
#define FIXEDGRADIENTPATCH_H

#include "PatchField.h"
#include "components/CmpTraits.h"

namespace SEM { namespace field {

template<typename T>
class FixedGradientPatch : public PatchField<T>
{
    typedef PatchField<T> baseType;

public:
    static std::string TYPE_NAME;
    
    using PatchField<T>::boundaryEdges;
    
//     FixedGradientPatch(iomanagment::RegistryFile::ref regFile, const mesh::Boundary &edges)
//         : baseType(regFile,edges)
//     {
//     }

    FixedGradientPatch( const mesh::Boundary &edges)
        : baseType(edges)
    {
    }

    FixedGradientPatch(const FixedGradientPatch<T>& other)
        : baseType(other)
    {
    }

    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    std::string typeName() const {return FixedGradientPatch<T>::TYPE_NAME;}
    std::string typeFamily() const {return NEUMANN;}
    
    // ----------------------------------------------------------------//
    //      EXPRESION TEMPLATES HANDLING
    // ----------------------------------------------------------------//
    /// copy constructor from expresion
    template<typename Exp>
    FixedGradientPatch(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }

    /// assign from single value
    FixedGradientPatch<T> & operator = (const T& singVal)
    {
        baseType::operator =(singVal);
        return *this;
    }

    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, FixedGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, FixedGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, FixedGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, FixedGradientPatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, FixedGradientPatch<T>,baseType)

};

template<typename T>
std::string FixedGradientPatch<T>::TYPE_NAME = "fixedGradient";

}//field
}//SEM


#endif // FIXEDGRADIENTPATCH_H
