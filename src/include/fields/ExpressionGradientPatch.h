#ifndef EXRPESSIONGRADIENTPATCH_H
#define EXRPESSIONGRADIENTPATCH_H

#include "ExpressionValuePatch.h"

namespace SEM { namespace field {

template<typename T>
class ExpressionGradientPatch : public ExpressionValuePatch<T>
{
    typedef ExpressionValuePatch<T> baseType;
public:
    static std::string TYPE_NAME;

    ExpressionGradientPatch(const mesh::Boundary &edges, boost::array<std::string,CmpTraits<T>::DIM_SIZE> & expressions)
    : baseType(edges,expressions)
    {
    }
    
    ExpressionGradientPatch(const mesh::Boundary &edges)
        : baseType(edges)
    {
    }

    ExpressionGradientPatch(const ExpressionGradientPatch<T>& other)
        : baseType(other)
    {
    }

    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    std::string typeName() const {return ExpressionGradientPatch<T>::TYPE_NAME;}
    std::string typeFamily() const {return NEUMANN;}
};

template<typename T>
std::string ExpressionGradientPatch<T>::TYPE_NAME = "expressionGradient";

}//field
}//SEM





#endif // EXRPESSIONGRADIENTPATCH_H
