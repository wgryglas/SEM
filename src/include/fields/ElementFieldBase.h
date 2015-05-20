#ifndef ELEMENTFIELDBASE_H
#define ELEMENTFIELDBASE_H

#include "mesh/Mesh.h"

namespace SEM { namespace field {

template<typename T,typename Derived>
struct ElementFieldTraits
{
    typedef numArray<T> return_type;
};
    
template<typename T, typename Derived> 
class ElementFieldBase 
{
    const mesh::Mesh & m_mesh;
public:
    typedef typename ElementFieldTraits<T,Derived>::return_type return_type;
    
    ElementFieldBase(const mesh::Mesh & mesh) : m_mesh(mesh)
    {
    }
    
    virtual ~ElementFieldBase(){};
    
    const mesh::Mesh & mesh() const 
    {
        return m_mesh;
    }
    
    size_t elementsNumber() const
    {
        return m_mesh.size();
    }
    
    return_type element(size_t e) const
    {
        return static_cast<const Derived*>(this)->element(e);
    }
    
};

#define SEM_MAKE_ELEMENT_FIELD_BINARY_FUNCTION(CLASS_NAME, OP)                                         \
template<typename T, typename Lhs, typename Rhs>                                                        \
class CLASS_NAME : public ElementFieldBase<T,CLASS_NAME<T,Lhs,Rhs> >                                    \
{                                                                                                       \
    const ElementFieldBase<T,Lhs> & m_lhs;                                                              \
    const ElementFieldBase<T,Rhs> & m_rhs;                                                              \
    typedef ElementFieldBase<T,CLASS_NAME<T,Lhs,Rhs> > baseType;                                        \
    typedef typename ElementFieldTraits< T, CLASS_NAME< T, Rhs, Lhs > >::return_type return_type;       \
public:                                                                                                 \
    CLASS_NAME(const ElementFieldBase<T,Lhs> & lhs,const ElementFieldBase<T,Rhs> & rhs)                 \
    : baseType(lhs.mesh()), m_lhs(lhs), m_rhs(rhs)                                                      \
    {                                                                                                   \
    }                                                                                                   \
    return_type element(size_t e) const                                                                 \
    {                                                                                                   \
        return  m_lhs.element(e) OP m_rhs.element(e);                                                   \
    }                                                                                                   \
};                                                                                                      \
template<typename T, typename Lhs, typename Rhs>                                                        \
CLASS_NAME<T,Lhs,Rhs>                                                                                   \
operator OP                                                                                             \
(const ElementFieldBase<T,Lhs> & lhs,const ElementFieldBase<T,Rhs> & rhs)                               \
{                                                                                                       \
    return CLASS_NAME<T,Lhs,Rhs>(lhs,rhs);                                                              \
}                                                                                                       

SEM_MAKE_ELEMENT_FIELD_BINARY_FUNCTION(ElementFieldSum, +)
SEM_MAKE_ELEMENT_FIELD_BINARY_FUNCTION(ElementFieldSub, -)
SEM_MAKE_ELEMENT_FIELD_BINARY_FUNCTION(ElementFieldMul, *)
SEM_MAKE_ELEMENT_FIELD_BINARY_FUNCTION(ElementFieldDiv, /)




#define SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(CLASS_NAME, FUNCTION_NAME)                                \
template<typename T, typename Arg>                                                                      \
class CLASS_NAME : public ElementFieldBase<T,CLASS_NAME<T,Arg> >                                        \
{                                                                                                       \
    const ElementFieldBase<T,Arg> & m_arg;                                                              \
    typedef ElementFieldBase<T,CLASS_NAME<T,Arg> > baseType;                                            \
    typedef typename ElementFieldTraits< T, CLASS_NAME< T,Arg > >::return_type return_type;             \
public :                                                                                                \
    CLASS_NAME(const ElementFieldBase<T,Arg> & arg)                                                     \
    : baseType(arg.mesh()), m_arg(arg)                                                                  \
    {                                                                                                   \
    }                                                                                                   \
    return_type element(size_t e) const                                                                 \
    {                                                                                                   \
        return FUNCTION_NAME(m_arg.element(e));                                                         \
    }                                                                                                   \
};                                                                                                      \
template<typename T, typename Arg>                                                                      \
CLASS_NAME<T,Arg> FUNCTION_NAME(const ElementFieldBase<T,Arg> & arg)                                    \
{                                                                                                       \
    return CLASS_NAME<T,Arg>(arg);                                                                      \
}                                                                                                       

SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldNegate, operator-)
SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldPositive, operator+)
SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldSin, sin)
SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldCos, cos)
SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldAbs, abs)
SEM_MAKE_ELEMENT_FIELD_UNARY_FUNCTION(ElementFieldExp, exp)


#define SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(CLASS_NAME,S_TYPE,OP)                    \
template<typename T, typename Field>                                                                    \
class CLASS_NAME##LeftOperand :                                                                         \
public ElementFieldBase< typename array::ConversionTraits<S_TYPE,T>::return_type, CLASS_NAME##LeftOperand<T,Field> > \
{                                                                                                       \
    typedef typename array::ConversionTraits<S_TYPE,T>::return_type type;                               \
    typedef ElementFieldBase<type,CLASS_NAME##LeftOperand<T,Field> > baseType;                          \
    typedef typename ElementFieldTraits< type, CLASS_NAME##LeftOperand< T,Field > >::return_type return_type;\
    const ElementFieldBase<T,Field> & m_field;                                                          \
    const S_TYPE m_val;                                                                                 \
public:                                                                                                 \
    CLASS_NAME##LeftOperand(const S_TYPE & val, const ElementFieldBase<T,Field> & field)                \
    : baseType(field.mesh()), m_val(val), m_field(field)                                                \
    {                                                                                                   \
    }                                                                                                   \
    return_type element(size_t e) const                                                                 \
    {                                                                                                   \
        return m_val OP m_field.element(e);                                                             \
    }                                                                                                   \
    S_TYPE singleValue() const { return m_val;}                                                         \
    const Field & field() const { return static_cast<const Field&>(m_field);}                           \
};                                                                                                      \
template<typename T, typename Field>                                                                    \
class CLASS_NAME##RightOperand :                                                                        \
public ElementFieldBase< typename array::ConversionTraits<S_TYPE,T>::return_type, CLASS_NAME##RightOperand<T,Field> > \
{                                                                                                       \
    typedef typename array::ConversionTraits<S_TYPE,T>::return_type type;                               \
    typedef ElementFieldBase<type,CLASS_NAME##RightOperand<T,Field> > baseType;                         \
    typedef typename ElementFieldTraits< type, CLASS_NAME##RightOperand< T,Field > >::return_type return_type;\
    const ElementFieldBase<T,Field> & m_field;                                                          \
    const S_TYPE m_val;                                                                                 \
public:                                                                                                 \
    CLASS_NAME##RightOperand(const S_TYPE & val, const ElementFieldBase<T,Field> & field)               \
    : baseType(field.mesh()), m_val(val), m_field(field)                                                \
    {                                                                                                   \
    }                                                                                                   \
    return_type element(size_t e) const                                                                 \
    {                                                                                                   \
        return m_field.element(e) OP m_val;                                                             \
    }                                                                                                   \
    S_TYPE singleValue() const { return m_val;}                                                         \
    const Field & field() const { return static_cast<const Field&>(m_field);}                           \
};                                                                                                      \
template<typename T, typename Field>                                                                    \
CLASS_NAME##LeftOperand<T,Field> operator OP(const S_TYPE& val, const ElementFieldBase<T,Field> & field)\
{                                                                                                       \
    return CLASS_NAME##LeftOperand<T,Field>(val,field);                                                 \
}                                                                                                       \
template<typename T, typename Field>                                                                    \
CLASS_NAME##RightOperand<T,Field> operator OP(const ElementFieldBase<T,Field> & field, const S_TYPE& val)\
{                                                                                                       \
    return CLASS_NAME##RightOperand<T,Field>(val,field);                                                \
}                                                                                                       
                        
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldScalarSum, Scalar, +)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldScalarSub, Scalar, -)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldScalarMul, Scalar, *)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldScalarDiv, Scalar, /)

SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldVectorSum, Vector, +)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldVectorSub, Vector, -)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldVectorMul, Vector, *)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldVectorDiv, Vector, /)

SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldTensorSum, Tensor, +)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldTensorSub, Tensor, -)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldTensorMul, Tensor, *)
SEM_MAKE_ELEMNT_FILED_SINGLE_VALUE_IN_BINARY_FUNCTIONS(ElementFieldTensorDiv, Tensor, /)



#define SEM_MAKE_ELEMENT_FIELD_ETARRAY_OPERATOR_SUPPORT(CLASS_NAME, OP)                                 \
template<typename T, typename Field, typename ET>                                                       \
class CLASS_NAME##LeftContinous : public ElementFieldBase<T,CLASS_NAME##LeftContinous<T,Field,ET> >     \
{                                                                                                       \
    const array::ET_Array<T,ET> & m_list;                                                               \
    const ElementFieldBase<T,Field> & m_field;                                                          \
    typedef ElementFieldBase<T,CLASS_NAME##LeftContinous<T,Field,ET> > baseType;                        \
public:                                                                                                 \
    using baseType::mesh;                                                                               \
    CLASS_NAME##LeftContinous(const array::ET_Array<T,ET> & list, const ElementFieldBase<T,Field> & field)\
    : baseType(field.mesh()), m_field(field), m_list(list)                                              \
    {                                                                                                   \
    }                                                                                                   \
    numArray<T> element(size_t e) const                                                                 \
    {                                                                                                   \
        const auto& mask = mesh()[e].indexVectorMask();                                                 \
        auto fieldElement = m_field.element(e);                                                         \
        numArray<T> result(mask.size());                                                                \
        for(size_t n=0; n<mask.size(); ++n)                                                             \
        {                                                                                               \
            result[n]=m_list.evalAt(mask[n]) OP fieldElement.evalAt(n);                                 \
        }                                                                                               \
        return result;                                                                                  \
    }                                                                                                   \
};                                                                                                      \
template<typename T, typename Field, typename ET>                                                       \
class CLASS_NAME##RightContinous : public ElementFieldBase<T,CLASS_NAME##RightContinous<T,Field,ET> >   \
{                                                                                                       \
    const array::ET_Array<T,ET> & m_list;                                                               \
    const ElementFieldBase<T,Field> & m_field;                                                          \
    typedef ElementFieldBase<T,CLASS_NAME##RightContinous<T,Field,ET> > baseType;                       \
    public:                                                                                             \
    using baseType::mesh;                                                                               \
    CLASS_NAME##RightContinous(const array::ET_Array<T,ET> & list, const ElementFieldBase<T,Field> & field)\
    : baseType(field.mesh()), m_field(field), m_list(list)                                              \
    {                                                                                                   \
    }                                                                                                   \
    numArray<T> element(size_t e) const                                                                 \
    {                                                                                                   \
        const auto& mask = mesh()[e].indexVectorMask();                                                 \
        auto fieldElement = m_field.element(e);                                                         \
        numArray<T> result(mask.size());                                                                \
        for(size_t n=0; n<mask.size(); ++n)                                                             \
        {                                                                                               \
            result[n]= fieldElement.evalAt(n) OP m_list.evalAt(mask[n]);                                \
        }                                                                                               \
        return result;                                                                                  \
    }                                                                                                   \
};                                                                                                      \
template<typename T, typename Field, typename ET>                                                       \
CLASS_NAME##LeftContinous<T,Field,ET> operator OP                                                       \
(const array::ET_Array<T,ET> & list, const ElementFieldBase<T,Field> & field)                           \
{                                                                                                       \
    return CLASS_NAME##LeftContinous<T,Field,ET>(list,field);                                           \
}                                                                                                       \
template<typename T, typename Field, typename ET>                                                       \
CLASS_NAME##RightContinous<T,Field,ET> operator OP                                                      \
(const ElementFieldBase<T,Field> & field,const array::ET_Array<T,ET> & list)                            \
{                                                                                                       \
    return CLASS_NAME##RightContinous<T,Field,ET>(list,field);                                          \
}                                                                                                       

SEM_MAKE_ELEMENT_FIELD_ETARRAY_OPERATOR_SUPPORT(ElementAndContinousFieldSum,+)
SEM_MAKE_ELEMENT_FIELD_ETARRAY_OPERATOR_SUPPORT(ElementAndContinousFieldDiff,-)
SEM_MAKE_ELEMENT_FIELD_ETARRAY_OPERATOR_SUPPORT(ElementAndContinousFieldMul,*)
SEM_MAKE_ELEMENT_FIELD_ETARRAY_OPERATOR_SUPPORT(ElementAndContinousFieldDiv,/)






} //field 
} //SEM

#endif // ELEMENTFIELDBASE_H
