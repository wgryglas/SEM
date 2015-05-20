#ifndef ET_Array_H
#define ET_Array_H

#include "components/Basic.h"

namespace SEM { namespace array {

///\class ET_Array
/// Base class for CTRP pattern. It's base class
/// which holds expression - all "derived expr."
/// derives it, so in fact any expression is
/// equivalent of this base class.
/// CTPR allows here(by static polymorphism) to
/// treat any DerivedExpression as ET_Array.
/// Note that, no virtual fun. is used in this
/// implementation.
/// -----------------------------------------------
/// Any container which would you like to be supported
/// by those expr. templates need to extend ET_Array
/// (this allows templates to figure out list element type,
/// and treat your container just like another expression)
/// like folowing:
/// template<typename T>
/// Contianer : public ET_Array<T,Container<T> >
/// {
///     ...
///      inline size_t size() const {...}-->return cont. size
///      inline T evalAt(size_t index) const {...}-->return value at index
///     ...
/// }
/// Contianer need also to implement (at least):
/// template<AnyType>
/// operator=(const ET_Array<T,AnyType>& compundExpr)
/// {
///     --> assigment from exrp. by evalAt method into
///         desired data container.
/// }
/// -----------------------------------------------
/// If you wan't to extend some list, which already
/// supports ET, then you are obligated only to place
/// 2 member function in your new class:
///
/// Copy constructor:
/// template<typename T,typename AnyExpr>
/// class_name(ET_Array<T,AnyExpr> & expr)
/// {
///   ---> assigment from ET_Array by evalAt and size methods
/// }
///
/// Assigment operator
/// template<typename T,typename AnyExpr>
/// class_name& operator=(ET_Array<T,AnyExpr> & expr)
/// {
///   ---> assigment from ET_Array by evalAt and size methods
/// }
///
/// That's all, now your new class shall work with any
/// ET.
///
/// If you would like to support also +=,-= etc., then
/// implement othere operators for example like below:
/// template<typename T,typename AnyExpr>
/// class_name& operator **** (ET_Array<T,AnyExpr> & expr)
/// {
///   ...
///    (*this)[i] **** expr.evalAt(i);
///   ...
/// }
/// (where **** is some sign of operator)
/// -----------------------------------------------
/// WARNING: Using single-value (of list element type)
/// in expressions is allowd(eg. List1+List2+1 or
/// List1+2+sin(List2)). But when you want do assigment
/// List = 5, then there is no solution to auto-convert
/// singleValue to "singleValueExpresion", which would
/// be handled by List, so if you want your list to
/// handle such a assigment then unfortunately you have
/// to implement this conversin inside your List class
/// definition. Definition of such operator is easy,
/// eg.:
/// class_name<T> & operator = (const T &sigleVal)
/// {
///   for(int i=0;i<size();++i)(*this)[i]=singleVal;
/// }
/// -----------------------------------------------
/// All derived expr. templ. implements above 2
/// inline methods(size,evalAt) and as
/// result of evalAt - returns evaluated value
/// according to it's specification.
/// -----------------------------------------------
/// Expr. function and ScalarExpr do no't take 2 arg., but
/// only 1, and basing on it implement apropriate return value.
/// Notice that ScalarExpr is also evaluated at "=" operator,
/// and behave just like all derived expresions- it allows
/// to keep all expressions in the same manner.
//////////////////////////////////////////////////////////////////
template<typename T, typename ET_DerivedExpresion >
class ET_Array
{
public:
    inline size_t expSize() const
    {
        return static_cast<const ET_DerivedExpresion*>(this)->expSize();
    }

    inline T evalAt(size_t index) const
    {
        return static_cast<const ET_DerivedExpresion*>(this)->evalAt(index);
    }
};

//////////////////////////////////////////////////////////////////
/// \macro MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR
/// Generator for template derived expressions with
/// linked binary operator function. Generated expr.
/// together with operator fun. behave as any of +-*/
/// operator.
//////////////////////////////////////////////////////////////////
#define MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR(OP,CLASS_NAME)                      \
    template<typename T, typename T1, typename T2>                              \
    class CLASS_NAME : public ET_Array<T,CLASS_NAME<T,T1,T2> >                  \
    {                                                                           \
        const ET_Array<T,T1> &m_t1;                                             \
        const ET_Array<T,T2> &m_t2;                                             \
    public:                                                                     \
        explicit CLASS_NAME(const ET_Array<T,T1> &t1, const ET_Array<T,T2> &t2) \
        : m_t1(t1), m_t2(t2) {}                                                 \
                                                                                \
       inline size_t expSize() const {return m_t1.expSize();}                   \
                                                                                \
       inline T evalAt(size_t index) const                                      \
       {                                                                        \
           return m_t1.evalAt(index) OP m_t2.evalAt(index);                     \
       }                                                                        \
    };                                                                          \
                                                                                \
    template<typename T, typename T1, typename T2>                              \
    inline CLASS_NAME<T,T1,T2>                                                  \
    operator OP (const ET_Array<T,T1> & t1, const ET_Array<T,T2> & t2)          \
    {                                                                           \
        return CLASS_NAME<T,T1,T2>(t1,t2);                                      \
    }                                                                           

MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR(+,ET_Sum)
MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR(-,ET_Diff)
MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR(*,ET_Mul)
MAKE_ET_SUPPORT_FOR_BINARY_OPERATOR(/,ET_Div)



//////////////////////////////////////////////////////////////////
/// \macro MAKE_ET_SUPPORT_FOR_FUNCTION
/// Generator for template derived expressions with
/// linked math function.Generated expression is like
/// unnary expression, which takes only one arg.-some list,
/// and then at "=" operator associated function by below
/// template classes evaluates function.
//////////////////////////////////////////////////////////////////
#define MAKE_ET_SUPPORT_FOR_FUNCTION(FUN_EXPR,FUN_NAME,FUN_DEF)     \
    template<typename T,typename C>                                 \
    class FUN_EXPR : public ET_Array<T,FUN_EXPR<T,C> >              \
    {                                                               \
        const ET_Array<T,C> & m_c;                                  \
    public:                                                         \
        explicit FUN_EXPR(const ET_Array<T,C> & c): m_c(c){}        \
                                                                    \
        inline size_t expSize() const {return m_c.expSize();}       \
                                                                    \
        inline T evalAt(size_t index) const                         \
        {                                                           \
            using namespace std;                                    \
            return FUN_DEF(m_c.evalAt(index));                      \
        }                                                           \
    };                                                              \
    template<typename T, typename C>                                \
    inline FUN_EXPR<T,C> FUN_NAME(const ET_Array<T,C> & c)          \
    {                                                               \
        return FUN_EXPR<T,C>(c);                                    \
    }                                                               

#include <cmath>
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Sin,  sin,  sin)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Cos,  cos,  cos)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Abs,  abs,  abs)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Exp,  exp,  exp)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Tan, tan, tan)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Atan, atan, atan)

MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Negate, operator-,-)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Positive, operator+,+)

MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Sinh,  sinh,  sinh)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Cosh,  cosh,  cosh)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Tanh,  tanh,  tanh)
MAKE_ET_SUPPORT_FOR_FUNCTION(ET_Atanh,  atanh, atanh)


//////////////////////////////////////////////////////////////////
/// \class ConversionTraits 
/// Class which provides compile time information
/// about result of operation between 2 base types
///------------------------------------------------
/// In case of "error::invalid use of incomplete type..." -this 
/// means that there is no known conversion between 2 provided 
/// types in equation - eg. Vector*Tensor--> there is no obvious
/// multiplication between those, so it must be done by some function
/// instead. 
//////////////////////////////////////////////////////////////////
template<typename T1,typename T2>
class ConversionTraits;

template<>
struct ConversionTraits<Scalar,Scalar>
{
    typedef Scalar return_type;
};
template<>
struct ConversionTraits<Scalar,Vector>
{
    typedef Vector return_type;
};
template<>
struct ConversionTraits<Vector,Scalar>
{
    typedef Vector return_type;
};
template<>
struct ConversionTraits<Vector,Vector>
{
    typedef Vector return_type;
};
template<>
struct ConversionTraits<Scalar,Tensor>
{
    typedef Tensor return_type;
};
template<>
struct ConversionTraits<Tensor,Scalar>
{
    typedef Tensor return_type;
};
template<>
struct ConversionTraits<Tensor,Tensor>
{
    typedef Tensor return_type;
};

//////////////////////////////////////////////////////////////////
/// \def macro SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS
/// Generates classes which defines operation between some kinde of
/// list expression with value type T2 and single value of type T1.
/// Returned expression list type is resolved from ConversionTraits
/// specialization. 
/// Below expressions solve issue when user would want to multiply
/// list of scalars with single scalar value, or list of scalar 
/// with single value vector varaible(returned expression would be
/// converted into vector type).
//////////////////////////////////////////////////////////////////
#define SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(CLASS_NAME,S_TYPE, OP)              \
template<typename T,typename ETDerived>                                                         \
class CLASS_NAME##LeftOperand                                                                   \
: public ET_Array<typename ConversionTraits<S_TYPE,T>::return_type,CLASS_NAME##LeftOperand<T,ETDerived> >    \
{                                                                                               \
    const S_TYPE m_singleVal;                                                                   \
    const ET_Array<T,ETDerived>& m_list;                                                        \
    typedef typename ConversionTraits<S_TYPE,T>::return_type return_type;                       \
public:                                                                                         \
    CLASS_NAME##LeftOperand(const S_TYPE& singleVal, const ET_Array<T,ETDerived> & list)        \
    : m_singleVal(singleVal), m_list(list)                                                      \
    {                                                                                           \
    }                                                                                           \
    size_t expSize() const { return m_list.expSize();}                                          \
    return_type evalAt(size_t i) const                                                          \
    {                                                                                           \
        return m_singleVal OP m_list.evalAt(i);                                                 \
    }                                                                                           \
    S_TYPE singleValue() const { m_singleVal;}                                                  \
    const ETDerived & listDerived() const                                                       \
    {                                                                                           \
        return static_cast<const ETDerived&>(m_list);                                           \
    }                                                                                           \
};                                                                                              \
template<typename T,typename ETDerived>                                                         \
class CLASS_NAME##RightOperand                                                                  \
: public ET_Array<typename ConversionTraits<S_TYPE,T>::return_type,CLASS_NAME##RightOperand<T,ETDerived> >    \
{                                                                                               \
    const S_TYPE m_singleVal;                                                                   \
    const ET_Array<T,ETDerived>& m_list;                                                        \
    typedef typename ConversionTraits<S_TYPE,T>::return_type return_type;                       \
public:                                                                                         \
    CLASS_NAME##RightOperand(const S_TYPE& singleVal, const ET_Array<T,ETDerived> & list)       \
    : m_singleVal(singleVal), m_list(list)                                                      \
    {                                                                                           \
    }                                                                                           \
    size_t expSize() const { return m_list.expSize();}                                          \
    return_type evalAt(size_t i) const                                                          \
    {                                                                                           \
        return m_list.evalAt(i) OP m_singleVal;                                                 \
    }                                                                                           \
    S_TYPE singleValue() const { m_singleVal;}                                                  \
    const ETDerived & listDerived() const                                                       \
    {                                                                                           \
    return static_cast<const ETDerived&>(m_list);                                               \
    }                                                                                           \
};                                                                                              \
template<typename T, typename ETDerived>                                                        \
CLASS_NAME##LeftOperand<T,ETDerived>                                                            \
operator OP                                                                                     \
(const S_TYPE & singleVal, const ET_Array<T,ETDerived> & list)                                  \
{                                                                                               \
    return CLASS_NAME##LeftOperand<T,ETDerived>(singleVal,list);                                \
}                                                                                               \
template<typename T, typename ETDerived>                                                        \
CLASS_NAME##RightOperand<T,ETDerived>                                                           \
operator OP                                                                                     \
(const ET_Array<T,ETDerived> & list, const S_TYPE & singleVal)                                  \
{                                                                                               \
    return CLASS_NAME##RightOperand<T,ETDerived>(singleVal,list);                               \
}                                                                                               


SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleScalarSum,Scalar,+)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleScalarSub,Scalar,-)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleScalarDiv,Scalar,/)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleScalarMul,Scalar,*)


SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleVectorSum,Vector,+)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleVectorSub,Vector,-)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleVectorDiv,Vector,/)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleVectorMul,Vector,*)

SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleTensorSum,Tensor,+)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleTensorSub,Tensor,-)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleTensorDiv,Tensor,/)
SEM_MAKE_ET_ARRAY_AND_SINGLE_VALUE_BINARY_OPERATORS(ETSingleTensorMul,Tensor,*)

// below is amibguous with scalar expression - in situation when  
// et_array inner type is Scalar and sing-val is also scalar -> it's 
// the same as ET_ScalarExpr. Due to fact, that any operation between
// some list and it's inner type can be not clear, so it is decided
// to remove below expression in faver of ET_ScalarExpr
// ///////////////////////////////////////////////////////////////////
// /// \class ET_SigValExpr
// /// \brief Expression representing single element value which
// /// will be at "=" operator processed for each
// /// element in list. This class allows to use
// /// single value (any T type) to be used with
// /// list supported by ET.
// /// \macro MAKE_ET_SUPPORT_FOR_SINGLE_VALUE
// /// \brief Generator for operators taking as lhs/rhs
// /// single value of T type (the same as list element type).
// /// This macro produce 2 operator functions which takes
// /// operatr OP(Expr, Value) or operatr OP( Value Expr)
// /// ----------------------------------------------------------
// /// This solution has one disadventage - ScalarExpr can't
// /// be used to directly assign to List, because ther will not
// /// exist conversion from valu(list elemnent type) to
// /// ET_ScalarExrp. eg. user can't write: List<double> v = 5.;
// /// To overcome this problem below, under this macro, there is
// /// introdueced special "operator =" which will do this conversion.
// /// (In fact this operator could be defined directly inside
// /// List, but to allow easier pluging in this ET solution to any
// /// type of list, it is introdueced here)
// //////////////////////////////////////////////////////////////////
// #define MAKE_ET_SUPPORT_FOR_SINGLE_VALUE(OP,NAME)               \
// template<typename T,typename OtherExp>                          \
// class NAME : public ET_Array<T, NAME<T,OtherExp> >              \
// {                                                               \
//     const OtherExp& m_exp;                                      \
//     T m_val;                                                    \
//     public:                                                     \
//         NAME(const T& val, const OtherExp& exp)                 \
//         :m_val(val),m_exp(exp){}                                \
//                                                                 \
//         inline size_t expSize() const {return m_exp.expSize();} \
//         inline T evalAt(size_t index) const                     \
//         {                                                       \
//             return m_exp.evalAt(index) OP m_val;                \
// }                                                               \
// };                                                              \
// template<typename T, typename T2>                               \
// inline NAME<T, ET_Array<T,T2> >                                 \
// operator OP (const T & t1, const ET_Array<T,T2> & t2)           \
// {                                                               \
//     return NAME<T, ET_Array<T,T2> >(t1,t2);                     \
// }                                                               \
// template<typename T, typename T2>                               \
// inline NAME<T, ET_Array<T,T2> >                                 \
// operator OP (const ET_Array<T,T2> & t2, const T & t1)           \
// {                                                               \
//     return NAME<T, ET_Array<T,T2> >(t1,t2);                     \
// }                                                               
// 
// MAKE_ET_SUPPORT_FOR_SINGLE_VALUE(+,ET_SigValSum)
// MAKE_ET_SUPPORT_FOR_SINGLE_VALUE(-,ET_SigValDiff)
// MAKE_ET_SUPPORT_FOR_SINGLE_VALUE(*,ET_SigValMul)
// MAKE_ET_SUPPORT_FOR_SINGLE_VALUE(/,ET_SigValDiv)

//////////////////////////////////////////////////////////////////
/// \class ET_ScalarExpr
/// \brief Expression representing scalar value which
/// will be at "=" operator processed for each
/// element in list. This class allows to use
/// scalar value to be used with
/// list supported by ET.
/// ----------------------------------------------------------
/// \macro MAKE_ET_SUPPORT_FOR_SCALAR_VALUE
/// \brief Generator for operators taking as lhs/rhs
/// single value of T type (the same as list element type).
/// This macro produce 2 operator functions which takes
/// operatr OP(Expr, Value) or operatr OP( Value Expr)
//////////////////////////////////////////////////////////////////
// #define MAKE_ET_SUPPORT_FOR_SCALAR_VALUE(OP,NAME)              \
// template<typename T,typename OtherExp>                         \
// class NAME : public ET_Array<T, NAME<T,OtherExp> >             \
// {                                                              \
//     const ET_Array<T,OtherExp>& m_exp;                         \
//     T m_val;                                                   \
// public:                                                        \
//     NAME(const T& val, const ET_Array<T,OtherExp>& exp)        \
//     :m_val(val),m_exp(exp){}                                   \
//                                                                \
//     inline size_t expSize() const {return m_exp.expSize();}    \
//     inline T evalAt(size_t index) const                        \
//     {                                                          \
//         return m_exp.evalAt(index) OP m_val;                   \
//     }                                                          \
//     T singValue() const { return m_val; }                      \
//     const OtherExp & arrayDerived() const                      \
//     {                                                          \
//         return static_cast<const OtherExp&>(m_exp);            \
//     }                                                          \
// };                                                             \
// template<typename T, typename T2>                              \
// inline NAME<T, T2>                                             \
// operator OP (const T & t1, const ET_Array<T,T2> & t2)          \
// {                                                              \
//     return NAME<T, T2>(t1,t2);                                 \
// }                                                              \
// template<typename T, typename T2>                              \
// inline NAME<T, T2>                                             \
// operator OP (const ET_Array<T,T2> & t2, const T & t1)          \
// {                                                              \
//     return NAME<T,T2>(t1,t2);                                  \
// }                                                              
// 
// MAKE_ET_SUPPORT_FOR_SCALAR_VALUE(+,ET_ScalarSum)
// MAKE_ET_SUPPORT_FOR_SCALAR_VALUE(-,ET_ScalarDiff)
// MAKE_ET_SUPPORT_FOR_SCALAR_VALUE(*,ET_ScalarMul)
// MAKE_ET_SUPPORT_FOR_SCALAR_VALUE(/,ET_ScalarDiv)


}//array
}//SEM




#endif // ET_Array_H
