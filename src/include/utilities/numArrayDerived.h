#ifndef NUMARRAYDERIVED_H
#define NUMARRAYDERIVED_H

#include "TypeDefs.h"

/// \macro ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE
/// \brief - mcro for easy generating possibility of use =,+=,.. witch ET_Array
///          in any derived from numArrayBase class
/// \param OP - operator sign
/// \param T_NAME - how data template type is called in your class definiton
/// \param CLASS_TYPE - full class type - including all templates arguments within <...>
///                     (if more then 1, use SINGLE_ARG macro to pass arg. with commas)
/// \param BASE_ARRAY_TYPE - some base class from what your class is deriving
///                         (not that BASE_ARRAY_TYPE must derive from numArrayBas, and
///                          all classes in this type hierarchy must also use this
///                          macro)
#define ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( OP, T_NAME, CLASS_TYPE,BASE_ARRAY_TYPE)\
        template<typename Exp>                                                      \
        CLASS_TYPE & operator OP (const SEM::array::ET_Array<T_NAME,Exp> &expr)     \
        {                                                                           \
            BASE_ARRAY_TYPE::operator OP(expr);                                     \
            return *this;                                                           \
        }                                                                 
/// \macro ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE
/// \brief generator assigment from single value to all elements in
///        class derived from numArrayBase.
/// \param T_NAME - how data template type is called in your class definiton
/// \param CLASS_TYPE - full class type - including all templates arguments within <...>
///                     (if more then 1, use SINGLE_ARG macro to pass arg. with commas)
/// \param BASE_ARRAY_TYPE - some base class from what your class is deriving
///                         (not that BASE_ARRAY_TYPE must derive from numArrayBas, and
///                          all classes in this type hierarchy must also use this
///                          macro)
#define ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(OP,T_NAME, CLASS_TYPE,BASE_ARRAY_TYPE)\
    CLASS_TYPE & operator OP(const T_NAME &singleVal)                      \
    {                                                                      \
        BASE_ARRAY_TYPE::operator OP(singleVal);                           \
        return *this;                                                      \
    }                                                                      

/// \macro CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE
/// \brief macro that shall configure all necessery assigment operators in
///        derived class from numArrayBase
/// \param OP - operator sign
/// \param T_NAME - how data template type is called in your class definiton
/// \param CLASS_TYPE - full class type - including all templates arguments within <...>
///                     (if more then 1, use SINGLE_ARG macro to pass arg. with commas)
/// \param BASE_ARRAY_TYPE - some base class from what your class is deriving
///                         (not that BASE_ARRAY_TYPE must derive from numArrayBas, and
///                          all classes in this type hierarchy must also use this
///                          macro)
#define CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(T_NAME, CLASS_TYPE, BASE_ARRAY_TYPE)                         \
    ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )  \
    ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(+=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) ) \
    ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(-=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) ) \
    ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(*=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) ) \
    ALLOW_ASSIGMENT_FROM_SINGLE_VAL_TO_DERIVED_FROM_NUMARRAYBASE(/=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) ) \
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )         \
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )         \
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )         \
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )         \
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=,T_NAME,SINGLE_ARG(CLASS_TYPE),SINGLE_ARG(BASE_ARRAY_TYPE) )


#endif // NUMARRAYDERIVED_H
