#ifndef NUMERICALLISTWRAP_H
#define NUMERICALLISTWRAP_H

#include "utilities/ListWrap.h"

namespace SEM {

////////////////////////////////////////////////////////
/// Macro for generating operators on List
/// workaround due to use of expressions
/// in numerical arrays -->expressions are evaluated
/// in final stage of assigment-->operation is evaluated
/// in one loop, not for each binary operation
////////////////////////////////////////////////////////
#define MAKE_LIST_WRAP_BINARY_OPERATOR(LIST_WRAP, OP)                                  \
    template<typename T, template<class> class TypeProvider>                           \
    typename TypeProvider<T>::type operator OP (const LIST_WRAP<T,TypeProvider>& lhs,  \
                                      const LIST_WRAP<T,TypeProvider> & rhs)           \
    {                                                                                  \
        typename TypeProvider<T>::type result(lhs.size());                             \
        for(int i=0;i<lhs.size(); i++)                                                 \
        {                                                                              \
            result[i]=lhs[i] OP rhs[i];                                                \
        }                                                                              \
        return result;                                                                 \
    }                                                                                  \
    template<typename T, template<class> class TypeProvider>                           \
    typename TypeProvider<T>::type operator OP (const LIST_WRAP<T,TypeProvider>& lhs,  \
                                               const typename TypeProvider<T>::type & rhs) \
    {                                                                                      \
        typename  TypeProvider<T>::type result(lhs.size());                                \
        for(int i=0;i<lhs.size(); i++)                                                     \
        {                                                                                  \
            result[i]=lhs[i] OP rhs[i];                                                    \
        }                                                                                  \
        return result;                                                                     \
    }                                                                                      \
    template<typename T, template<class> class TypeProvider>                               \
    typename TypeProvider<T>::type operator OP (const typename TypeProvider<T>::type & lhs,\
                                               const LIST_WRAP<T,TypeProvider>& rhs)       \
    {                                                                                      \
        typename TypeProvider<T>::type result(lhs.size());                                 \
        for(int i=0;i<lhs.size(); i++)                                                     \
        {                                                                                  \
            result[i]=lhs[i] OP rhs[i];                                                    \
        }                                                                                  \
        return result;                                                                     \
    }                                                                                      \

//////////////////////////////////////////////////////
/// Macro to build unary operators as class
/// member variable
/////////////////////////////////////////////////////
#define MAKE_LIST_WRAP_UNARY_OPERATOR(LIST_WRAP, OP)                                                 \
    typename TypeProvider<T>::type & operator OP (const LIST_WRAP<T,TypeProvider> & other)           \
    {                                                                                                \
        for(int i=0; i<size();i++)                                                                   \
            m_data[i] OP other[i];                                                                   \
        return m_data;                                                                               \
    }                                                                                                \
    typename TypeProvider<T>::type & operator OP (const typename TypeProvider<T>::type & other)      \
    {                                                                                                \
        for(int i=0; i<size();i++)                                                                   \
            m_data[i] OP other[i];                                                                   \
        return m_data;                                                                               \
    }                                                                                                \

template <typename T, template<class> class TypeProvider>
class NumericalListWrap : public ListWrap<T,TypeProvider>
{

public:
    typedef ListWrap<T,TypeProvider> baseType;
    typedef typename baseType::Container Container;

    using baseType::m_data;
    using baseType::size;


    NumericalListWrap(size_t size=0):baseType(size){}
    NumericalListWrap(const NumericalListWrap<T,TypeProvider> &other):baseType(other){}
    NumericalListWrap(const Container & c): baseType(c){}

    virtual ~NumericalListWrap(){}

    NumericalListWrap<T,TypeProvider> & operator=(const NumericalListWrap<T,TypeProvider> & other)
    {
        m_data = other.m_data;
        return *this;
    }

    NumericalListWrap<T,TypeProvider> & operator=(const Container & data)
    {
        m_data = data;
        return *this;
    }

    MAKE_LIST_WRAP_UNARY_OPERATOR(NumericalListWrap,+=)
    MAKE_LIST_WRAP_UNARY_OPERATOR(NumericalListWrap,-)
    MAKE_LIST_WRAP_UNARY_OPERATOR(NumericalListWrap,*)
    MAKE_LIST_WRAP_UNARY_OPERATOR(NumericalListWrap,/)
};

MAKE_LIST_WRAP_BINARY_OPERATOR(NumericalListWrap,+)
MAKE_LIST_WRAP_BINARY_OPERATOR(NumericalListWrap,-)
MAKE_LIST_WRAP_BINARY_OPERATOR(NumericalListWrap,*)
MAKE_LIST_WRAP_BINARY_OPERATOR(NumericalListWrap,/)

} //SEM

#endif // NUMERICALLISTWRAP_H
