#ifndef numArray_H
#define numArray_H

#include <vector>
#include <iostream>


#include "ET_Array.h"
#include "iomanagment/ReadWriteStream.h"
#include "utilities/TypeDefs.h"
#include "numArrayBase.h"

namespace SEM {

///////////////////////////////////////////////////////////////////////////////
/// \class numArray
/// template container, desigend to work with any numerical valus.
/// Derives from numArrayBase, and implements special interface
/// and some operators for converting expr. into this class data.
/// ---------------------------------------------------------------------------
/// numArrayBase is base interface how any numArray can be used. Other
/// known derived classes from numArrayBase are simple mappers from some data,
/// which allow to use mapped data in the same way as oryginal data, but in mor
/// ---------------------------------------------------------------------------
/// Beacause it derives from ET_Array so it can be used as
/// just as any expresion, can handle any combination of +-*/ in
/// one single loop (less time spend on allocation/dealocation of
/// memory in temproal objects).
/// ---------------------------------------------------------------------------
/// This class derives also from std::vector, so it can be also
/// used as std::vector class;
/// ---------------------------------------------------------------------------
/// Because numArray is single parameter template, have operator [],
/// size and resize method, then can be easly used in
/// iomanagment::read/write method --> easy input output data.
///////////////////////////////////////////////////////////////////////////////
template<typename T>
class numArray;

template<typename T>
class BaseTraits<numArray<T> >
{
public:
    typedef T value_type;
    typedef typename std::vector<T>::reference reference;
    typedef typename std::vector<T>::const_reference const_reference;
};


template<typename T>
class numArray : public std::vector<typename BaseTraits<numArray<T> >::value_type>, public numArrayBase<numArray<T> >
{
    REFERENCE_TYPE(numArray<T>)
    typedef std::vector<T> dataBaseType;
    typedef numArrayBase< numArray<T> > numArrayBaseType;
public:
    typedef typename BaseTraits<numArray<T> >::value_type value_type;
    typedef typename BaseTraits<numArray<T> >::reference reference;
    typedef typename BaseTraits<numArray<T> >::const_reference const_reference;

    /// \brief default constructor
    explicit numArray(size_t size=0): dataBaseType(size)
    {
    }

    /// \brief numArray initializing constructor
    /// \param size size of array
    /// \param init initial value set for all elements
    numArray(size_t size, const value_type& init) : dataBaseType(size,init)
    {
    }

    /// \brief - copy constructor from other numArray
    /// (method is neccessery because expression assigment
    ///  do not allow resizement)
    numArray(const numArray<T> & other): dataBaseType(other)
    {
    }

    /// \brief - copy constructor from other base type of this class-general numArray
    /// (method is neccessery because expression assigment
    ///  do not allow resizement)
    template<typename OtherDerived>
    numArray(const numArrayBase<OtherDerived> &other): dataBaseType(other.size())
    {
        for(size_t i=0;i<other.size(); ++i)
            (*this)[i]=other[i];
    }

    virtual ~numArray(){}

    /// \brief - assigment operator from other vector
    /// (in fact this method is not neccessery because
    /// numArray is also ET_Array, so it is handled
    /// below, but this method would be a litle be faster)
    numArray<T> & operator = (const numArray<T>& other)
    {
        dataBaseType::operator =(other);
        return *this;
    }

    //-----------------------------------------------------------------------//
    //                   numArrayBase IMPLEMENTATION                         //
    //-----------------------------------------------------------------------//
    /// data base class implements automaticly numArrayBese interface
    using dataBaseType::operator [];
    using dataBaseType::size;
    using dataBaseType::resize;

    //------------------------------------------------------------------------//
    //                      EXPRESSION HANDLING
    //------------------------------------------------------------------------//
    /// \brief copy constructor from expression
    /// \param expr - any kinde of ET
    template<typename Expr>
    numArray(const array::ET_Array<value_type,Expr> & expr): dataBaseType(expr.expSize())
    {
        for(size_t i=0; i<expr.expSize(); ++i)
        {
            (*this)[i]=expr.evalAt(i);
        }
    }

    //-----------------------------------------------------------------------//
    //                   configure operators                                 //
    //-----------------------------------------------------------------------//
    /// \brief Include assigment operators form expression templates to
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, numArray<T>, numArrayBaseType )

//    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =,T,numArray<T>,SINGLE_ARG(numArrayBase<T,numArray<T> >))
//    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=,T,numArray<T>,SINGLE_ARG(numArrayBase<T,numArray<T> >))
//    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=,T,numArray<T>,SINGLE_ARG(numArrayBase<T,numArray<T> >))
//    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=,T,numArray<T>,SINGLE_ARG(numArrayBase<T,numArray<T> >))
//    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=,T,numArray<T>,SINGLE_ARG(numArrayBase<T,numArray<T> >))
};




//template<typename T>
//class numArray : public std::vector<T>, public numArrayBase<T,numArray<T> > //public ET_Array<T,numArray<T> >
//{
//    REFERENCE_TYPE(numArray<T>)

//    typedef std::vector<T> baseVector;

//public:

//    /// \brief default constructor
//    numArray(size_t size=0): baseVector(size)
//    {
//    }

//    /// \brief - copy constructor from other vector
//    /// (in fact this method is not neccessery because
//    /// numArray is also ET_Array, so it is handled
//    /// below, but this method would be a litle be faster)
//    numArray(const numArray<T> & other): baseVector(other)
//    {
//    }

//    /// \brief - assigment operator from other vector
//    /// (in fact this method is not neccessery because
//    /// numArray is also ET_Array, so it is handled
//    /// below, but this method would be a litle be faster)
//    numArray<T> & operator = (const numArray<T>& other)
//    {
//        baseVector::operator =(other);
//        return *this;
//    }

//    /// \brief operator = allow assigment one value to all elements
//    ///                   (it apears here, because there is no solution
//    ///                    to do that as expr. template)
//    /// \param singleValue -value to be assigned to all elements
//    /// \return ref. to this object
//    numArray<T> & operator = (const T& singleValue)
//    {
//        for(int i=0; i<baseVector::size(); ++i) (*this)[i]=singleValue;
//        return *this;
//    }

//    /// \brief slice - generate sub-vector from this vector
//    /// \param start - start index
//    /// \param end   - end index
//    /// \return numArraySlice ->sub-vector manipulator for this class data
//    numArraySlice<T,numArrayRangeMapper> slice(size_t start, size_t end)
//    {
//        return numArraySlice<T,numArrayRangeMapper>( numArrayRangeMapper<T>(*this, start, end) );
//    }

//    /// \brief slice - generate sub-vector from this vector
//    /// \param mapping - collection of indexes refering to this class data
//    /// \return numArraySlice ->sub-vector manipulator for this class data
//    numArraySlice<T,numArrayIndexMapper> slice(const VectorToVectorMap& mapping)
//    {
//        return numArraySlice<T,numArrayIndexMapper>( numArrayIndexMapper<T>(*this, mapping) );
//    }

//    //------------------------------------------------------------------------//
//    //                      numArrayBase implementation                       //
//    //------------------------------------------------------------------------//
//    /// \brief whole numArrayBase interface is implemented automaticly
//    ///        by std::vector
//    using baseVector::operator[];
//    using baseVector::size;

//    //------------------------------------------------------------------------//
//    //                      EXPRESSION HANDLING
//    //------------------------------------------------------------------------//
//    /// \brief copy constructor from expression
//    /// \param expr - any kinde of ET
//    template<typename Expr>
//    numArray(const ET_Array<T,Expr> & expr): baseVector(expr.expSize())
//    {
//        for(size_t i=0; i<expr.expSize(); ++i)
//        {
//            (*this)[i]=expr.evalAt(i);
//        }
//    }

//    /// \macro MAKE_NUM_EXPR_UNARY_OP
//    /// Generator for unary operators function which will meke
//    /// assigment of expression evaluation into this container
//    #define MAKE_NUM_EXPR_UNARY_OP(OP)                                      \
//        template<typename Expr>                                             \
//        numArray<T> & operator OP(const ET_Array<T,Expr> & expr)            \
//        {                                                                   \
//            using namespace iomanagment;                                    \
//            if(baseVector::size()!=expr.expSize())                          \
//            ErrorInFunction<<"incoherent size of expression with numArray\n"\
//                           <<"while evaluating expression into numArray"    \
//                           <<endProgram;                                    \
//                                                                            \
//            for(size_t i=0; i<baseVector::size(); ++i)                      \
//            {                                                               \
//                (*this)[i] OP expr.evalAt(i);                               \
//            }                                                               \
//            return *this;                                                   \
//        }                                                                   \

//    MAKE_NUM_EXPR_UNARY_OP(=)
//    MAKE_NUM_EXPR_UNARY_OP(+=)
//    MAKE_NUM_EXPR_UNARY_OP(-=)
//    MAKE_NUM_EXPR_UNARY_OP(*=)

//    /// \brief ~numArray -virtual destructor for this container
//    virtual ~numArray(){}

//};






}//SEM




#endif // numArray_H
