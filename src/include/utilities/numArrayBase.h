#ifndef NUMARRAYBASE_H
#define NUMARRAYBASE_H

#include "ET_Array.h"
#include "utilities/TypeDefs.h"
#include "numArraySlices.h"

namespace SEM {

//TODO ----> add iterator do numArrayBase, to be able to use
//           std algorithms in each numArrayBase, not only
//           in data storing implementation named numArray

/// \brief BaseTraits
///  provider for data types used by arrays.
///  Required due to problems with circular ref. of types in CRTP
///  when using typedefs for numArrayBase (typedefs was required due
///  to fact that std::vector with T=bool uses not bool& but special
///  memory treatment, and std::vector<bool>::reference is somthing
///  diffrent. This is why we coudn't use for operator [] directly
///  T&, but had to make typedef to std::vector<bool>::reference.
///
///  This is widely used solution with CRTP:
///  Declar baseTriats,
///  Define base class body,
///  Declare derived class,
///  define Traits specialization for  the derived class
///  Define derived class body
template<typename Derived>
class BaseTraits<numArrayBase<Derived> >
{
public:
    typedef typename BaseTraits< Derived >::value_type value_type;
    typedef typename BaseTraits< Derived >::reference reference;
    typedef typename BaseTraits< Derived >::const_reference const_reference;
};

template<typename Derived>
class BaseTraits<const numArrayBase<Derived> >
{
public:
    typedef typename BaseTraits< Derived >::value_type value_type;
    typedef typename BaseTraits< Derived >::const_reference reference;
    typedef typename BaseTraits< Derived >::const_reference const_reference;
};


/// \class numArrayBase
/// \brief Base interface to work with any array data directly
///        or indirectly by mappers extension of this class.
///        All classes deriving from this derives this interface,
///        so everything deriving from here can be used in the
///        same manner --> the same treamnt of array and any kind of
///        mappedArray or subArray.
template<typename DerivedData>
class numArrayBase : public array::ET_Array<typename BaseTraits<DerivedData>::value_type, numArrayBase<DerivedData> >
{
public:
    typedef numArrayBase<DerivedData> type;
    typedef typename BaseTraits<DerivedData>::value_type value_type;
    typedef typename BaseTraits<DerivedData>::reference reference;
    typedef typename BaseTraits<DerivedData>::const_reference const_reference;

    /// \brief slice - generate sub-vector from this vector
    /// \param start - start index
    /// \param end   - end index
    /// \return numArraySlice ->sub-vector manipulator for this class data
    typedef numArrayRangeMapped<DerivedData> rangeMapped;
    rangeMapped slice(const size_t &start, const size_t &end)
    {
        return rangeMapped( static_cast<DerivedData&>(*this), start,end );
    }

    typedef numArrayRangeMapped<const DerivedData> const_rangeMapped;
    const_rangeMapped slice(const size_t &start, const size_t &end) const
    {
        return const_rangeMapped( static_cast<const DerivedData&>(*this), start,end );
    }


    /// \brief slice - generate sub-vector from this vector
    /// \param mapping - collection of indexes refering to this class data
    /// \return numArraySlice ->sub-vector manipulator for this class data
    typedef numArrayIndexMapped<DerivedData> indexMapped;
    indexMapped slice(const VectorToVectorMap & mapping)
    {
        return indexMapped( static_cast<DerivedData&>(*this), mapping );
    }

    typedef numArrayIndexMapped<const DerivedData> const_indexMapped;
    const_indexMapped slice(const VectorToVectorMap & mapping) const
    {
        return const_indexMapped( static_cast<const DerivedData&>(*this), mapping );
    }

    /// \brief sliceArray2D - generate sub-array2d which refers to data stored in this array
    /// \param sizes - list of 2D indexes which defines mapping from this array data to 2D array
    /// \return sub-array2d which is manipulator of data in this class object.
    typedef numArray2DData1DIndexMapped<DerivedData> toArray2DMapped;
    toArray2DMapped sliceArray2D(const MatrixToVectorMap & mapping)
    {
        return toArray2DMapped( static_cast<DerivedData&>(*this), mapping );
    }

    typedef numArray2DData1DIndexMapped<const DerivedData> const_toArray2DMapped;
    const_toArray2DMapped sliceArray2D(const MatrixToVectorMap & mapping) const
    {
        return const_toArray2DMapped( static_cast<const DerivedData&>(*this), mapping );
    }


    /// \brief sliceArray2D - generate sub-array2d which refers to data stored in this array
    /// \param sizes - collection of subsequent rows column sizes
    /// \return sub-array2d which is manipulator of data in this class object.
    typedef numArray2DData1DSplited<DerivedData> toArray2DSplitted;
    toArray2DSplitted sliceArray2D(const std::vector<size_t> & sizes)
    {
        return toArray2DSplitted(static_cast<DerivedData&>(*this),sizes);
    }

    typedef numArray2DData1DSplited<const DerivedData> const_toArray2DSplitted;
    const_toArray2DSplitted sliceArray2D(const std::vector<size_t> & sizes) const
    {
        return const_toArray2DSplitted(static_cast<const DerivedData&>(*this),sizes);
    }

    //-----------------------------------------------------------------------//
    //                      numArrayBase INTERFACE                           //
    //-----------------------------------------------------------------------//
    inline reference operator [](const size_t &index) { return static_cast<DerivedData&>(*this)[index]; }
    inline const_reference operator [](const size_t &index) const { return static_cast<const DerivedData&>(*this)[index]; }
    inline size_t size() const { return static_cast<const DerivedData&>(*this).size(); }

    //------------------------------------------------------------------------//
    //                      ET_Array IMPLEMENTATION                           //
    //------------------------------------------------------------------------//
    inline value_type evalAt(const size_t &index) const {return static_cast<const DerivedData&>(*this)[index];}

    inline size_t expSize() const {return static_cast<const DerivedData&>(*this).size();}

    // ---------------------------------------------------------------------//
    //                      Expressions handling
    // ---------------------------------------------------------------------//
    /// single value assigment
#define MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(OP)                      \
    numArrayBase<DerivedData>& operator OP(const value_type &singVal)   \
    {                                                                   \
        for(int i=0; i< size(); ++i)                                    \
            (*this)[i] OP singVal;                                      \
                                                                        \
        return *this;                                                   \
    }                                                                   \
    
    MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(=)
    MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(-=)
    MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(+=)
    MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(*=)
    MAKE_NUMARRAYBASE_SINGLE_VAL_ASSIGMENT(/=)

    /// \macro MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR
    /// \brief To generate unary operators working with ET
    #define MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR(OP)                                                         \
        template<typename Exp>                                                                              \
        numArrayBase<DerivedData>& operator OP(const array::ET_Array<value_type,Exp> & expr)  \
        {                                                                                     \
            if(expr.expSize()!=size())                                          \
                ErrorInFunction                                                 \
                <<"Not consistient expression size "<<expr.expSize()            \
                <<"\n with numArrayBase size "<<size()                          \
                <<"\n in assigment operator"<<iomanagment::endProgram;          \
                                                                                \
            for(int i=0; i<size();++i)                                          \
            {                                                                   \
                 (*this)[i] OP expr.evalAt(i);                                  \
            }                                                                   \
            return *this;                                                       \
        }                                                                       \

    MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR( =)
    MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR(+=)
    MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR(-=)
    MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR(*=)
    MAKE_NUMARRAYBASE_ET_UNARY_OPERATOR(/=)

};



}//SEM




#endif // NUMARRAYBASE_H
