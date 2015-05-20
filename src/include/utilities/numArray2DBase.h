#ifndef NUMARRAY2DBASE_H
#define NUMARRAY2DBASE_H

#include <stdlib.h>

#include "ET_Array.h"
#include "TypeDefs.h"
#include "numArrayBase.h"
#include "numArraySlices.h"

namespace SEM {

/// \class numArray2DBase
/// \brief Base interface to work with any array 2D data directly
///        or indirectly by mappers - extension of this class.
///        All classes deriving from this derives this interface,
///        so everything deriving from here can be used in the
///        same manner --> the same treamnt of array2D and any kind of
///        mappedArray2D. Real data storage shall be implemented as
///        derived from this class as well(example should be numArray2D)

template<typename Derived>
class BaseTraits2D;

template<typename DerivedData>
class numArray2DBase : public array::ET_Array<typename BaseTraits2D<DerivedData>::inner_type, numArray2DBase<DerivedData> >
{

public:
    typedef numArray2DBase<DerivedData> type;
    typedef typename BaseTraits2D<DerivedData>::inner_type inner_type;
    typedef typename BaseTraits2D<DerivedData>::inner_reference inner_reference;
    typedef typename BaseTraits2D<DerivedData>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<DerivedData>::value_type value_type;
    typedef typename BaseTraits2D<DerivedData>::reference reference;
    typedef typename BaseTraits2D<DerivedData>::const_reference const_reference;

    //--------------------------------------------------//
    //          numArray2DBase interface                //
    //--------------------------------------------------//
    /// \brief operator [] - inner data type getter
    /// \param index       - 1st dimension index
    /// \return  - somthing like row from 2D data
    inline reference operator[](size_t& index)
    {
        return static_cast<DerivedData&>(*this)[index];
    }

    /// \brief operator [] - inner data type const getter
    /// \param index       - 1st dimension index
    /// \return  - somthing like row from 2D data with const specification type
    inline const_reference operator[](const size_t& index) const
    {
        return static_cast<const DerivedData&>(*this)[index];
    }

    /// \brief size - getter for 1st dimension size
    /// \return  size of first dimension, eg. number of rows in 2D data
    inline size_t size() const { return static_cast<const DerivedData*>(this)->size();}

    typedef value_type rowType;
    inline reference row(const size_t &index) { return static_cast<DerivedData&>(*this)[index];}

    typedef const value_type const_rowType;
    inline const_reference row(const size_t &index) const { return static_cast<const DerivedData&>(*this)[index];}


    typedef numArrayData2DMapped<DerivedData> columnType;
    inline columnType column(const size_t &index)
    {
        using namespace iomanagment;

        for(size_t i=0; i<size(); ++i)
            if( row(i).size() <= index )
            {
                ErrorInFunction<<"Can't obtain column from 2D array, because "
                               <<index<<" don't cut all elements in non-squared 2D array"
                               <<endProgram;
            }

        VectorToMatrixMap map(size());

        for(size_t i=0; i<size(); ++i)
        {
            map[i][0]=i;
            map[i][1]=index;
        }

        return columnType(static_cast<DerivedData&>(*this),map);
    }

    typedef numArrayData2DMapped<const DerivedData> const_columnType;
    inline const_columnType column(const size_t &index) const
    {
        using namespace iomanagment;

        for(size_t i=0; i<size(); ++i)
            if( row(i).size() <= index )
            {
                ErrorInFunction<<"Can't obtain column from 2D array, because "
                               <<index<<" don't cut all elements in non-squared 2D array"
                               <<endProgram;
            }

        VectorToMatrixMap map(size());

        for(size_t i=0; i<size(); ++i)
        {
            map[i][0]=i;
            map[i][1]=index;
        }

        return const_columnType(static_cast<const DerivedData&>(*this),map);
    }


    typedef numArrayDiagonalFrom2D<DerivedData> diagonalType;
    inline diagonalType diagonal()
    {
        return diagonalType(static_cast<DerivedData&>(*this));
    }

    typedef numArrayDiagonalFrom2D<const DerivedData> const_diagonalType;
    inline const_diagonalType diagonal() const
    {
        return const_diagonalType(static_cast<const DerivedData&>(*this));
    }


    /// \brief slice method to get sub-arra2d from this object.
    ///        sub-array2D defined by range
    /// \param start1D - start 1D index
    /// \param end1D   - end 1D index
    /// \param start2D - start 2D index
    /// \param end2D   - end 2D index
    /// \return "sub-numArray2D" type
    typedef numArray2DRangeMapped<DerivedData> rangeMapped;
    inline rangeMapped slice(const size_t &start1D,const size_t &end1D,const size_t &start2D,const size_t &end2D)
    {
        return rangeMapped(static_cast<DerivedData&>(*this),start1D,end1D,start2D,end2D);
    }

    typedef numArray2DRangeMapped<const DerivedData> const_rangeMapped;
    inline const_rangeMapped slice(const size_t &start1D,const size_t &end1D,const size_t &start2D,const size_t &end2D) const
    {
        return const_rangeMapped(static_cast<const DerivedData&>(*this),start1D,end1D,start2D,end2D);
    }

    /// \brief slice method to get sub-arra2d from this object.
    ///        sub-array2D defined by mapping from array2D to array2D.
    /// \param map - mapping 2D array of pair values(indexes i,j from this data object)
    /// \return "sub-numArray2D" type
    typedef numArray2DIndexMapped<DerivedData> indexMapped ;
    inline numArray2DIndexMapped<DerivedData> slice(const MatrixToMatrixMap &map)
    {
        return numArray2DIndexMapped<DerivedData>(static_cast<DerivedData&>(*this),map);
    }

    typedef numArray2DIndexMapped<DerivedData> const_indexMapped;
    inline const_indexMapped slice(const MatrixToMatrixMap &map) const
    {
        return const_indexMapped(static_cast<const DerivedData&>(*this),map);
    }

    /// \brief sliceArray - make 1D array, which elements referes to data in this 2D array
    /// \param map        - mapping from 2D data to 1D
    /// \return           - kind of reference array, here called slice
    typedef numArrayData2DMapped<DerivedData> mappedArray;
    inline mappedArray sliceArray(const VectorToMatrixMap &map)
    {
        return mappedArray(static_cast<DerivedData&>(*this),map);
    }

    typedef numArrayData2DMapped<const DerivedData> const_mappedArray;
    inline const_mappedArray sliceArray(const VectorToMatrixMap &map) const
    {
        return const_mappedArray(static_cast<const DerivedData&>(*this),map);
    }


    //-----------------------------------------------------//
    //              ET_Array implementation                //
    //-----------------------------------------------------//
    inline inner_type evalAt(const size_t &index) const
    {
        using namespace iomanagment;

        size_t curr=0;
        size_t i=0;
       
        while(index>=curr+(*this)[i].size() )
        {
            if( i == size() )//static_cast<const DerivedData*>(this)->size()
            {
                ErrorInFunction<<"Index out of size, while getting value from numArray2D \n"
                               <<"for expression evaluation"<<endProgram;
            }
            curr+=(*this)[i].size();
            ++i;
        }

        return static_cast<const DerivedData&>(*this)[i][index-curr];
    }

    inline size_t expSize() const
    {
        size_t s=0;
        for(size_t i=0; i<size(); ++i)
        {
            s+=(*this)[i].size();
        }
        return s;
    }

    //-----------------------------------------------------//
    //             EXPRESSIONS HANDLING                    //
    //-----------------------------------------------------//
    /// \brief copy constructor from expression -not allowd
    ///     -->no inforamtion about second dimmensions

    /// \brief operator = - assigment form single value to all elements
    ///                     which can't be directly done by ET_Array library
    /// \param singVal    - value to assign
    /// \return           - this object
    numArray2DBase<DerivedData> & operator =(const inner_type &singVal)
    {
        for(size_t i=0; i<size(); ++i)
            for(size_t j=0; j<(*this)[i].size(); ++j)
            {
                (*this)[i][j]=singVal;
            }

        return *this;
    }

    /// \macro MAKE_NUMARRAY2D_EXPR_UNARY_OP
    /// Generator for unary operators function which will meke
    /// assigment of expression evaluation into this container
    #define MAKE_NUMARRAY2D_EXPR_UNARY_OP(OP)                                                   \
        template<typename Expr>                                                                 \
        numArray2DBase<DerivedData> & operator OP(const array::ET_Array<inner_type,Expr> & expr)\
        {                                                                                       \
            using namespace iomanagment;                                                        \
            if(expr.expSize()!=expSize())                                                       \
            {                                                                                   \
                ErrorInFunction<<"incoherent size between numArray2D("<<expSize()<<")"          \
                               <<"and exprssion("<<expr.expSize()<<")"<<endProgram;             \
            }                                                                                   \
            size_t curr=0;                                                                      \
            for(size_t i=0; i<size(); ++i)                                                      \
            {                                                                                   \
                for(size_t j=0; j<(*this)[i].size(); ++j)                                       \
                {                                                                               \
                    (*this)[i][j] OP expr.evalAt(curr);                                         \
                    ++curr;                                                                     \
                }                                                                               \
            }                                                                                   \
            return *this;                                                                       \
        }                                                                                       \

    MAKE_NUMARRAY2D_EXPR_UNARY_OP( =)
    MAKE_NUMARRAY2D_EXPR_UNARY_OP(+=)
    MAKE_NUMARRAY2D_EXPR_UNARY_OP(-=)
    MAKE_NUMARRAY2D_EXPR_UNARY_OP(*=)

};
}//SEM

#endif // NUMARRAY2DBASE_H
