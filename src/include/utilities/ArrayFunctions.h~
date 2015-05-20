#ifndef ARRAYFUNCTIONS_H
#define ARRAYFUNCTIONS_H

#include "ET_Array.h"
#include "numArray2DBase.h"
#include "numArray.h"

namespace SEM { namespace array {
    
    template<typename T,typename Derived>
    T avg(const ET_Array<T,Derived>& list);
    
    template<typename T, typename Derived>
    T sum(const ET_Array<T,Derived>& list);
    
    template<typename T, typename Derived>
    T min(const ET_Array<T,Derived>& list);
    
    template<typename T, typename Derived>
    T max(const ET_Array<T,Derived>& list);
    
//     template<typename T, typename DerivedList, typename DerivedMatrix>
//     numArray<T> matrixMul(const numArray2DBase<DerivedMatrix> & matrix, const ET_Array<T,DerivedList> &vector )
//     {
//         numArray<T> result(matrix.size());
//         for(size_t r=0; r < matrix.size(); ++r)
//         {
//             result = 
//         }
//     }
    
    /// \class matrixMul - Expression representing matrix multiplication.
    /// This class when assigned to some kind of array converts automaticly
    /// into result of this operation
    template<typename DerivedArray, typename DerivedMatrix>
    class MatrixMultiplication : public  ET_Array<typename BaseTraits2D<DerivedMatrix>::inner_type,MatrixMultiplication<DerivedArray,DerivedMatrix> >
    {
        const numArrayBase<DerivedArray> &m_array;
        const numArray2DBase<DerivedMatrix> &m_matrix;
        typedef typename BaseTraits2D<DerivedMatrix>::inner_type T;
        
    public :
        MatrixMultiplication(const numArray2DBase<DerivedMatrix> &matrix,const numArrayBase<DerivedArray> &array):
            m_array(array), m_matrix(matrix)
        {
        }
        
        size_t expSize() const { return m_matrix.size();}
        T evalAt(size_t index) const { return sum(m_matrix.row(index)*m_array);}
    };
    
    template<typename DerivedArray, typename DerivedMatrix>
    MatrixMultiplication<DerivedArray, DerivedMatrix>
    matrixMul(const numArray2DBase<DerivedMatrix> &matrix,const numArrayBase<DerivedArray> &array)
    {
        return MatrixMultiplication<DerivedArray,DerivedMatrix>(matrix,array);
    }
    
    
    
    //------------------------------------------//
    //       implementations                    
    //------------------------------------------//
    template<typename T, typename Derived>
    T avg(const ET_Array<T,Derived>& list)
    {
        return sum(list)/list.expSize();
    }
    
    template<typename T, typename Derived>
    T sum(const ET_Array<T,Derived>& list)
    {
        T result=0;
        for(size_t i=0; i<list.expSize(); ++i)
        {
            result+=list.evalAt(i);
        }
        
        return result;
    }
    
    template<typename T, typename Derived>
    T min(const ET_Array<T,Derived>& list)
    {
        T result = list.evalAt(0);
        for(size_t i=1; i<list.expSize(); ++i)
            if(list.evalAt(i) < result)
                result = list.evalAt(i);
            
            return result;
    }
    
    template<typename T, typename Derived>
    T max(const ET_Array<T,Derived>& list)
    {
        T result = list.evalAt(0);
        for(size_t i=1; i<list.expSize(); ++i)
            if(list.evalAt(i) > result)
                result = list.evalAt(i);
            
            return result;
    }
    
    
    
}//array
}//SEM




#endif // ARRAYFUNCTIONS_H
