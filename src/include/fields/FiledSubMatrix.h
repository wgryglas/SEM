#ifndef FILEDSUBMATRIX_H
#define FILEDSUBMATRIX_H

#include "utilities/Reference.h"
#include "utilities/TypeDefs.h"

namespace SEM { namespace field {

template<typename T, typename C>
class FieldSubMatrix
{
    REFERENCE_TYPE(FieldSubMatrix<T COMMA C>)

    C & m_data;
    Matrix<int>::type m_mapping;

public:
    FieldSubMatrix(C & data, const Matrix<int>::type & mapping)
        : m_data(data), m_mapping(mapping)
    {
    }

    FieldSubMatrix(const FieldSubMatrix& other)
        : m_data(other.m_data), m_mapping(other.m_mapping)
    {
    }


    int rows() const { return m_mapping.size(); }
    int cols() const { return m_mapping[0].size(); }

    T& operator ()(int i, int j)
    {
        return m_data[m_mapping[i][j]];
    }

    const T& operator ()(int i, int j) const
    {
        return m_data[m_mapping[i][j]];
    }

#define MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(OP)         \
    template<typename Matrix>                           \
    FieldSubMatrix<T,C> & operator OP (const Matrix& m) \
    {                                                   \
        for(int i=0; i<m_mapping.size(); ++i)           \
            for(int j=0; j<m_mapping[i].size(); ++j)    \
                m_data[m_mapping[i][j]] = m(i,j);       \
                                                        \
        return *this;                                   \
    }                                                   \

    MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(=)
    MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(+=)
    MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(-=)
    MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(*=)
    MAKE_FIELDSUBMATRIX_UNNARY_OPERATOR(/=)
};


//#define MAKE_FIELDSUBMATRIX_BINARY_OPERATOR(OP)                                     \
//template<typename T,typename C, typename Matrix>                                    \
//typename valMatrix<T>::type operator OP ( FieldSubMatrix<T,C> & lhs, Matrix & rhs)  \
//{                                                                                   \
//    typename valMatrix<T>::type retMat(lhs.rows(), lhs.cols());                     \
//    for(int i=0; i<lhs.rows(); ++i)                                                 \
//        for(int j=0; j<lhs.cols(); ++j)                                             \
//            retMat(i,j)=lhs(i,j) OP rhs(i,j);                                       \
//                                                                                    \
//    return retMat;                                                                  \
//}                                                                                   \
//template<typename T,typename C, typename Matrix>                                    \
//typename valMatrix<T>::type operator OP ( Matrix & lhs, FieldSubMatrix<T,C> & rhs)  \
//{                                                                                   \
//    typename valMatrix<T>::type retMat( rhs.rows(), rhs.cols() );                   \
//    for(int i=0; i<lhs.rows(); ++i)                                                 \
//        for(int j=0; j<lhs.cols(); ++j)                                             \
//            retMat(i,j)=lhs(i,j) OP rhs(i,j);                                       \
//                                                                                    \
//    return retMat;                                                                  \
//}                                                                                   \

//MAKE_FIELDSUBMATRIX_BINARY_OPERATOR(+)
//MAKE_FIELDSUBMATRIX_BINARY_OPERATOR(-)
//MAKE_FIELDSUBMATRIX_BINARY_OPERATOR(*)
//MAKE_FIELDSUBMATRIX_BINARY_OPERATOR(/)




}//field
}//SEM





#endif // FILEDSUBMATRIX_H
