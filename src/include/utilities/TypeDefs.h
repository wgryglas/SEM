#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#include <vector>

#include "boost/array.hpp"

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Core"

#include "components/Basic.h"
#include "Reference.h"

namespace SEM {

/// \typedef basic type for index
///(if space will be bigger here wi can switch it to long int)
typedef int Index;

/// Workaround typedef for defining new
/// template type from more complex template
/// (in c++11 there is direct utility for such a typedefs)
template <typename T>
struct Matrix
{
    typedef std::vector<std::vector<T> > type;
    REFERENCE_TYPE(type)
};

template<typename T>
struct Matrix4
{
    typedef std::vector<std::vector<std::vector<std::vector<T> > > > type;
    REFERENCE_TYPE(type)
};

#define ALLOW_EIGEN_EXPRESION_TEMPLATES_IN_DERIVED_CALSS(CLASS_NAME, CLASS_TYPE, BASE_TYPE)            \
    template<typename OtherDerived>                                                                    \
    CLASS_NAME(const Eigen::MatrixBase<OtherDerived> & other):BASE_TYPE(other)                         \
    {                                                                                                  \
    }                                                                                                  \
    template<typename OtherDerived>                                                                    \
    CLASS_TYPE & operator =(const Eigen::MatrixBase<OtherDerived> & other)                             \
    {                                                                                                  \
        BASE_TYPE::operator=(other);                                                                   \
        return *this;                                                                                  \
    }                                                                                                  \

/// \typedef mapping from matrix into matrix
typedef Matrix<boost::array<Index,2> >::type MatrixToMatrixMap;

/// \typedef mapping from vector into matrix
typedef std::vector< boost::array<Index,2> > VectorToMatrixMap;

/// \typedef mapping form matrix to vector
typedef Matrix<Index>::type MatrixToVectorMap;

/// \typedf mapping from one vector into another
typedef std::vector<Index> VectorToVectorMap;


//struct EquationMatrix
//{
//    typedef Eigen::SparseMatrix<Scalar> type;
//    REFERENCE_TYPE(Eigen::SparseMatrix<Scalar>)
//};

//struct EquationVector
//{
//    typedef Eigen::Matrix<Scalar,Eigen::Dynamic,1> type;
//    REFERENCE_TYPE(type)
//};

struct EMCoefficients
{
    typedef std::vector<Eigen::Triplet<Scalar> > type;
    REFERENCE_TYPE(type)
};

struct EMCoefff
{
    typedef Eigen::Triplet<Scalar> type;
    REFERENCE_TYPE(type)
};


///////////////////////////////////////////////////////////
/// Predefined containers type providers to use
/// inside class like List<T,TypeProvider>
///////////////////////////////////////////////////////////
#define MAKE_TYPE_PROVIDER(NAME, TYPE_EXPRESSED_IN_T)                      \
    template<typename T> struct NAME { typedef TYPE_EXPRESSED_IN_T type;}; \


MAKE_TYPE_PROVIDER(stdVectorType, std::vector<T>)
MAKE_TYPE_PROVIDER(stdListType, std::list<T>)



}//SEM

#endif // TYPEDEFS_H
