#ifndef VECTOR_UTILS_H
#define VECTOR_UTILS_H

#include <stdlib.h>

#include "components/Basic.h"
#include "utilities/Utilities.h"
#include "utilities/ET_Array.h"
#include "utilities/numArrayBase.h"
#include "utilities/numArraySlices.h"
#include "utilities/numArray.h"
#include "iomanagment/InfoStream.h"


namespace SEM {
    
template<typename T, size_t Dim>
struct BaseTraits<VectorTX<T,Dim> >
{
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
};

template<typename T, size_t Dim>
struct BaseTraits<const VectorTX<T,Dim> >
{
    typedef T value_type;
    typedef const T& reference;
    typedef const T& const_reference;
};


template<typename ArrayType>
class numArrayVectorComponent;

template<typename ArrayType>
struct BaseTraits<numArrayVectorComponent<ArrayType> >
{
    /// typedefs to extract type T held by VectorXT<T,Dim> element, which is stored in ArrayType by asumption
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::value_type value_type;
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::reference reference;
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::const_reference const_reference;
};

template<typename ArrayType>
struct BaseTraits<numArrayVectorComponent<const ArrayType> >
{
    /// typedefs to extract type T held by VectorXT<T,Dim> element, which is stored in ArrayType by asumption
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::value_type value_type;
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::const_reference reference;
    typedef typename BaseTraits<typename BaseTraits<ArrayType>::value_type>::const_reference const_reference;
};

/// \class numArrayVectorComponent
/// \brief class which behaves as single list to operate on Vector element component.
/// Shall be used when user have some kinde of numArray<Vector> (or its slice) and 
/// want to operate on one component of Vector (eg. Vector.x() ) like with list consiting
/// of this component elements
template<typename ArrayType>
class numArrayVectorComponent : public numArrayBase<numArrayVectorComponent<ArrayType> >
{
    typedef numArrayVectorComponent<ArrayType> type;
    typedef numArrayBase<type> baseType;
    
    ArrayType & m_refData;
    size_t m_dim;
    
public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;
    
    numArrayVectorComponent(ArrayType &data, const int& dim)
    : m_refData(data),m_dim(dim)
    {
    }
    
    numArrayVectorComponent(const type &other)
    : m_refData(other.m_refData),m_dim(other.m_dim)
    {
    }
    
    //-----------------------------------------------------------------------//
    //                   numArrayBase IMPLEMENTATION                         //
    //-----------------------------------------------------------------------//
    inline const_reference operator[](const size_t &loc) const
    {
        if(loc >= m_refData.size())
            ErrorInFunction<<"Index out of range in numArrayVectorComponent\n"
            <<"asked for "<<loc<<" element and size is"<<m_refData.size()<<"\n"
            <<"-class mapping 1D data to another 1D data"
            <<iomanagment::endProgram;
        
        return m_refData[loc][m_dim];
    }
    
    inline reference operator[](const size_t &loc)
    {
        if(loc >= m_refData.size())
            ErrorInFunction<<"Index out of range in numArrayVectorComponent\n"
            <<"asked for "<<loc<<" element and size is"<<m_refData.size()<<"\n"
            <<"-class mapping 1D data to another 1D data"
            <<iomanagment::endProgram;
        
        return m_refData[loc][m_dim];
    }
    
    inline size_t size() const {return m_refData.size();}
    //-----------------------------------------------------------------------//
    //                   configure operators                                 //
    //-----------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, type, baseType )
};

template<typename DerivedArray>
numArrayVectorComponent<numArrayBase<DerivedArray> > xCmps(numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<numArrayBase<DerivedArray> >(vecArray,0);
}
template<typename DerivedArray>
numArrayVectorComponent<const numArrayBase<DerivedArray> > xCmps(const numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<const numArrayBase<DerivedArray> >(vecArray,0);
}

template<typename DerivedArray>
numArrayVectorComponent<numArrayBase<DerivedArray> > yCmps(numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<numArrayBase<DerivedArray> >(vecArray,1);
}
template<typename DerivedArray>
numArrayVectorComponent<const numArrayBase<DerivedArray> > yCmps(const numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<const numArrayBase<DerivedArray> >(vecArray,1);
}

template<typename DerivedArray>
numArrayVectorComponent<numArrayBase<DerivedArray> > zCmps(numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<numArrayBase<DerivedArray> >(vecArray,2);
}
template<typename DerivedArray>
numArrayVectorComponent<const numArrayBase<DerivedArray> > zCmps(const numArrayBase<DerivedArray> & vecArray)
{
    return numArrayVectorComponent<const numArrayBase<DerivedArray> >(vecArray,2);
}



//----------------------------------------------------------------------------//
// functions to operate on any array/expression with Vector type            
//----------------------------------------------------------------------------//
namespace array {

// This issue is solved by new apporach - all possible singleValue-list operations
// are defined by 1 template class now :)
    
// --------------- scalar - vector list  operations -------------------------------- //
// / \def MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING 
// / macro to generate operation handling between single scalar value and vector list
// / macro generates class defining operation and 2 operator functions
// #define MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING(CLASS_NAME, OP)                          \
// template<typename T, size_t Dim, typename Array>                                        \
// struct CLASS_NAME : public ET_Array<VectorTX<T,Dim>,CLASS_NAME<T,Dim,Array> >           \
// {                                                                                       \
//     CLASS_NAME(const Scalar & scalarVal, const ET_Array<VectorTX<T,Dim>,Array> &list)   \
//     : m_scalarValue(scalarVal), m_list(list)                                            \
//     {                                                                                   \
//     }                                                                                   \
//     size_t expSize() const { return m_list.expSize();}                                  \
//     VectorTX<T,Dim> evalAt(size_t pos) const                                            \
//     {                                                                                   \
//         return m_scalarValue OP m_list.evalAt(pos);                                     \
//     }                                                                                   \
// private:                                                                                \
//     Scalar m_scalarValue;                                                               \
//     const ET_Array<VectorTX<T,Dim>,Array> & m_list;                                     \
// };                                                                                      \
// template<typename T, size_t Dim, typename Array>                                        \
// CLASS_NAME<T,Dim,Array> operator OP                                                     \
// (const Scalar & scalarVal, const ET_Array<VectorTX<T,Dim>,Array> &list)                 \
// {                                                                                       \
//     return CLASS_NAME<T,Dim,Array>(scalarVal,list);                                     \
// }                                                                                       \
// template<typename T, size_t Dim, typename Array>                                        \
// CLASS_NAME<T,Dim,Array> operator OP                                                     \
// (const ET_Array<VectorTX<T,Dim>,Array> &list,const Scalar & scalarVal)                  \
// {                                                                                       \
//     return CLASS_NAME<T,Dim,Array>(scalarVal,list);                                     \
// }                                                                                       \
// 
// MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING(ScalarVectorListAdd,  +)
// MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING(ScalarVectorListSub,  -)
// MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING(ScalarVectorListMult, *)
// MAKE_VECTOR_LIST_SCALAR_VALUE_HANDLING(ScalarVectorListDiv,  /)
    
    
//----------------- dot product ---------------------------------------------//    
template<typename T, size_t Dim, typename Lhs, typename Rhs>
struct DotProdExpression : public ET_Array<T,DotProdExpression<T,Dim,Lhs,Rhs> >
{
    DotProdExpression(const ET_Array<VectorTX<T,Dim>,Lhs> &lhs,const ET_Array<VectorTX<T,Dim>,Rhs> &rhs)
    : m_lhs(lhs),m_rhs(rhs)
    {
        if(lhs.expSize() != rhs.expSize())
        {
            ErrorInFunction<<" Inconsistient sizes between 2 expressions "
                            <<" provided to calculate dotProd on them"
                            <<" first expr size = "<<lhs.expSize() 
                            <<" second expr size ="<<rhs.expSize() 
                            <<iomanagment::endProgram;
        }
    }
    
    T evalAt(size_t pos) const 
    { 
        return dotProd(m_lhs.evalAt(pos), m_rhs.evalAt(pos));
    }
    
    size_t expSize() const { return m_lhs.expSize();}
    
private:
    const ET_Array<VectorTX<T,Dim>,Lhs> &m_lhs;
    const ET_Array<VectorTX<T,Dim>,Rhs> &m_rhs;
};
    
template<typename T, size_t Dim,typename Derived1, typename Derived2>
DotProdExpression<T,Dim,Derived1,Derived2> dotProd(const ET_Array<VectorTX<T,Dim>, Derived1>& list1, const ET_Array<VectorTX<T,Dim>, Derived2>& list2)
{
    return DotProdExpression<T,Dim,Derived1,Derived2>(list1,list2);
}

template<typename T,size_t Dim, typename Derived>
struct DotProdSingleValExpression : public ET_Array<T,DotProdSingleValExpression<T,Dim,Derived> >
{
    DotProdSingleValExpression(const ET_Array<VectorTX<T,Dim>,Derived> &list, const VectorTX<T,Dim> & singVal)
    : m_list(list), m_singVal(singVal)
    {
    }
    
    T evalAt(size_t pos) const 
    { 
        return dotProd(m_list.evalAt(pos), m_singVal);
    }
    
    size_t expSize() const { return m_list.expSize();}
        
private:
    VectorTX<T,Dim> m_singVal;
    const ET_Array<VectorTX<T,Dim>,Derived> &m_list;
};

template<typename T,size_t Dim, typename Derived>
DotProdSingleValExpression<T,Dim,Derived> dotProd(const ET_Array<VectorTX<T,Dim>,Derived> &list, const VectorTX<T,Dim> & singVal)
{
    return DotProdSingleValExpression<T,Dim,Derived>(list,singVal);
}

template<typename T,size_t Dim, typename Derived>
DotProdSingleValExpression<T,Dim,Derived> dotProd(const VectorTX<T,Dim> & singVal,const ET_Array<VectorTX<T,Dim>,Derived> &list)
{
    return DotProdSingleValExpression<T,Dim,Derived>(list,singVal);
}


//----------------- cross product 2d ---------------------------------------------//
template<typename T, typename Lhs, typename Rhs>
struct CrossProd2DExpression : public ET_Array<T,CrossProd2DExpression<T,Lhs,Rhs> >
{
    CrossProd2DExpression(const ET_Array<VectorTX<T,2>,Lhs> &lhs,const ET_Array<VectorTX<T,2>,Rhs> &rhs)
    : m_lhs(lhs), m_rhs(rhs)
    {
        if(lhs.expSize() != rhs.expSize())
        {
            ErrorInFunction<<" Inconsistient sizes between 2 expressions "
            <<" provided to calculate crossProd on them"
            <<" first expr size = "<<lhs.expSize() 
            <<" second expr size ="<<rhs.expSize() 
            <<iomanagment::endProgram;
        }
    }
    
    T evalAt(size_t pos) const { return crossProd(m_lhs.evalAt(pos),m_rhs.evalAt(pos));}
    size_t expSize() const { return m_lhs.expSize(); }
  
private:
    const ET_Array<VectorTX<T,2>,Lhs> &m_lhs;
    const ET_Array<VectorTX<T,2>,Rhs> &m_rhs;
};

template <typename T, typename Array1, typename Array2>
CrossProd2DExpression<T,Array1,Array2> crossProd(const ET_Array<VectorTX<T,2>, Array1>& list1, const ET_Array<VectorTX<T,2>, Array2>& list2)
{
    return CrossProd2DExpression<T,Array1,Array2>(list1,list2);
}

template<typename T,typename Derived>
struct CrossProd2DSingValueExpression : public ET_Array<T,CrossProd2DSingValueExpression<T,Derived> >
{
    CrossProd2DSingValueExpression(const ET_Array<VectorTX<T,2>,Derived> &list, const VectorTX<T,2> & singVal)
    : m_list(list), m_singVal(singVal)
    {
    }
    
    T evalAt(size_t pos) const 
    { 
        return crossProd(m_list.evalAt(pos), m_singVal);
    }
    
    size_t expSize() const { return m_list.expSize();}
    
private:
    VectorTX<T,2> m_singVal;
    const ET_Array<VectorTX<T,2>,Derived> &m_list;
};

template <typename T, typename Array>
CrossProd2DSingValueExpression<T,Array> crossProd(const ET_Array<VectorTX<T,2>, Array>& list, const VectorTX<T,2>& singVal)
{
    return CrossProd2DSingValueExpression<T,Array>(list,singVal);
}

//inverted order
template <typename T, typename Array>
CrossProd2DSingValueExpression<T,Array> crossProd(const VectorTX<T,2>& singVal,const ET_Array<VectorTX<T,2>, Array>& list)
{
    return CrossProd2DSingValueExpression<T,Array>(list,-singVal);
}


//------------- cross product 3d ------------------------------------------------------//
template<typename T, typename Lhs, typename Rhs>
struct CrossProd3DExpression : public ET_Array<VectorTX<T,3>,CrossProd3DExpression<T,Lhs,Rhs> >
{
    CrossProd3DExpression(const ET_Array<VectorTX<T,3>,Lhs> &lhs,const ET_Array<VectorTX<T,3>,Rhs> &rhs)
    : m_lhs(lhs), m_rhs(rhs)
    {
        if(lhs.expSize() != rhs.expSize())
        {
            ErrorInFunction<<" Inconsistient sizes between 2 expressions "
            <<" provided to calculate crossProd on them"
            <<" first expr size = "<<lhs.expSize() 
            <<" second expr size ="<<rhs.expSize() 
            <<iomanagment::endProgram;
        }
    }
    
    VectorTX<T,3> evalAt(size_t pos) const { return crossProd(m_lhs.evalAt(pos),m_rhs.evalAt(pos));}
    size_t expSize() const { return m_lhs.expSize(); }
    
private:
    const ET_Array<VectorTX<T,3>,Lhs> &m_lhs;
    const ET_Array<VectorTX<T,3>,Rhs> &m_rhs;
};

template <typename T, typename Array1, typename Array2>
CrossProd3DExpression<T,Array1,Array2> crossProd(const ET_Array<VectorTX<T,3>, Array1>& list1, const ET_Array<VectorTX<T,3>, Array2>& list2)
{
    return CrossProd3DExpression<T,Array1,Array2>(list1,list2);
}

template<typename T, typename Array>
struct CrossProd3DSingValueExpression : public ET_Array< VectorTX<T,3>, CrossProd3DSingValueExpression<T,Array> >
{
    CrossProd3DSingValueExpression(const ET_Array<VectorTX<T,3>,Array> &list, const VectorTX<T,3> & singVal)
    : m_list(list), m_singVal(singVal)
    {
    }
    
    VectorTX<T,3> evalAt(size_t pos) const 
    { 
        return crossProd(m_list.evalAt(pos), m_singVal);
    }
    
    size_t expSize() const { return m_list.expSize();}
    
private:
    VectorTX<T,3> m_singVal;
    const ET_Array<VectorTX<T,3>,Array> &m_list;
};

template<typename T, typename Array>
CrossProd3DSingValueExpression<T,Array> crossProd(const ET_Array<VectorTX<T,3>,Array> &list, const VectorTX<T,3> & singVal)
{
    return CrossProd3DSingValueExpression<T,Array>(list,singVal);
}

//inverted order
template<typename T,typename Array>
CrossProd3DSingValueExpression<T,Array> crossProd(const VectorTX<T,3> & singVal,const ET_Array<VectorTX<T,3>,Array> &list)
{
    return CrossProd3DSingValueExpression<T,Array>(list,-singVal);
}

//-------------------- magnituted of vector ------------------------//
template<typename T, size_t Dim, typename Array>
struct VectorMagnitude : public ET_Array<T,VectorMagnitude<T,Dim,Array> >
{
    VectorMagnitude(const ET_Array<VectorTX<T,Dim>,Array> & list): m_list(list){}
    size_t expSize() const { return m_list.expSize();}
    T evalAt(size_t pos) const { return m_list.evalAt(pos).mag();}
private:
    const ET_Array<VectorTX<T,Dim>,Array> & m_list;
};

template<typename T, size_t Dim, typename Array>
VectorMagnitude<T,Dim,Array> mag(const ET_Array<VectorTX<T,Dim>,Array> & list)
{
    return VectorMagnitude<T,Dim,Array>(list);
}

//-------------------- squared magnituted of vector ------------------------//
template<typename T, size_t Dim, typename Array>
struct VectorMagnitudeSqrt : public ET_Array<T,VectorMagnitudeSqrt<T,Dim,Array> >
{
    VectorMagnitudeSqrt(const ET_Array<VectorTX<T,Dim>,Array> & list): m_list(list){}
    size_t expSize() const { return m_list.expSize();}
    T evalAt(size_t pos) const 
    { 
        VectorTX<T,Dim> element = m_list.evalAt(pos);
        T magS = 0;
        for(unsigned int i=0; i<Dim; ++i)
            magS = element[i]*element[i];
        
        return  magS;
    }
private:
    const ET_Array<VectorTX<T,Dim>,Array> & m_list;
};

template<typename T, size_t Dim, typename Array>
VectorMagnitudeSqrt<T,Dim,Array> magSqrt(const ET_Array<VectorTX<T,Dim>,Array> & list)
{
    return VectorMagnitudeSqrt<T,Dim,Array>(list);
}



//-------------- absolut values of vectors in list ----------------------//
template<typename T, size_t Dim, typename Array>
struct ET_AbsFromArrayOfVectors : public ET_Array<VectorTX<T,Dim>,ET_AbsFromArrayOfVectors<T,Dim,Array> >
{
    typedef ET_Array<VectorTX<T,Dim>,Array> ArrayType;
    
    ET_AbsFromArrayOfVectors(const ArrayType &array):
    m_array(array)
    {
    }
    
    size_t expSize() const  { return m_array.expSize();}
    VectorTX<T,Dim> evalAt(size_t index) const 
    {
        VectorTX<T,Dim> value = m_array.evalAt(index);
        for(size_t i=0;i<Dim; i++)
            value[i]=std::abs(value[i]);
        
        return value;
    }
    
private:
    const ArrayType &m_array;
};

template<typename T, size_t Dim, typename Array>
ET_AbsFromArrayOfVectors<T,Dim,Array> abs(const ET_Array<VectorTX<T,Dim>,Array> & arrayExpr)
{
    return ET_AbsFromArrayOfVectors<T,Dim,Array>(arrayExpr);
}




}//array



} //SEM


#endif //VECTOR_UTILS_H