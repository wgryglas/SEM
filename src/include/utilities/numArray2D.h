#ifndef NUMARRAY2D_H
#define NUMARRAY2D_H

#include <vector>
#include <iostream>
#include <string>

#include "numArray.h"
#include "numArray2DBase.h"
#include "iomanagment/InfoStream.h"
#include "iomanagment/ReadWriteStream.h"

namespace SEM {

template<typename T>
class numArray2D;

template<typename T>
struct BaseTraits2D<numArray2D<T> >
{
    typedef T inner_type;
    typedef typename BaseTraits<numArray<T> >::reference inner_reference;
    typedef typename BaseTraits<numArray<T> >::const_reference const_inner_reference;

    typedef numArray<T> value_type;
    typedef numArray<T>& reference;
    typedef const numArray<T>& const_reference;
};


template<typename T>
class numArray2D : public std::vector<typename BaseTraits2D<numArray2D<T> >::value_type>, public numArray2DBase<numArray2D<T> >
{
    REFERENCE_TYPE(numArray2D<T>)
private:
    typedef std::vector<typename BaseTraits2D<numArray2D<T> >::value_type> baseDataType;
    typedef numArray2DBase<numArray2D<T> > array2DBaseType;
    typedef numArray2D<T> type;
public:
    typedef typename BaseTraits2D<type>::inner_type inner_type;
    typedef typename BaseTraits2D<type>::inner_reference inner_reference;
    typedef typename BaseTraits2D<type>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<type>::value_type value_type;
    typedef typename BaseTraits2D<type>::reference reference;
    typedef typename BaseTraits2D<type>::const_reference const_reference;


    numArray2D(const size_t &size1=0): baseDataType(size1)
    {
    }

    numArray2D(const size_t &size1, const size_t &size2): baseDataType(size1,value_type(size2))
    {
    }
    
    numArray2D(const size_t &size1, const size_t &size2, const T value): baseDataType(size1, value_type(size2,value))
    {
    }

    /// \brief numArray2D copy constructor from other
    /// \param other
    numArray2D(const numArray2D<T> & other): baseDataType(other)
    {
    }

    /// \brief numArray2D - copy constructor
    /// \param other - other array2DBase - can be mapped type
    template<typename DerivedType>
    numArray2D(const numArray2DBase<DerivedType> & other):baseDataType(other.size())
    {
        for(int i=0; i<size();++i)
        {
            (*this)[i].resize(other[i].size());
            (*this)[i] = other[i];
        }
    }

    /// \brief operator = assigment operator form other numArray2D
    /// \param other
    /// \return
    numArray2D<T> & operator =(const numArray2D<T>& other)
    {
        checkOtherarray2DSize(other);
        baseDataType::operator =(other);
        return *this;
    }

    /// \brief resize2D -method to resize 2 dimmension as squared array
    /// \param size1  -1st dimmension size
    /// \param size2  -2st dimmension size
    void resize2D(size_t size1, size_t size2)
    {
        baseDataType::resize(size1);

        for(int i=0; i<size1; ++i)
            (*this)[i].resize(size2);
    }

    /// \brief resize2D -method to resize 2 dimmension as squared array and set value 
    /// \param size1  -1st dimmension size
    /// \param size2  -2st dimmension size
    void resize2D(size_t size1, size_t size2, const T& val)
    {
        baseDataType::resize(size1);
        
        for(int i=0; i<size1; ++i)
            (*this)[i].resize(size2,val);
    }
    

    //--------------------------------------------------//
    //          numArray2DBase implementation           //
    //--------------------------------------------------//
    using baseDataType::operator [];
    using baseDataType::size;
//    const_reference operator[](size_t index) const { return baseDataType::operator[](index);}
//    reference operator[](size_t index) { return baseDataType::operator[](index);}
//    size_t size() const { return baseDataType::size();}

    //--------------------------------------------------//
    //          std::vector utils                       //
    //--------------------------------------------------//
    using baseDataType::resize;
    using baseDataType::clear;
    
    //--------------------------------------------------//
    //              EXPRESSION HANDLING                 //
    //--------------------------------------------------//
    /// \brief copy constructor from expression not allowed
    ///        because 2nd dimension typed won't be resolved.

    /// add operators for assigment from any expression templates array
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(inner_type,type, array2DBaseType)

private:
    template <typename OtherDerivedData>
    void checkOtherarray2DSize(const numArray2DBase<OtherDerivedData> & other) const
    {
        using namespace iomanagment;

        std::string msg = "Other numArray2DBase object has \n"
                          "incoherent size with this class object, \n"
                          "can't do assigment from other";

        if(size()==other.size())
        {
            for(int i=0; i<size();++i)
            {
                if((*this)[i].size()!=other[i].size())
                {
                    ErrorInFunction<<msg<<endProgram;
                }
            }
        }
        else
            ErrorInFunction<<msg<<endProgram;
    }
};

namespace iomanagment {

/// \brief Specialization writHelper structure to allow using write with numArray2D<T>
///        it can't be used now because of this class definition - it's not nasted, so
///        compiler can't figure out that TEMPLATE numArray2D<T> is also std::vector<numArray<T> >
///        which is nasted(nasted strucutres read/write functions can handle). Template resolver
///        don't take into acount type hierachy
template<typename T>
struct writeHelper<numArray2D<T> >
{
    static void writeValue(const numArray2D<T>& val, std::ostream& s)
    {
        write<std::vector<typename numArray2D<T>::value_type> >(val,s);
    }
};

/// \brief Specialization readHelper structure to allow using read with numArray2D<T>
///        it can't be used now because of this class definition - it's not nasted, so
///        compiler can't figure out that TEMPLATE numArray2D<T> is also std::vector<numArray<T> >
///        which is nasted(nasted strucutres read/write functions can handle). Template resolver
///        don't take into acount type hierachy
template<typename T>
struct readHelper<SEM::numArray2D<T> >
{
    static void readValue(std::istream& s,numArray2D<T>& val)
    {
        read<std::vector<typename numArray2D<T>::value_type> >(s,val);
    }
};



}//iomanagment
}//SEM



#endif // NUMARRAY2D_H
