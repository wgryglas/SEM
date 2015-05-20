#ifndef numArraySlice_H
#define numArraySlice_H

#include <vector>

#include "numArrayDerived.h"

namespace SEM
{
/// \file numArraySlices - file containing simple
///       implementations of numArrayBase and numArray2DBase.
///       Those implementations are aimed to be special kind
///       of reference to real data storage numArray. Such a
///       reference would allow to treat sliced numArray in
///       the same manner as it would be real data storege.
///       References are made by index mapping from real array
///       (or even mapping from mapping - making slice form another
///        slice) element indexes. Mapping is performd by defining
///       some mask with int values or by range -start-end pairs.
///
///       Generaly it can be very useful if somen would want to
///       operate only on part of great array, and refer in subArray
///       by it's local indexes - eg. using mes solution values
///       stored in global array by referencing in local element
///       indexing (even ref. to elemetal values like 2Darray if
///       you would make slice2D from 1D array).
///
///       Below classes are implementation of numArrayBase and
///       numArray2DBase. But below classes are also used by
///       numArrayBase so, in fact base class will return it's
///       derived class(when perform slice method), so all below
///       classes need to be forward declared to allow compiler
///       include this file without circular reference with numArrayBase.h


/// \brief forward declaration of numArrayBase
///        to be able to prepare classes which
///        would be used as slice of numArrayBase
///        (in fact, it would be inmplementation
///         of numArrayBase)
template<typename DerivedData>
class numArrayBase;


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
class BaseTraits;

//--------------------------------------------------------------------------------------------//
//                                  1D SLICES
//--------------------------------------------------------------------------------------------//

//------------------------------------------------------------------------------//
//          RANGE SPECIFICATION MAPPING ARRAY INTO ANOTHER SUB-ARRAY 
//------------------------------------------------------------------------------//
/// \class numArrayRangeMapped
/// \brief Implementation of numArrayBase - some kind of array wrapper
///        which refers to data elements and behaves in the same way
///        as ordinary array. Reference to data is made by range of
///        array, so this class can be called as "subArray from-to"
template<typename DataType>
class numArrayRangeMapped;

template<typename DataType>
struct BaseTraits<numArrayRangeMapped<DataType> >
{
    typedef typename BaseTraits<DataType >::value_type value_type;
    typedef typename BaseTraits<DataType >::reference reference;
    typedef typename BaseTraits<DataType >::const_reference const_reference;
};

template<typename DataType>
struct BaseTraits<numArrayRangeMapped<const DataType> >
{
    typedef typename BaseTraits<DataType >::value_type value_type;
    typedef typename BaseTraits<DataType >::const_reference reference;
    typedef typename BaseTraits<DataType >::const_reference const_reference;
};

template<typename DataType>
class numArrayRangeMapped : public numArrayBase<numArrayRangeMapped<DataType> >
{
   typedef numArrayRangeMapped<DataType> type;
   typedef numArrayBase<type > baseType;

   DataType & m_refData;
   size_t m_start;
   size_t m_size;

public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;

   numArrayRangeMapped(DataType &data, const size_t & start, const size_t & end)
       : m_refData(data),m_start(start), m_size(end-start+1)
   {
   }

   numArrayRangeMapped(const type &other)
       : m_refData(other.m_refData), m_start(other.m_start), m_size(other.m_size)
   {
   }

//   /// direct copy constructor is not allowed, because it's kind of reference
//   /// so it can't change it's data reference at assigment, instead it must
//   /// assigne values from other data to this data elenets. This kind of stuff
//   /// is performed in operators below generated by macro
//   /// CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE
//   void copyMapping(const type &other)
//   {
//        m_start = other.m_start;
//        m_start = other.m_size;
//   }

   //-----------------------------------------------------------------------//
   //                   numArrayBase IMPLEMENTATION                         //
   //-----------------------------------------------------------------------//
   inline const_reference operator[](const size_t &loc) const
   {
       if(loc >= m_size)
           ErrorInFunction<<"Index out of range in RangeMapper\n"
                        <<"asked for "<<loc<<" element and size is"<<m_size<<"\n"
                        <<"-class mapping 1D data to another 1D data"
                        <<iomanagment::endProgram;

       return m_refData[m_start+loc];
   }

   inline reference operator[](const size_t &loc)
   {
       if(loc >= m_size)
           ErrorInFunction<<"Index out of range in RangeMapper\n"
                        <<"asked for "<<loc<<" element and size is "<<m_size<<"\n"
                        <<"-class mapping 1D data to another 1D data"
                        <<iomanagment::endProgram;

       return m_refData[m_start+loc];
   }

   inline size_t size() const {return m_size;}
   //-----------------------------------------------------------------------//
   //                   configure operators                                 //
   //-----------------------------------------------------------------------//
   CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, type, baseType )
};

//------------------------------------------------------------------------------//
//          INDEX MAPPING ARRAY INTO ANOTHER SUB-ARRAY
//------------------------------------------------------------------------------//

/// \class numArrayIndexMapped
/// \brief Implementation of numArrayBase - some kind of array wrapper
///        which refers to data elements and behaves in the same way
///        as ordinary array. Referance to data made by index mapping,
///        so this class is subArray, but access to original data
///        elements is mapped by indexes defined as input to this class.
template<typename DataType>
class numArrayIndexMapped;

template<typename DataType>
struct BaseTraits<numArrayIndexMapped<DataType> >
{
    typedef typename BaseTraits<DataType >::value_type value_type;
    typedef typename BaseTraits<DataType >::reference reference;
    typedef typename BaseTraits<DataType >::const_reference const_reference;
};

template<typename DataType>
struct BaseTraits<numArrayIndexMapped<const DataType> >
{
    typedef typename BaseTraits<DataType >::value_type value_type;
    typedef typename BaseTraits<DataType >::const_reference reference;
    typedef typename BaseTraits<DataType >::const_reference const_reference;
};

template<typename DataType>
class numArrayIndexMapped : public numArrayBase<numArrayIndexMapped<DataType> >
{
    typedef numArrayIndexMapped<DataType> type;
    typedef numArrayBase<type> baseType;

    DataType & m_refData;
    const VectorToVectorMap m_map;

public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;

    numArrayIndexMapped(DataType &data, const VectorToVectorMap & mapping)
        : m_refData(data),m_map(mapping)
    {
    }

    numArrayIndexMapped(const type &other)
        : m_refData(other.m_refData), m_map(other.m_map)
    {
    }

//    /// direct copy constructor is not allowed, because it's kind of reference
//    /// so it can't change it's data reference at assigment, instead it must
//    /// assigne values from other data to this data elenets. This kind of stuff
//    /// is performed in operators below generated by macro
//    /// CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE
//    void copyMapping(const type &other)
//    {
//         m_map = other.m_map;
//    }

    //-----------------------------------------------------------------------//
    //                   numArrayBase IMPLEMENTATION                         //
    //-----------------------------------------------------------------------//
    inline const_reference operator[](const size_t &loc) const
    {
        return m_refData[m_map[loc]];
    }

    inline reference operator[](const size_t &loc)
    {
        return m_refData[m_map[loc]];
    }

    inline size_t size() const { return m_map.size(); }
    //-----------------------------------------------------------------------//
    //                   configure operators                                 //
    //-----------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, type, baseType )
};




//------------------------------------------------------------------------------//
//          INDEX MAPPING 2D-ARRAY INTO ARRAY
//------------------------------------------------------------------------------//

template<typename Data2D>
class numArrayData2DMapped;

template<typename Data2D>
struct BaseTraits<numArrayData2DMapped<Data2D> >
{
    typedef typename BaseTraits<typename Data2D::value_type >::value_type value_type;
    typedef typename BaseTraits<typename Data2D::value_type >::reference reference;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference const_reference;
};

template<typename Data2D>
struct BaseTraits<numArrayData2DMapped<const Data2D> >
{
    typedef typename BaseTraits<typename Data2D::value_type >::value_type value_type;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference reference;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference const_reference;
};


/// \brief numArrayData2DMapped
/// class which maps some numArray2D data into array 1D data.
template<typename Data2D >
class numArrayData2DMapped : public numArrayBase<numArrayData2DMapped<Data2D> >
{
    typedef numArrayData2DMapped<Data2D> type;
    typedef numArrayBase<numArrayData2DMapped<Data2D> > baseType;

    Data2D & m_refData;
    VectorToMatrixMap m_map;
    
    //disallow asigment from other in this way (only by expression)
    numArrayData2DMapped<Data2D> operator = (const numArrayData2DMapped<Data2D>&);
    
public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;

    numArrayData2DMapped(Data2D & data, const VectorToMatrixMap & mapping)
        : m_refData(data), m_map(mapping)
    {
    }
    
    numArrayData2DMapped(const numArrayData2DMapped<Data2D> & other)
    : m_refData(other.m_refData), m_map(other.m_map)
    {
    }
    
    //-----------------------------------------------------------------------//
    //                   numArrayBase IMPLEMENTATION                         //
    //-----------------------------------------------------------------------//
    inline const_reference operator[](const size_t &loc) const
    {
        return m_refData[ m_map[loc][0] ][ m_map[loc][1] ];
    }

    inline reference operator[](const size_t &loc)
    {
        return m_refData[ m_map[loc][0] ][ m_map[loc][1] ];
    }

    inline size_t size() const { return m_map.size(); }
    //-----------------------------------------------------------------------//
    //                   configure operators                                 //
    //-----------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, type, baseType )
};


//------------------------------------------------------------------------------//
//          DIAGONAL ARRAY FROM ARRAY-2D
//------------------------------------------------------------------------------//
template<typename Data2D>
class numArrayDiagonalFrom2D;

template<typename Data2D>
struct BaseTraits<numArrayDiagonalFrom2D<Data2D> >
{
    typedef typename BaseTraits<typename Data2D::value_type >::value_type value_type;
    typedef typename BaseTraits<typename Data2D::value_type >::reference reference;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference const_reference;
};

template<typename Data2D>
struct BaseTraits<numArrayDiagonalFrom2D<const Data2D> >
{
    typedef typename BaseTraits<typename Data2D::value_type >::value_type value_type;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference reference;
    typedef typename BaseTraits<typename Data2D::value_type >::const_reference const_reference;
};

template<typename Data2D>
class numArrayDiagonalFrom2D : public numArrayBase<numArrayDiagonalFrom2D<Data2D> >
{
    typedef numArrayDiagonalFrom2D<Data2D> type;
    typedef numArrayBase<type> baseType;
    Data2D & m_data;
public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;
    
    numArrayDiagonalFrom2D(Data2D & data)
    :m_data(data)
    {
        using namespace iomanagment;
        size_t sMin = data.size();
        for(size_t i=0; i<sMin; ++i)
            if(data[i].size() < sMin)
            {
                ErrorInFunction<<"Can't make diagonal subArray from 2D data, where "<<i
                               <<"column size is smaller then number of columns"
                               <<endProgram;
            }
    }
    
    ~numArrayDiagonalFrom2D()
    {
    }
    //---------------------------------------------------------------------------//
    //              numArrayBase implemementation
    //---------------------------------------------------------------------------//
    reference operator [](size_t index) { return m_data[index][index];}
    const_reference operator [](size_t index) const { return m_data[index][index];}
    size_t size() const {return m_data.size(); }
    //---------------------------------------------------------------------------//
    //               numArryBase operators config.
    //---------------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(value_type,type,baseType)
};




//--------------------------------------------------------------------------------------------//
//                                  2D SLICES
//--------------------------------------------------------------------------------------------//

/// \brief forward declaration of numArray2DBase
template<typename DerivedData>
class numArray2DBase;

/// \brief forward declaration of BaseTraits2D,
///        simmilar description as BaseTraits
template<typename Derived>
class BaseTraits2D;

/// \class numArray2DIndexMapped
/// \brief Implementation of numArray2DBase - some kind of array2D wrapper
///        which refers to data elements and behaves in the same way
///        as ordinary array2D. Reference to data is made by ranges of
///        array2D, so this class can be called as "subArray2D from-to"
///        Sub element, get by operator [] is in fact 1D arrayIndexMapped
///        type.
template<typename Data>
class numArray2DIndexMapped;

template<typename Data>
struct BaseTraits2D<numArray2DIndexMapped<Data> >
{
    typedef numArrayData2DMapped<Data> value_type;
    typedef numArrayData2DMapped<Data>& reference;
    typedef const numArrayData2DMapped<Data>& const_reference;

    typedef typename BaseTraits2D<Data>::inner_type inner_type;
    typedef typename BaseTraits2D<Data>::inner_reference inner_reference;
    typedef typename BaseTraits2D<Data>::const_inner_reference const_inner_reference;
};

template<typename Data>
struct BaseTraits2D<numArray2DIndexMapped<const Data> >
{
    typedef numArrayData2DMapped<Data> value_type;
    typedef numArrayData2DMapped<Data>& reference;
    typedef const numArrayData2DMapped<Data>& const_reference;

    typedef typename BaseTraits2D<Data>::inner_type inner_type;
    typedef typename BaseTraits2D<Data>::const_inner_reference inner_reference;
    typedef typename BaseTraits2D<Data>::const_inner_reference const_inner_reference;
};


template<typename Data>
class numArray2DIndexMapped : public numArray2DBase<numArray2DIndexMapped<Data> >
{
    typedef numArray2DIndexMapped<Data> type;
    typedef numArray2DBase<type> baseType;

    std::vector< typename BaseTraits2D<type>::value_type* > m_mapping;

public:
    typedef typename BaseTraits2D<type>::inner_type inner_type;
    typedef typename BaseTraits2D<type>::inner_reference inner_reference;
    typedef typename BaseTraits2D<type>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<type>::value_type value_type;
    typedef typename BaseTraits2D<type>::reference reference;
    typedef typename BaseTraits2D<type>::const_reference const_reference;

    //------------------------------------------------------//
    //              numArray2DBase implementation           //
    //------------------------------------------------------//
    reference operator[](const size_t &index) { return *m_mapping[index]; }
    const_reference operator[](const size_t &index) const { return *m_mapping[index]; }
    size_t size() const {return m_mapping.size(); }

    numArray2DIndexMapped(Data& data, const MatrixToMatrixMap& map)
        : m_mapping(map.size())
    {
        for(size_t i=0; i<map.size(); ++i)
        {
            m_mapping[i] = new value_type(data,map[i]);
        }
    }

    ~numArray2DIndexMapped()
    {
        for(typename std::vector<value_type*>::iterator it=m_mapping.begin(); it!=m_mapping.end(); ++it)
            delete *it;

        m_mapping.clear();
    }

    //--------------------------------------------------------//
    //       configure derived class from numArray2DBase      //
    //--------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(inner_type,type,baseType)
};

/// \class numArrayIndexMapped
/// \brief Implementation of numArrayBase - some kind of array wrapper
///        which refers to data elements and behaves in the same way
///        as ordinary array. Referance to data made by index mapping,
///        so this class is subArray, but access to original data
///        elements is mapped by indexes defined as input to this class.
///        Sub element, get by operator [] is in fact 1D arrayRangeMapped
///        type.
template<typename Data>
class numArray2DRangeMapped;

template<typename Data>
struct BaseTraits2D<numArray2DRangeMapped<Data> >
{
    typedef typename BaseTraits2D<Data>::inner_type inner_type;
    typedef typename BaseTraits2D<Data>::inner_reference inner_reference;
    typedef typename BaseTraits2D<Data>::const_inner_reference const_inner_reference;

    typedef numArrayRangeMapped<typename BaseTraits2D<Data>::value_type> value_type;
    typedef numArrayRangeMapped<typename BaseTraits2D<Data>::value_type>& reference;
    typedef const numArrayRangeMapped<typename BaseTraits2D<Data>::value_type>& const_reference;
};

template<typename Data>
struct BaseTraits2D<numArray2DRangeMapped<const Data> >
{
    typedef typename BaseTraits2D<Data>::inner_type inner_type;
    typedef typename BaseTraits2D<Data>::const_inner_reference inner_reference;
    typedef typename BaseTraits2D<Data>::const_inner_reference const_inner_reference;

    typedef numArrayRangeMapped<typename BaseTraits2D<Data>::value_type> value_type;
    typedef numArrayRangeMapped<typename BaseTraits2D<Data>::value_type>& reference;
    typedef const numArrayRangeMapped<typename BaseTraits2D<Data>::value_type>& const_reference;
};


template<typename Data>
class numArray2DRangeMapped : public numArray2DBase<numArray2DRangeMapped<Data> >
{
    typedef numArray2DRangeMapped<Data> type;
    typedef numArray2DBase<type > baseType;

    std::vector<typename BaseTraits2D<type>::value_type*> m_mapping;

public:
    typedef typename BaseTraits2D<type>::inner_type inner_type;
    typedef typename BaseTraits2D<type>::inner_reference inner_reference;
    typedef typename BaseTraits2D<type>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<type>::value_type value_type;
    typedef typename BaseTraits2D<type>::reference reference;
    typedef typename BaseTraits2D<type>::const_reference const_reference;
    //------------------------------------------------------//
    //              numArray2DBase implementation           //
    //------------------------------------------------------//
    reference operator[](const size_t &index) { return *m_mapping[index]; }
    const_reference operator[](const size_t &index) const { return *m_mapping[index]; }
    size_t size() const {return m_mapping.size(); }

    numArray2DRangeMapped( Data& data,
            const size_t &start1D, const size_t &end1D,
            const size_t &start2D, const size_t &end2D
            ): m_mapping(end1D-start1D+1)
    {
        if(data.size()<end1D)
        {
            ErrorInFunction<<"mapping inconsitient with provided data in"<<iomanagment::endProgram;
        }


        for(size_t i=0; i<m_mapping.size(); ++i)
        {
            if(data[i].size()<end2D)
            {
                ErrorInFunction<<"mapping inconsitient with provided data in"
                               <<iomanagment::endProgram;
            }
            m_mapping[i] = new value_type(data[start1D+i],start2D,end2D) ;
        }
    }

    ~numArray2DRangeMapped()
    {
        for(typename std::vector<value_type*>::iterator it=m_mapping.begin(); it!=m_mapping.end(); ++it)
            delete *it;

        m_mapping.clear();
    }

    //--------------------------------------------------------//
    //       configure derived class from numArray2DBase      //
    //--------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(inner_type,type,baseType)
};





template<typename Data1D>
class numArray2DData1DIndexMapped;

template<typename Data1D>
struct BaseTraits2D<numArray2DData1DIndexMapped<Data1D> >
{
    typedef typename BaseTraits<Data1D>::value_type inner_type;
    typedef typename BaseTraits<Data1D>::reference inner_reference;
    typedef typename BaseTraits<Data1D>::const_reference const_inner_reference;

    typedef numArrayIndexMapped<Data1D> value_type;
    typedef numArrayIndexMapped<Data1D>& reference;
    typedef const numArrayIndexMapped<Data1D>& const_reference;
};

template<typename Data1D>
struct BaseTraits2D<numArray2DData1DIndexMapped<const Data1D> >
{
    typedef typename BaseTraits<Data1D>::value_type inner_type;
    typedef typename BaseTraits<Data1D>::const_reference inner_reference;
    typedef typename BaseTraits<Data1D>::const_reference const_inner_reference;

    typedef numArrayIndexMapped<Data1D> value_type;
    typedef numArrayIndexMapped<Data1D>& reference;
    typedef const numArrayIndexMapped<Data1D>& const_reference;
};

template<typename Data1D>
class numArray2DData1DIndexMapped : public numArray2DBase<numArray2DData1DIndexMapped<Data1D> >
{
    typedef numArray2DData1DIndexMapped<Data1D> type;
    typedef numArray2DBase<type> baseType;

    std::vector<typename BaseTraits2D<type>::value_type*> m_mapping;

public:
    typedef typename BaseTraits2D<type>::inner_type inner_type;
    typedef typename BaseTraits2D<type>::inner_reference inner_reference;
    typedef typename BaseTraits2D<type>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<type>::value_type value_type;
    typedef typename BaseTraits2D<type>::reference reference;
    typedef typename BaseTraits2D<type>::const_reference const_reference;

    numArray2DData1DIndexMapped(Data1D & data, const MatrixToVectorMap & mapping)
        :m_mapping(mapping.size())
    {
        for(size_t i=0; i<mapping.size(); ++i)
            m_mapping[i] = new value_type(data,mapping[i]);
    }

    ~numArray2DData1DIndexMapped()
    {
        for(typename std::vector<value_type*>::iterator it=m_mapping.begin(); it!=m_mapping.end(); ++it)
            delete (*it);
        m_mapping.clear();
    }

    //---------------------------------------------------------------------------//
    //              numArrayBase implemementation
    //---------------------------------------------------------------------------//
    reference operator [](size_t index) { return *m_mapping[index];}
    const_reference operator [](size_t index) const { return *m_mapping[index];}
    size_t size() const {return m_mapping.size(); }
    //---------------------------------------------------------------------------//
    //               numArryBase operators config.
    //---------------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(inner_type,type,baseType)
};




template<typename Data1D>
class numArray2DData1DSplited;

template<typename Data1D>
struct BaseTraits2D<numArray2DData1DSplited<Data1D> >
{
    typedef typename BaseTraits<Data1D>::value_type inner_type;
    typedef typename BaseTraits<Data1D>::reference inner_reference;
    typedef typename BaseTraits<Data1D>::const_reference const_inner_reference;

    typedef numArrayRangeMapped<Data1D> value_type;
    typedef numArrayRangeMapped<Data1D>& reference;
    typedef const numArrayRangeMapped<Data1D>& const_reference;
};

template<typename Data1D>
struct BaseTraits2D<numArray2DData1DSplited<const Data1D> >
{
    typedef typename BaseTraits<Data1D>::value_type inner_type;
    typedef typename BaseTraits<Data1D>::const_reference inner_reference;
    typedef typename BaseTraits<Data1D>::const_reference const_inner_reference;

    typedef numArrayRangeMapped<Data1D> value_type;
    typedef numArrayRangeMapped<Data1D>& reference;
    typedef const numArrayRangeMapped<Data1D>& const_reference;
};


template<typename Data1D>
class numArray2DData1DSplited : public numArray2DBase<numArray2DData1DSplited<Data1D> >
{
    typedef numArray2DData1DSplited<Data1D> type;
    typedef numArray2DBase<type> baseType;

    std::vector<typename BaseTraits2D<type>::value_type*> m_mapping;
public:
    typedef typename BaseTraits2D<type>::inner_type inner_type;
    typedef typename BaseTraits2D<type>::inner_reference inner_reference;
    typedef typename BaseTraits2D<type>::const_inner_reference const_inner_reference;
    typedef typename BaseTraits2D<type>::value_type value_type;
    typedef typename BaseTraits2D<type>::reference reference;
    typedef typename BaseTraits2D<type>::const_reference const_reference;

    numArray2DData1DSplited(Data1D & data, const std::vector<size_t> & sizes)
        :m_mapping(sizes.size())
    {
        size_t start = 0;
        for(size_t i=0; i<sizes.size(); ++i)
        {
            m_mapping[i] = new value_type(data,start,start+sizes[i]-1);
            start+=sizes[i];
        }
    }


    ~numArray2DData1DSplited()
    {
        for(typename std::vector<value_type*>::iterator it=m_mapping.begin(); it!=m_mapping.end(); ++it)
            delete (*it);
        m_mapping.clear();
    }
    //---------------------------------------------------------------------------//
    //              numArrayBase implemementation
    //---------------------------------------------------------------------------//
    reference operator [](size_t index) { return *m_mapping[index];}
    const_reference operator [](size_t index) const { return *m_mapping[index];}
    size_t size() const {return m_mapping.size(); }
    //---------------------------------------------------------------------------//
    //               numArryBase operators config.
    //---------------------------------------------------------------------------//
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE(inner_type,type,baseType)
};




}//SEM

#endif // numArraySlice_H
