#include "components/Basic.h"
#include "utilities/numArray.h"
#include "utilities/Utilities.h"
#include "utilities/VectorUtils.h"
#include <cmath>

#ifndef CMPTRAITS_H
#define CMPTRAITS_H

namespace SEM {

/** **********************************************************
 * \struct CmpTraits - component traits provider.
 * Forward decalaration of basic components traits.
 * Type informations are provided by T explicit specialization
 * *********************************************************** */
template<typename T> 
class CmpTraits;

/** **********************************************************
 * \struct CmpTraits<Scalar> - specialization for Scalar type
 * *********************************************************** */
template<>
struct CmpTraits<Scalar>
{
    typedef const Scalar& const_Cmp;
    typedef Scalar& Cmp;
    typedef numArray<Scalar> CmpArray;
    typedef numArray<Scalar>& CmpArrayRef;
    typedef const numArray<Scalar>& const_CmpArrayRef;
    static const unsigned int DIM_SIZE =1;
    
    static Scalar zero() { return 0.;}
    
    static size_t dim()
    {
        return 1;
    }
    
    static const_Cmp component(const Scalar &val, const size_t& dim) 
    { 
        SEM_UNUSED(dim)
        return val;
    }
    
    static Cmp component(Scalar &val, const size_t& dim)
    { 
        SEM_UNUSED(dim)
        return val;
    }
    
    
    static const_CmpArrayRef cmpArray(const numArray<Scalar>& array, const size_t& dim)
    {
        SEM_UNUSED(dim)
        return array;
    }
    
    static CmpArrayRef cmpArray(numArray<Scalar>& array, const size_t& dim)
    {
        SEM_UNUSED(dim)
        return array;
    }
    
    static std::string fieldTypeName() 
    {
        return "scalar";
    }
    
    static bool isValid(const Scalar & value)
    {
        return ! std::isnan(value);
    }
    
};

/** **********************************************************
 * \struct CmpTraits<Vector> - specialization for Vector type
 * *********************************************************** */
template<size_t DIM>
struct CmpTraits<VectorTX<Scalar,DIM> >
{
    typedef const Scalar& Cmp;
    typedef numArrayVectorComponent<numArray<VectorTX<Scalar,DIM> > > CmpArray;
    typedef numArrayVectorComponent<numArray<VectorTX<Scalar,DIM> > > CmpArrayRef;
    typedef numArrayVectorComponent<const numArray<VectorTX<Scalar,DIM> > > const_CmpArrayRef;
   
    static const unsigned int DIM_SIZE = DIM;
    
    static VectorTX<Scalar,DIM> zero() { return VectorTX<Scalar,DIM>();}
    
    static int dim()
    {
        return DIM;
    }
    
    static const Scalar & component(const VectorTX<Scalar,DIM>& val, const size_t& dim)
    {
        return val[dim];
    }
    
    static Scalar & component(VectorTX<Scalar,DIM>& val, const size_t& dim)
    {
        return val[dim];
    }
    
    
    static const_CmpArrayRef cmpArray(const numArray<VectorTX<Scalar,DIM> >& list, const size_t& dim)
    {
        return const_CmpArrayRef(list,dim);
    }
    
    static CmpArrayRef cmpArray(numArray<VectorTX<Scalar,DIM> >& list, const size_t &dim)
    {
        return CmpArrayRef(list,dim);
    }
    
    static std::string fieldTypeName() 
    {
        return "vector";
    }
    
    static bool isValid(const VectorTX<Scalar,DIM> & value)
    {
        for(unsigned int d=0; d<DIM; ++d)
            if( std::isnan(value[d]) )
                return false;
        
        return true;
    }
    
};

}//SEM


#endif // CMPTRAITS_H
