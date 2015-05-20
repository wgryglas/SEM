
#ifndef NUMARRAYEIGEN_H
#define NUMARRAYEIGEN_H

#include "Eigen/Dense"

#include "utilities/ET_Array.h"
#include "utilities/numArrayBase.h"


namespace SEM {

template<typename T>
class numArrayEigen;

template<typename T>
struct BaseTraits<numArrayEigen<T> >
{
    typedef T value_type;
    typedef T& reference;
    typedef const T& const_reference;
};


template<typename T>
class numArrayEigen : public numArrayBase<numArrayEigen<T> >
{
    
public:
    typedef Eigen::Matrix<double,Eigen::Dynamic,1> EigenData;
private:
    typedef numArrayEigen<T> type;
    typedef numArrayBase<numArrayEigen<T> > baseType;
    
    EigenData eData;
    
public:
    typedef typename BaseTraits<type>::value_type value_type;
    typedef typename BaseTraits<type>::reference reference;
    typedef typename BaseTraits<type>::const_reference const_reference;
    
    
    numArrayEigen(const size_t& size=0)
    : eData(size)
    {
    }
    
    numArrayEigen(const Eigen::Matrix<double,Eigen::Dynamic,1> & eValues)
    :eData(eValues)
    {
    }
    
    template<typename OtherExpr> 
    numArrayEigen(const array::ET_Array<T,OtherExpr> &other)
        :eData(other.expSize())
    {
        for(int i=0;i<size(); ++i)
            eData[i]=other.evalAt(i);
    }
    
    reference operator [](size_t index) { return eData[index];}
    const_reference operator [](size_t index) const { return eData[index];}
    size_t size() const {return eData.size(); }
    
    value_type evalAt(size_t index) const { return eData[index]; }
    size_t expSize() const { return eData.size();}

    
    void resize(size_t size) { eData.resize(size); }
    EigenData & eigenData() { return eData; }
    const EigenData & eigenData() const { return eData; }

    ~numArrayEigen(){}
    
    
    CONFIGURE_ASSIGMENT_OPERATORS_IN_DERIVED_FROM_NUMARRAYBASE( value_type, type, baseType )

};

}//SEM
#endif // NUMARRAYEIGEN_H
