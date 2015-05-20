#ifndef DISCONTINOUSFIELD_H
#define DISCONTINOUSFIELD_H

#include "ElementFieldBase.h"
#include "utilities/ET_Array.h"
#include "GeometricField.h"

namespace SEM { namespace field {

template<typename T>
class DiscontinousField;
    
template<typename T>
struct ElementFieldTraits<T,DiscontinousField<T> >
{
    typedef const numArray<T>& return_type;
};
    
template<typename T> 
class DiscontinousField :
    public ElementFieldBase<T,DiscontinousField<T> >, 
    public std::vector<numArray<T> >
{
    typedef ElementFieldBase<T,DiscontinousField<T> > ElementFieldType;
    typedef std::vector<numArray<T> > ArrayType;
    
public:
    using ElementFieldType::mesh;
    using ElementFieldType::elementsNumber;
    
    DiscontinousField(const mesh::Mesh & mesh)
    : ElementFieldType(mesh), ArrayType(mesh.size())
    {
        for(size_t e=0; e<mesh.size(); ++e)
            (*this)[e].resize(mesh[e].indexVectorMask().size());
    }
    
    DiscontinousField(const DiscontinousField& other)
    : ElementFieldType(other.mesh()),ArrayType(other)
    {
    }
    
    DiscontinousField<T> & operator=(const DiscontinousField& other)
    {
        ArrayType::operator=(other);
        return *this;
    }
    
    template<typename OtherDerived>
    DiscontinousField(const ElementFieldBase<T,OtherDerived> & other)
    : ElementFieldType(other.mesh()), ArrayType(other.elementsNumber())
    {
        for(size_t e=0; e<elementsNumber(); ++e)
        {
            (*this)[e].resize(mesh()[e].indexVectorMask().size());
            (*this)[e] = other.element(e);
        }
    }
    
    DiscontinousField(const GeometricField<T> & gField)
    : ElementFieldType(gField.mesh()), ArrayType(mesh().size())
    {
        for(size_t e=0; e<elementsNumber(); ++e)
        {
            (*this)[e].resize(mesh()[e].indexVectorMask().size());
            (*this)[e] = gField.element(e);
        }
    }
    
    template<typename DerivedArray>
    DiscontinousField(const array::ET_Array<T,DerivedArray> & expr,const mesh::Mesh& m)
    : ElementFieldType(m),ArrayType(m.size())
    {
        for(size_t e=0; e<mesh().size(); ++e)
        {
            const std::vector<int> & mask=mesh()[e].indexVectorMask();
            (*this)[e].resize(mask.size());
            
            for(size_t n=0; n<mask.size(); ++n)
            {
                (*this)[e][n] = expr.evalAt(mask[n]);
            }
        }
    }
    
    
    
    virtual ~DiscontinousField(){}
    
    const numArray<T> & element(size_t e) const 
    {
        return this->at(e); 
    }
    
    numArray<T> & element(size_t e)
    { 
        return this->at(e); 
    }
    
    //--------------------------------------------------------------
    //  ElementFieldBase assigment
    //-------------------------------------------------------------
    
#define SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(OP)                                       \
    template<typename OtherDerived>                                                         \
    DiscontinousField<T> & operator OP (const ElementFieldBase<T,OtherDerived> & other)     \
    {                                                                                       \
        for(size_t e=0; e<elementsNumber(); ++e)                                            \
        {                                                                                   \
            (*this)[e] OP other.element(e);                                                 \
        }                                                                                   \
        return *this;                                                                       \
    }                                                                                       

    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(+=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(-=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(*=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUS_FIELD(/=)
    
    #define SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(OP)                        \
    DiscontinousField<T> & operator OP (const T & val)                                      \
    {                                                                                       \
        for(size_t e=0; e<elementsNumber(); ++e)                                            \
        {                                                                                   \
            (*this)[e] OP val;                                                              \
        }                                                                                   \
        return *this;                                                                       \
    }                                                                                       
    
    SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(=)
    SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(+=)
    SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(-=)
    SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(*=)
    SEM_ALLOW_SINGLE_VAL_ASSIGMENT_TO_DISCONTINOUS_FIELD(/=)
    
    
    
    
#define SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(OP)                  \
    template<typename ET>                                                           \
    DiscontinousField<T> & operator OP (const array::ET_Array<T,ET> & continous)    \
    {                                                                               \
        for(size_t e=0; e<elementsNumber(); ++e)                                    \
        {                                                                           \
            const std::vector<int> &map =mesh()[e].indexVectorMask();               \
            auto& elementVals = element(e);                                         \
            for(size_t n=0; n<elementVals.size(); ++n )                             \
            {                                                                       \
                elementVals[n] OP continous.evalAt(map[n]);                         \
            }                                                                       \
        }                                                                           \
    }                                                                               
  
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(+=)  
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(-=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(*=)
    SEM_ALLOW_ASSIGMENT_TO_DISCONTINOUSFIELD_FROM_ET_ARRAY(/=)
};

}//field
}//SEM




#endif // DISCONTINOUSFIELD_H
