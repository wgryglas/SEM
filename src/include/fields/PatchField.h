#ifndef PATCHFIELD_H
#define PATCHFIELD_H

#include <map>
#include <string>

#include "iomanagment/RegistryFile.h"
#include "iomanagment/Register.h"

#include "components/Basic.h"
#include "utilities/Reference.h"

#include "mesh/Boundary.h"

#include "ContinousField.h"

#include "utilities/numArray2D.h"

#include "time/TimeChangeListner.h"

namespace SEM { namespace field {

extern const std::string DIRICHLET;
extern const std::string NEUMANN;
    
template<typename T>
class PatchField :  public numArray2D<T>, public time::TimeChangeListner//public ContinousField<T>,
{
    REFERENCE_TYPE(PatchField<T>)

    //typedef ContinousField<T> baseType;
    typedef numArray2D<T> baseType;
    
    const mesh::Boundary &m_edges;

public:
    //using baseType::setRegistryFile;

//     PatchField(iomanagment::RegistryFile::ref regFile, const mesh::Boundary &edges)
//         : baseType(edges.size()), m_edges(edges)
//     {
//         setRegistryFile(regFile);
//     }

    PatchField(const mesh::Boundary &edges)
        : baseType(edges.size()), m_edges(edges)
    {
        for(size_t e=0; e<edges.size(); ++e)
            (*this)[e].resize(edges[e].gNodesIds().size());
    }
    
    PatchField(const mesh::Boundary &edges, size_t s1, size_t s2)
    : baseType(s1,s2), m_edges(edges)
    {
    }


    PatchField(const PatchField& other)
        : baseType(other), m_edges(other.m_edges)
    {
    }

    virtual ~PatchField(){}

    const mesh::Boundary & boundaryEdges() const { return m_edges;}

    std::string name() const { return m_edges.name(); }
    
    // ----------------------------------------------------------------//
    //      PATCH INTERFACE TO IMPLEMENT
    // ----------------------------------------------------------------//
    virtual std::string typeName() const =0;
    virtual std::string typeFamily() const=0;
    
    virtual void updateValues()
    {
    }
    
    void timeChanged(Scalar time)
    {
        SEM_UNUSED(time);
    }

    virtual void read(const iomanagment::Dictionary &fieldDict)
    {
        fieldDict.subDictionary("boundary").subDictionary(name()).entry("value")>>*((baseType*)this);
    }

    virtual void write(iomanagment::Dictionary &fieldDict)
    {
        iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());

        iomanagment::DictEntry *valEntry =new iomanagment::DictEntry("value");
        *valEntry<<*((baseType*)this);

        iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
        bDict->add(typeEntry);
        bDict->add(valEntry);

        if(!fieldDict.hasSubDictionary("boundary"))
            fieldDict.add(new iomanagment::Dictionary("boundary"));

        fieldDict.subDictionary("boundary").add(bDict);
    }


    // ----------------------------------------------------------------//
    //      EXPRESION TEMPLATES HANDLING
    // ----------------------------------------------------------------//
    /// copy constructor from expresion
    template<typename Exp>
    PatchField(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }

    /// assign from single value
    PatchField<T> & operator = (const T& singVal)
    {
        baseType::operator =(singVal);
        return *this;
    }

    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, PatchField<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, PatchField<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, PatchField<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, PatchField<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, PatchField<T>,baseType)

};



}//field
}//SEM

#endif // PATCHFIELD_H
