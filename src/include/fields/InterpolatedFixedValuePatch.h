#ifndef INTERPOLATED_FIXED_VALUE_PATCH_H
#define INTERPOLATED_FIXED_VALUE_PATCH_H

#include "fields/FixedValuePatch.h"
#include "utilities/LinearInterpolator.h"

namespace SEM { namespace field {

template<typename T>
class InterpolatedFixedValuePatch : public FixedValuePatch<T>
{
    typedef FixedValuePatch<T> baseType;
    
    LinearInterpolator<Vector,T> m_interpolator;
    
public:
    static std::string TYPE_NAME;
    
    InterpolatedFixedValuePatch(const mesh::Boundary &edges) 
    : baseType(edges)
    {
    }
    
    InterpolatedFixedValuePatch(const InterpolatedFixedValuePatch<T>& other)
        : baseType(other)
    {
    }
    
    std::string typeName() const {return InterpolatedFixedValuePatch<T>::TYPE_NAME;}
    std::string typeFamily() const {return DIRICHLET;}
    
    // ----------------------------------------------------------------//
    //      EXPRESION TEMPLATES HANDLING
    // ----------------------------------------------------------------//
    /// copy constructor from expresion
    template<typename Exp>
    InterpolatedFixedValuePatch(const array::ET_Array<T,Exp> & exp): baseType(exp)
    {
    }
    
    /// asigment operators from expressions
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, T, InterpolatedFixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, T, InterpolatedFixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, T, InterpolatedFixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, T, InterpolatedFixedValuePatch<T>,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, T, InterpolatedFixedValuePatch<T>,baseType)
    
    using baseType::boundaryEdges;
    using baseType::name;
    virtual void read(const iomanagment::Dictionary &fieldDict)
    {
        m_interpolator.update( fieldDict.subDictionary("boundary").subDictionary( name() ).entry("value") );
        
        for(size_t e=0; e<boundaryEdges().size(); ++e)
        {
                auto bNodes = boundaryEdges().nodes(e);
            
               for(size_t n=0; n<bNodes.size(); ++n)
               {
                    (*this)[e][n] = m_interpolator(bNodes[n]);
               }
        }
    }

    virtual void write(iomanagment::Dictionary &fieldDict)
    {
        iomanagment::DictEntry *typeEntry = new iomanagment::DictEntry("type",typeName());

        iomanagment::DictEntry *valEntry = new iomanagment::DictEntry("value");

        m_interpolator.writeValues(*valEntry);
        
        iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
        bDict->add(typeEntry);
        bDict->add(valEntry);

        if(!fieldDict.hasSubDictionary("boundary"))
            fieldDict.add(new iomanagment::Dictionary("boundary"));

        fieldDict.subDictionary("boundary").add(bDict);
    }

    
    
};
    
    template<typename T>
    std::string InterpolatedFixedValuePatch<T>::TYPE_NAME = "interpolatedFixedValue";

} //field
} //SEM

#endif // INTERPOLATED_FIXED_VALUE_PATCH_H
