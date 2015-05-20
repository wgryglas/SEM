#ifndef EXRPESSIONVALUEPATCH_H
#define EXRPESSIONVALUEPATCH_H

#include "boost/array.hpp"

#include "PatchField.h"
#include "components/CmpTraits.h"
#include "utilities/DynamicExpression.h"

namespace SEM { namespace field {

    
template<typename T>
class ExpressionValuePatch : public PatchField<T>
{
private:
    typedef PatchField<T> baseType;
    
    boost::array<DynamicExpression,CmpTraits<T>::DIM_SIZE> m_expressions;
    
    boost::array<std::string,CmpTraits<T>::DIM_SIZE> m_exprStr;
    
public:
    static std::string TYPE_NAME;
    using baseType::boundaryEdges;
    using baseType::name;

    ExpressionValuePatch(const mesh::Boundary &edges, boost::array<std::string,CmpTraits<T>::DIM_SIZE> & expresions)
    : baseType(edges), m_exprStr(expresions), m_expressions(build2DExpressions(expresions))
    {
        evaluateExpr();
    }
    
    ExpressionValuePatch(const mesh::Boundary &edges)
        : baseType(edges)
    {
    }

    ExpressionValuePatch(const ExpressionValuePatch<T>& other)
        : baseType(other)
    {
    }

    // ----------------------------------------------------------------//
    //      PATCH INTERFACE
    // ----------------------------------------------------------------//
    
    void timeChanged(Scalar time)
    {
        for(size_t d=0; d< CmpTraits<T>::DIM_SIZE; ++d)
        {
            if( m_expressions[d].isTransient() )
            {
                evaluateExpr();
                break;
            }
        }
     }
    
    std::string typeName() const {return ExpressionValuePatch<T>::TYPE_NAME;}
    std::string typeFamily() const {return DIRICHLET;}
    
    void read(const iomanagment::Dictionary &fieldDict)
    {
        fieldDict.subDictionary("boundary").subDictionary(name()).entry("value")>>m_exprStr;
        
        for(size_t d=0; d<CmpTraits<T>::DIM_SIZE; ++d)
            m_expressions[d].setExpression(m_exprStr[d]);
        
        evaluateExpr();
    }
    
    void evaluateExpr()
    {
        //set current solving time to expression
        for(DynamicExpression & exp : m_expressions)
            exp.setT( Case::time().time());
        
        for(size_t d=0; d<CmpTraits<T>::DIM_SIZE; ++d)
        {
            if(m_expressions[d].isSpatial()) //if spatial then do expression evaluation for each node
            {
                for(size_t e=0; e<boundaryEdges().size(); ++e)
                {
                    auto edgeNodes = boundaryEdges().nodes(e);
                    for(size_t n=0; n<edgeNodes.size(); ++n)
                    {
                        m_expressions[d].setX(edgeNodes[n].x());
                        m_expressions[d].setY(edgeNodes[n].y());
                        
                        CmpTraits<T>::component((*this)[e][n],d) = m_expressions[d].value();
                    }
                }
            }
            else //if not spatial then do expression evaluation only once 
            {
                Scalar val = m_expressions[d].value();
                for(size_t e=0; e<boundaryEdges().size(); ++e)
                    CmpTraits<T>::cmpArray((*this)[e],d) = val;
            }
        }
    }
    
    void write(iomanagment::Dictionary &fieldDict)
    {
        iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
        
        iomanagment::DictEntry *valEntry =new iomanagment::DictEntry("value");
        *valEntry<<m_exprStr;
        
        iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
        bDict->add(typeEntry);
        bDict->add(valEntry);
        
        if(!fieldDict.hasSubDictionary("boundary"))
            fieldDict.add(new iomanagment::Dictionary("boundary"));
        
        fieldDict.subDictionary("boundary").add(bDict);
    }
    
};

template<typename T>
std::string ExpressionValuePatch<T>::TYPE_NAME = "expressionValue";


}//field
}//SEM





#endif // EXRPESSIONVALUEPATCH_H
