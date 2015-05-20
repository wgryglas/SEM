#ifndef WEAKDIVW_H
#define WEAKDIVW_H

#include "fields/ElementFieldBase.h"
#include "fields/DiscontinousField.h"
#include "solver/DiscretOperator.h"
#include "utilities/ArrayFunctions.h"
namespace SEM {

    class WeakDivWDiscretOperator : public las::DiscretOperator<Scalar,WeakDivWDiscretOperator>
{
    const field::DiscontinousField<Vector> m_field;
public:
    
    template<typename ElDerived>
    WeakDivWDiscretOperator(const field::ElementFieldBase<Vector,ElDerived> & elField)
    : m_field(elField)
    {
    }
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(WeakDivWDiscretOperator,Scalar)
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &elements, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
        for(size_t e=0; e < elements.size(); ++e)    
        {
            numArray<Scalar> elResult;
            elements[e].weakDivW(m_field.element(e), elResult);
            auto localRhs = rhsVector[0].slice( elements[e].indexVectorMask() );
            assign( elResult, localRhs);
        }
    }
};

// namespace weak {
//     
//     
// }


}

#endif // WEAKDIVW_H
