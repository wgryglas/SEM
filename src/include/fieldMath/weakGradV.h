#ifndef WEAKGRADV_H
#define WEAKGRADV_H


#include "solver/DiscretOperator.h"
#include "solver/EquationPart.h"
#include "components/Basic.h"
#include "utilities/VectorUtils.h"
#include "components/CmpTraits.h"
#include "fields/ElementFieldBase.h"
#include "fields/DiscontinousField.h"

namespace SEM { 
/** *************************************************************************
 * \class WeakGradVDiscretOperator
 * Operator to calculate expression: Ipq = integrate( field*Grad(Phi_pq) )
 * where field is defined just by any template array expression. Due to 
 * that fact values of such expression(or explicitly some array deriving from
 * expression) is evaluated(in case of simple array it's unfortuanately copied)
 * to local storage -> due to fact of using expression values must be evaluated,
 * can't be stored as reference.
 ************************************************************************* **/
class WeakGradVDiscretOperator : public SEM::las::DiscretOperator<Vector,WeakGradVDiscretOperator>
{
    typedef las::DiscretOperator<Vector, WeakGradVDiscretOperator> baseType;
    field::DiscontinousField<Scalar> m_field;
    
public:
    template<typename ElementDerived>
    WeakGradVDiscretOperator(const field::ElementFieldBase<Scalar,ElementDerived> & expression):
        m_field(expression)
    {
    }
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(WeakGradVDiscretOperator,Vector)
    
    template<typename Assigner>
    void buildExplicit(const SEM::mesh::Mesh &mesh, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
    {
        for(size_t e=0; e < mesh.size(); ++e) //iterate over elements
        {
            const numArray<Scalar>& localField = m_field.element(e);
            numArray<Vector> localResult;
            mesh[e].weakGradV(localField,localResult);
            
            for(size_t dim=0; dim<CmpTraits<Vector>::dim(); ++dim)
            {
                auto cmpLocalResult = CmpTraits<Vector>::cmpArray(localResult,dim);
                auto cmpLocalRhs = rhsVector[dim].slice(mesh[e].indexVectorMask());
                assign(cmpLocalResult, cmpLocalRhs);
            }
        }
    }
};

}//SEM

#endif // WEAKGRADV_H
