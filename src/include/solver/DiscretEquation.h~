
#ifndef DISCRETEQUATION_H
#define DISCRETEQUATION_H

#include <boost/shared_ptr.hpp>

#include "components/CmpTraits.h"
#include "DiscretOperator.h"

namespace SEM { namespace las {

template<typename T, typename Implicit, typename Explicit>
class DiscretEquation
{
    boost::shared_ptr<DiscretOperator<T, Implicit> > m_implicitPart;
    boost::shared_ptr<DiscretOperator<T, Explicit> > m_explicitPart;
    
    SEM::field::GeometricField<T> &m_solveField;
    
public:
    DiscretEquation(boost::shared_ptr<DiscretOperator<T, Implicit> > implicitPart, boost::shared_ptr<DiscretOperator<T, Explicit> > explicitPart):
        m_implicitPart(implicitPart), m_explicitPart(explicitPart), m_solveField(*implicitPart->field())
    {
    }

    SEM::field::GeometricField<T> &field() { return m_solveField;}
    
    // Compiler generated copy,assigne 

    virtual ~DiscretEquation(){}
    
    void buildEquation(SEMMatrix &matrix, SEMVector &rhsVector) const
    {
        const mesh::Mesh &mesh = m_solveField.mesh();
        
        Scalar unknownNumber= mesh.spectralNodes().size();
        size_t dim = CmpTraits<T>::dim();
        
        matrix.resize(dim);
        rhsVector.resize(dim);
        for(size_t d=0; d<dim; ++d)
        {
            matrix[d].resize(mesh.size());
            for(size_t e=0; e<mesh.size(); ++e)   
            {
                size_t s =mesh[e].indexVectorMask().size();
                matrix[d][e].resize2D(s,s,0.);
            }
            rhsVector[d].resize(unknownNumber,0.);
        }
        
        AddEqAssigment assigner;
        m_implicitPart->buildImplicit(mesh,matrix,rhsVector,assigner);
        m_explicitPart->buildExplicit(mesh,rhsVector,assigner);
    }
};
    
template<typename T, typename Implicit, typename Explicit>
DiscretEquation<T, Implicit, Explicit> 
operator ==
(boost::shared_ptr<DiscretOperator<T,Implicit> > implicitPart, boost::shared_ptr< DiscretOperator<T,Explicit> > explicitPart)
{
    return DiscretEquation<T,Implicit, Explicit>(implicitPart,explicitPart);
}


template<typename T> 
class DummyExplicitOperator : public DiscretOperator<T,DummyExplicitOperator<T> >
{
public:
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(DummyExplicitOperator<T>,T)

    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &mesh, SEMVector &rhsVector, const AssigmentBase<Assigner> & assigner) const
    {
       //Do nothing;
    }
    
};

    



}//las
}//SEM

#endif // DISCRETEQUATION_H
