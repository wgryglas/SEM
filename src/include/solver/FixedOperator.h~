
#include "boost/shared_ptr.hpp"

#include "DiscretOperator.h"


namespace SEM {
namespace las {

template<typename T, typename Operator>
class FixedOperator : public DiscretOperator<T, FixedOperator<T, Operator> >
{
    boost::shared_ptr<DiscretOperator<T, Operator> > m_operator;
    las::SEMMatrix m_matrix;
    bool build;
    
    void buildMatrix(const mesh::Mesh& mesh, SEM::las::SEMMatrix &matrix)
    {
        Scalar unknownNumber= mesh.spectralNodes().size();
        size_t dim = CmpTraits<T>::dim();
        
        matrix.resize(dim);
        SEM::las::SEMVector rhsVector;
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
        m_operator->buildImplicit(mesh,matrix,rhsVector,assigner);
        
        build = true;
    }
    
public:
    FixedOperator(boost::shared_ptr<DiscretOperator<T,Operator> > op): m_operator(op), build(false)
    {
        buildMatrix(op->field()->mesh(), m_matrix);
    }
    
    SEM::field::GeometricField<T> *field() { return m_operator->field();}
    
    template<typename Assigner>
    void buildImplicit(const mesh::Mesh& mesh, SEM::las::SEMMatrix &matrix, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assigner) const
    {
        for(size_t dim=0; dim<m_matrix.size(); ++dim)
        {
            for(size_t e=0; e<m_matrix[dim].size(); ++e)
            {
                assigner(m_matrix[dim][e], matrix[dim][e]);
            }
        }
    }
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assigner) const
    {
        //TODO
    }
};

template<typename T, typename Operator>
boost::shared_ptr<DiscretOperator<T, FixedOperator<T,Operator> > > fixOperator(boost::shared_ptr<DiscretOperator<T,Operator> > op)
{
       DiscretOperator<T, FixedOperator<T,Operator> > * fOp = new FixedOperator<T,Operator>(op);
       boost::shared_ptr<DiscretOperator<T, FixedOperator<T,Operator> > > shared_ptr_fixedOp(fOp);
       return shared_ptr_fixedOp;
}


}
}