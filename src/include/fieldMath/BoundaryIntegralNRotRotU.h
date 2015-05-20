
#ifndef BOUNDARYINTEGRALNROTROTU_H
#define BOUNDARYINTEGRALNROTROTU_H

#include "fields/GeometricField.h"
#include "solver/DiscretOperator.h"
#include "mesh/Boundary.h"
#include "mesh/BoundaryEdge.h"
#include "elements/RealElement.h"
#include "fields/PatchField.h"
#include "fieldMath/TimeExtrapolation.h"
#include "fieldMath/BIntegral_ddtUdotN.h"

namespace SEM {

template<typename Acceptor>
class BoundaryIntegralNRotRotU : public las::DiscretOperator<Scalar,BoundaryIntegralNRotRotU<Acceptor> >
{
    const field::GeometricField<Vector> &m_vel;
    const Scalar m_coeff;
    bool m_extrapolate;
    Acceptor m_accept;
public:
    BoundaryIntegralNRotRotU(const field::GeometricField<Vector> &velocity, Acceptor acceptor, Scalar coeff=1., bool extrapolat=false)
    : m_vel(velocity), m_coeff(coeff), m_extrapolate(extrapolat), m_accept(acceptor)
    {
    }
    

    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(BoundaryIntegralNRotRotU,Scalar)
    
    /// \brief interface for evaluating explicite operator. Result is assigned only to rhsVector
    /// \param mesh - mesh associated with solving domain
    /// \param rhsVector - rhsVectors for each field entity component
    /// \param assigner - funcator designed to assigne values into rhsVector -->this allows to assigne directly computded 
    ///                   values into rhsVector, without temprary vector
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
        for( const mesh::Boundary & boundary : mesh.boundaryMesh() )
        {
            if(! m_accept( m_vel.patch( boundary.name() ).typeFamily() ) )
                continue;
            
            for(size_t edge=0; edge<boundary.size(); ++edge)
            {
                Vector normal = boundary.normal(edge);
                const mesh::RealElement & element=boundary.element(edge);
                const std::vector<int> & map = element.indexVectorMask();
                
                numArray<Vector> elementExtrapVelocity(map.size());
                
                if(m_extrapolate)
                {
                    elementExtrapVelocity = extInTime(m_vel, map);
                }
                else
                {
                    elementExtrapVelocity= m_vel.slice(map);
                }
                
                //elementExtrapVelocity = 3.*m_vel.oldField(0).slice(map) -3.*m_vel.oldField(1).slice(map) + m_vel.oldField(2).slice(map);
                

                auto localRhs = rhsVector[0].slice(map);
                assign(m_coeff*element.boundInt_NxRotRotU(elementExtrapVelocity, boundary[edge].edgeId() ),localRhs);
            }
        }
    }
};


namespace weak{
    
boost::shared_ptr<las::DiscretOperator<Scalar,BoundaryIntegralNRotRotU<AcceptAllBoundary> > >    
bInt_rotUdotNxGradPhi(const field::GeometricField<Vector> & velocity,Scalar coeff=1., bool extrapolate=false);
    
boost::shared_ptr<las::DiscretOperator<Scalar,BoundaryIntegralNRotRotU<CheckSingleBoundary> > >
bInt_rotUdotNxGradPhi(const field::GeometricField<Vector> & velocity,const std::string & patchTypeFamily, Scalar coeff=1., bool extrapolate=false );

}//weak


}//SEM

#endif // BOUNDARYINTEGRALNROTROTU_H
