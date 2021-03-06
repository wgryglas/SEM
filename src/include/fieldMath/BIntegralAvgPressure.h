#ifndef BINTEGRALAVGPRESSURE_H
#define BINTEGRALAVGPRESSURE_H

#include "solver/DiscretOperator.h"

#include "fields/GeometricField.h"

#include "fields/AvgPressure.h"
#include "fields/TCAvgPressure.h"
#include "fields/TCAvgPressure_Stab.h"

namespace SEM {

class BIntegralAvgPressure : public las::DiscretOperator<Vector,BIntegralAvgPressure>
{
    const field::GeometricField<Scalar> &m_pressure;
    const field::GeometricField<Vector> &m_velocity;
    
public:
    BIntegralAvgPressure(const field::GeometricField<Scalar> & pressureField, const field::GeometricField<Vector> & velocityField);
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(BIntegral_ddtUdotN,Vector)
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
        for(const typename field::GeometricField<Scalar>::Patch* p: m_pressure.boundaryFields())
        {
            if( p->typeName()== field::TCAvgPressure::TYPE_NAME || p->typeName()== field::AvgPressure::TYPE_NAME)
            {
                
                const field::AvgPressure* avgPatch = dynamic_cast<const field::AvgPressure*>(p);
                
                if(avgPatch == nullptr)
                {
                    iomanagment::Warning<< "can't cast patch with typeName equal to "<<p->typeName()<<" to "<<field::AvgPressure::TYPE_NAME<<" instance"<<iomanagment::endInfo;
                    continue;
                }
                
                const mesh::Boundary & boundary = p->boundaryEdges();   
                
                Scalar avgPressure = avgPatch->avgPressure();
                
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const numArray<Scalar> &bCoeffs= boundary[e].neumanBCPreValue();
                    const std::vector<int> & edgeMap = boundary[e].gNodesIds();
                    Vector normal = boundary.normal(e);
                    
                    auto xCmpRhs = rhsVector[0].slice(edgeMap);
                    auto yCmpRhs = rhsVector[1].slice(edgeMap);
                    assign( avgPressure*normal.x()*bCoeffs, xCmpRhs);
                    assign( avgPressure*normal.y()*bCoeffs, yCmpRhs);
                    
                    //!!!!!!!! Shoudln't appear stabilization part -nu*div(U) in above equation?
                    
                }
            }
            else if( p->typeName() == field::TCAvgPressure_Stab::TYPE_NAME )
            {
                const field::TCAvgPressure_Stab* avgPatch = dynamic_cast<const field::TCAvgPressure_Stab*>(p);
                
                if(avgPatch == nullptr)
                {
                    iomanagment::Warning<< "can't cast patch with typeName equal to "<<p->typeName()<<" to "<<field::TCAvgPressure_Stab::TYPE_NAME<<" instance"<<iomanagment::endInfo;
                    continue;
                }
                
                const mesh::Boundary & boundary = p->boundaryEdges();   
                
                Scalar avgPressure = avgPatch->avgPressure();
                
                Scalar ref = avgPatch->Uref()*avgPatch->delta();
                
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const numArray<Scalar> &bCoeffs= boundary[e].neumanBCPreValue();
                    const std::vector<int> & edgeMap = boundary[e].gNodesIds();
                    Vector normal = boundary.normal(e);
                    
                    
                    auto xCmpRhs = rhsVector[0].slice(edgeMap);
                    auto yCmpRhs = rhsVector[1].slice(edgeMap);
                    
                    auto lU = m_velocity.slice(edgeMap);
                    
                    using array::tanh;
                    using array::magSqrt;
                    numArray<Scalar> amplitude = avgPressure - 0.25*magSqrt(lU)*( 1.-tanh( dotProd(lU,normal)/ref ) );
                    
                    assign( amplitude*normal.x()*bCoeffs, xCmpRhs);
                    assign( amplitude*normal.y()*bCoeffs, yCmpRhs);
                }
                
                
                
            }
        }
    }

};

boost::shared_ptr<las::DiscretOperator<Vector,BIntegralAvgPressure> >
bInt_avgPressure(const field::GeometricField<Scalar> & pressureField,  const field::GeometricField<Vector> & velocityField);


}//SEM

#endif // BINTEGRALAVGPRESSURE_H
