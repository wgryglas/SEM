
#ifndef BINTEGRAL_DDTUDOTN_H
#define BINTEGRAL_DDTUDOTN_H


#include <vector>
#include <string>
#include <algorithm>
#include <functional>

#include "solver/DiscretOperator.h"
#include "fields/GeometricField.h"
#include "fieldMath/ddt.h"
#include "mesh/Mesh.h"
#include "utilities/Utilities.h"
#include "utilities/VectorUtils.h"

namespace SEM {


struct CheckSingleBoundary
{
    CheckSingleBoundary(const std::string & family);
    
    bool operator()(const std::string & toCheck) const;
private:
    std::string m_family;
};

struct AcceptAllBoundary
{
    bool operator()(const std::string & toCheck) const;
};

    
template<typename BoundaryAcceptor>
class BIntegral_ddtUdotN : public las::DiscretOperator<Scalar,BIntegral_ddtUdotN<BoundaryAcceptor> >
{
    const field::GeometricField<Vector> & m_field;
    BoundaryAcceptor m_accept;
public:
    
    BIntegral_ddtUdotN(const field::GeometricField<Vector> & field, BoundaryAcceptor acceptor)
    : m_field(field), m_accept(acceptor)
    {
    }
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(BIntegral_ddtUdotN,Vector)
  
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
//         std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(Case::solutionControl().timeDiscretization());
        
//         if(m_field.oldFieldsNumber() >= timeCoeffs.size())
//             buildSelected(mesh,rhsVector, timeCoeffs, assign);
//         else
//             buildEuler(mesh,rhsVector,assign);
//         buildSelected(mesh,rhsVector,timeCoeffs,assign);
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(m_field);
        
        Scalar dt = Case::time().timeStep();
        
//         for(const typename field::GeometricField<Vector>::Patch* p : m_field.boundaryFields())
//         {
//             if(m_accept(p->typeFamily()))
//             {
//                 const mesh::Boundary & boundary = p->boundaryEdges();
//                 
//                 for(size_t e=0; e<boundary.size(); ++e)
//                 {
//                     Vector normal = boundary.normal(e);
//                     
//                     const std::vector<int> & edgeMap = boundary[e].gNodesIds();
//                     
//                     //Calculate edge ddt(U)*normal
//                     numArray<Vector> ddtField = (timeCoeffs[0]/dt)*m_field.slice(edgeMap);
//                     for(unsigned int c=1; c<timeCoeffs.size(); ++c)
//                     {
//                         ddtField += (timeCoeffs[c]/dt) * m_field.oldField(c-1).slice(edgeMap);
//                     }
//                     
//                     // note: value*neumanBCPreValue  is equal to bInt(value,phi_pq)
//                     numArray<Scalar> weakValues = array::dotProd(ddtField,normal)*boundary[e].neumanBCPreValue();
//                     
//                     // Assign (ddt(U)*normal,phi) into edge nodes coresponding equations in rhs vector
//                     auto edgeRhs = rhsVector[0].slice(edgeMap);
//                     assign(weakValues, edgeRhs);
//                 }
//             }
//         }
        
        for(const typename field::GeometricField<Vector>::Patch* p : m_field.boundaryFields())
        {
            if(m_accept(p->typeFamily()))
            {
                const mesh::Boundary & boundary = p->boundaryEdges();
                
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    Vector normal = boundary.normal(e);
                    
                    const std::vector<int> & edgeMap = boundary[e].gNodesIds();
                    const std::vector<double> & bcCoeffs=boundary[e].neumanBCPreValue();
                    
                    for(size_t n=0; n<edgeMap.size(); ++n)
                    {
                        Scalar value = dotProd(normal,timeCoeffs[0]/dt*m_field[edgeMap[n]]);
                        
                        for(unsigned int c=1; c<timeCoeffs.size(); ++c)
                        {
                            value += dotProd(normal, timeCoeffs[c]/dt * m_field.oldField(c-1)[edgeMap[n]] );
                        }
                        
                        value *= bcCoeffs[n];
                        
                        assign(value,rhsVector[0][edgeMap[n]]);
                    }
                }
            }
        }
        
    }
    
// private:
//     template<typename Assigner>
//     void buildEuler(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
//     {
//         Scalar dt = Case::time().timeStep();
//         
//         for(const typename field::GeometricField<Vector>::Patch* p : m_field.boundaryFields())
//         {
//             if(m_accept(p->typeFamily()))
//             {
//                 const mesh::Boundary & boundary = p->boundaryEdges();
//                 
//                 for(size_t e=0; e<boundary.size(); ++e)
//                 {
//                     Vector normal = boundary.normal(e);
//                     const std::vector<int> & edgeMap = boundary[e].gNodesIds();
//                     
//                     // note: value*neumanBCPreValue  is equal to bInt(value,phi_pq)
//                     
//                     numArray<Scalar> weakValues = array::dotProd(( m_field.slice(edgeMap) - m_field.oldField(0).slice(edgeMap) )/dt,normal)*boundary[e].neumanBCPreValue();
//                     
//                     // Assign (ddt(U)*normal,phi) into edge nodes coresponding equations in rhs vector
//                     auto edgeRhs = rhsVector[0].slice(edgeMap);
//                     assign(weakValues, edgeRhs);
//                 }
//             }
//         }
//     }
//     
//     template<typename Assigner>
//     void buildSelected(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const std::vector<Scalar> & timeCoeffs, const SEM::las::AssigmentBase<Assigner> & assign) const
//     {
//         Scalar dt = Case::time().timeStep();
//         
//         for(const typename field::GeometricField<Vector>::Patch* p : m_field.boundaryFields())
//         {
//             if(m_accept(p->typeFamily()))
//             {
//                 const mesh::Boundary & boundary = p->boundaryEdges();
//                 
//                 for(size_t e=0; e<boundary.size(); ++e)
//                 {
//                     Vector normal = boundary.normal(e);
//                     const std::vector<int> & edgeMap = boundary[e].gNodesIds();
//                     
//                     //Calculate edge ddt(U)*normal
//                     numArray<Vector> ddtField = (timeCoeffs[0]/dt)*m_field.slice(edgeMap);
//                     for(unsigned int c=1; c<timeCoeffs.size(); ++c)
//                     {
//                         ddtField += (timeCoeffs[c]/dt) * m_field.oldField(c-1).slice(edgeMap);
//                     }
//                     
//                     // note: value*neumanBCPreValue  is equal to bInt(value,phi_pq)
//                     numArray<Scalar> weakValues = array::dotProd(ddtField,normal)*boundary[e].neumanBCPreValue();
//                     
//                     // Assign (ddt(U)*normal,phi) into edge nodes coresponding equations in rhs vector
//                     auto edgeRhs = rhsVector[0].slice(edgeMap);
//                     assign(weakValues, edgeRhs);
//                 }
//             }
//         } 
//     }
    
};

boost::shared_ptr<las::DiscretOperator<Scalar,BIntegral_ddtUdotN<AcceptAllBoundary> > >
bInt_ddtUdotN(const field::GeometricField<Vector> &field);

boost::shared_ptr<las::DiscretOperator<Scalar,BIntegral_ddtUdotN<CheckSingleBoundary> > >
bInt_ddtUdotN(const field::GeometricField<Vector> &field, const std::string & boundaryFamily);




}//SEM

#endif // BINTEGRAL_DDTUDOTN_H
