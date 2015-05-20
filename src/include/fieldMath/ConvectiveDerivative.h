#ifndef CONVECTIVEDERIVATIVE_H
#define CONVECTIVEDERIVATIVE_H


#include "fields/ElementFieldBase.h"
#include "fields/GeometricField.h"
#include "utilities/ArrayFunctions.h"
#include "solver/DiscretOperator.h"

namespace SEM { namespace field {

template<typename T>
class ConvectiveDerivative : public ElementFieldBase<T,ConvectiveDerivative<T> >
{
    const GeometricField<Vector> & m_velocity;
    GeometricField<T> & m_field;
    typedef ElementFieldBase<T,ConvectiveDerivative<T> > baseType;
public:
    ConvectiveDerivative(const GeometricField<Vector> &velocity,GeometricField<T> & field)
    : baseType(field.mesh()), m_velocity(velocity), m_field(field)
    {
    }
    
    const GeometricField<Vector> & velocity() const { return m_velocity;}
    
    GeometricField<T> & field() const { return m_field;}
    
    numArray<T> element(size_t e) const 
    {
        numArray2D<Scalar> derivMatrix = m_field.mesh()[e].convDerivMatrix(m_velocity.element(e));
        
        numArray<T> result(derivMatrix.size());
        const std::vector<int> & elMap = m_field.mesh()[e].indexVectorMask();
        for(size_t d=0; d<CmpTraits<T>::dim(); ++d )
        {
            CmpTraits<T>::cmpArray(result,d) = array::matrixMul( derivMatrix, CmpTraits<T>::cmpArray(m_field,d).slice(elMap) );
        }
        return result;
    }
};

template<typename T>
class ConvectiveDerivativeOperator : public las::DiscretOperator<T,ConvectiveDerivativeOperator<T> >
{
    const GeometricField<Vector> & m_velocity;
    GeometricField<T> & m_field;
    Scalar m_coeff;
public:
    ConvectiveDerivativeOperator(const GeometricField<Vector> &velocity,GeometricField<T> & field, Scalar coeff=1.)
    : m_velocity(velocity), m_field(field), m_coeff(coeff)
    {
    }
    
    
    SEM::field::GeometricField<T> *field() 
    {
        return &m_field;
    }
    
    const GeometricField<Vector> & velocity() { return m_velocity;}
    
    
    template<typename Assigner>
    void buildImplicit(const mesh::Mesh& elements, SEM::las::SEMMatrix &matrix, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
        for(size_t e=0; e< elements.size(); ++e)
        {
            numArray2D<Scalar> derivMat = elements[e].convDerivMatrix(m_velocity.element(e));
            auto massVec = elements[e].massMatrix().sliceArray(elements[e].localNodesInMatrix());
            
            //calculated weak form
            for(size_t n=0; n < derivMat.size(); ++n)
            {
                derivMat[n] *=m_coeff*massVec[n];
            }

            //assign derivMatrix into matrices in each dimmension
            for(size_t d=0; d<CmpTraits<T>::dim(); ++d)
            {
                assign(derivMat,matrix[d][e]);
            }
        }
    }
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &elements, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assign) const
    {
        for(size_t e=0; e< elements.size(); ++e)
        {
            numArray2D<Scalar> derivMat = elements[e].convDerivMatrix(m_velocity.element(e));
            auto massVec = elements[e].massMatrix().sliceArray(elements[e].localNodesInMatrix());
            
            //calculated weak form
            for(size_t n=0; n < derivMat.size(); ++n)
            {
                derivMat[n] *=m_coeff * massVec[n];
            }
            
            const std::vector<int> & elMap = elements[e].indexVectorMask();
            
            //assign derivMatrix into matrices in each dimmension
            for(size_t d=0; d<CmpTraits<T>::dim(); ++d)
            {
                auto localRhs=rhsVector[d].slice(elMap);
                
                auto cmpField=CmpTraits<T>::cmpArray(m_field,d).slice(elMap);
                
                assign(array::matrixMul(derivMat,cmpField),localRhs);
            }
        }
    }
    
};


template<typename T>
ConvectiveDerivative<T> cDeriv(const GeometricField<Vector> &velocity,GeometricField<T> & field)
{
    return ConvectiveDerivative<T>(velocity,field);
}


}//field
}//SEM

#endif // CONVECTIVEDERIVATIVE_H
