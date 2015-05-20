
#ifndef WEAKFORMDEFINITION_H
#define WEAKFORMDEFINITION_H

#include "solver/DiscretOperator.h"
#include "fields/GeometricField.h"
#include "fields/DiscontinousField.h"
#include "fieldMath/TimeDerivativeDiscreteOperator.h"
#include "iomanagment/case.h"
#include "fieldMath/ddt.h"
#include "fieldMath/weakGradV.h"
#include "fieldMath/SourceDiscreteOperator.h"
#include "fieldMath/BoundaryIntegralGeometricField.h"
#include "fieldMath/LaplaceDiscreteOperator.h"
#include "fieldMath/GradientOperator.h"
#include "fieldMath/weakDivW.h"
#include "fieldMath/ConvectiveDerivative.h"

namespace SEM { namespace weak {
 

class phiFunction {};
class grad_phi {};

const phiFunction phi;


grad_phi grad(phiFunction p) { return grad_phi();}

// ----------------- ddt ------------------------------------//

template<typename T>
struct ddtForm : public SEM::array::ET_Array<T,ddtForm<T> >
{
    ddtForm(field::GeometricField<T> & f) 
    : m_field(f), m_dt(Case::time().timeStep()),
      m_coeffs(selectDiscretizationScheme(f)) 
    {
    }
    
    Scalar dt() const {return m_dt;}
    std::vector<Scalar> coeffs() const { return m_coeffs;}
    
    field::GeometricField<T> & field()  { return m_field;}
    
    size_t expSize() const { return m_field.size(); }
    
    T evalAt(size_t node) const 
    {
        T val = m_coeffs[0]/m_dt*m_field[node];
        for(unsigned int i=1; i< m_coeffs.size(); ++i)
        {
            val+= (m_coeffs[i]/m_dt) * m_field.oldField(i-1)[node];
        }
        return val;
    }
    
private:
    field::GeometricField<T> & m_field;
    Scalar m_dt;
    std::vector<Scalar> m_coeffs;
};

template<typename T>
ddtForm<T> ddt(field::GeometricField<T> & field) { return ddtForm<T>(field);}


/// BILINEAR FORM OF DDT AND TEST FUNCTION
template<typename T>
boost::shared_ptr<las::DiscretOperator<T,TimeDerivativeDiscreteOperator<T> > > 
a(ddtForm<T> fieldForm, phiFunction phiFunction)
{
    return boost::shared_ptr<las::DiscretOperator<T,TimeDerivativeDiscreteOperator<T> > >
    (
        new TimeDerivativeDiscreteOperator<T>(fieldForm.field(),fieldForm.dt(), fieldForm.coeffs() )
    );
}

template<typename T>
struct ddnForm 
{
    ddnForm(SEM::field::GeometricField<T> & f): m_field(f)
    {
    }
    
    ddnForm(const ddnForm<T> & other): m_field(other.field()){}
    
    SEM::field::GeometricField<T> & field() const { return m_field;}
private:
    SEM::field::GeometricField<T> & m_field;
};

template<typename T>
ddnForm<T> ddn(field::GeometricField<T> &f)
{
    return ddnForm<T>(f);
}

}//weak    

namespace las{
    
template<typename T>
class BilinearFieldPhi : public DiscretOperator<T,BilinearFieldPhi<T> >
{
    field::DiscontinousField<T> m_field;
    
public:
    template<typename DerivedField>
    BilinearFieldPhi(const field::ElementFieldBase<T,DerivedField> & f)
    : m_field(f)
    {
    }
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &elements, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
    {
        unsigned int dimSize = CmpTraits<T>::dim();
        
        for(size_t e=0; e<m_field.elementsNumber(); ++e)
        {
            auto localMassVector = elements[e].massMatrix().sliceArray(elements[e].localNodesInMatrix());                
            
            for(unsigned int dim=0; dim<dimSize; ++dim)
            {
                auto rhsElementCmp = rhsVector[dim].slice(elements[e].indexVectorMask());

                auto elementFieldCmp = CmpTraits<T>::cmpArray(m_field.element(e),dim);
                
                assign(elementFieldCmp*localMassVector, rhsElementCmp);
            }
        }
    }
    
    DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(BilinearFieldPhi<T>,T)
};

template<typename T> 
class BilinearImplFieldPhi : public DiscretOperator<T,BilinearImplFieldPhi<T> >
{
    field::GeometricField<T> & m_field;
    Scalar m_coeff;
public:
    BilinearImplFieldPhi(field::GeometricField<T> & f, Scalar coeff)
    : m_field(f), m_coeff(coeff)
    {
    }
    
    field::GeometricField<T>* field() const { return &m_field;}
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh &elements, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
    {
        unsigned int dimSize = CmpTraits<T>::dim();
        
        for(size_t e=0; e<elements.size(); ++e)
        {
            numArray2D<Scalar>::const_mappedArray localMassVector(elements[e].massMatrix(), elements[e].localNodesInMatrix());
            
            for(unsigned int dim=0; dim<dimSize; ++dim)
            {
                auto rhsElementCmp = rhsVector[dim].slice(elements[e].indexVectorMask());
                
                auto cmpField = CmpTraits<T>::cmpArray(m_field,dim);
                auto localCmpField = cmpField.slice(elements[e].indexVectorMask());
                
                assign(m_coeff*localCmpField*rhsElementCmp, rhsElementCmp);
            }
        }
    }
    
    template<typename Assigner>
    void buildImplicit(const mesh::Mesh &mesh, las::SEMMatrix & matrix, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> &assign) const 
    {
        const mesh::Mesh & elements = m_field.mesh();
        
        size_t dimSize = CmpTraits<T>::dim();
        
        //Iterate over elements
        for(size_t e=0; e<elements.size(); ++e)
        {
            numArray2D<Scalar>::const_mappedArray massVector(elements[e].massMatrix(), elements[e].localNodesInMatrix());
            
            for(size_t dim =0; dim < dimSize; ++dim)
            {
                numArray2D<Scalar>::diagonalType matrixDiag = matrix[dim][e].diagonal();
                
                //Assign discretization elements to local matrix
                assign( massVector*m_coeff, matrixDiag );
            }
        }
    }
    
};



}//las


namespace weak {

/// BILINEAR FORM OF ANY (EXPLICIT) FIELD AND TEST FUNCTION
template<typename T, typename ElementDerived>
boost::shared_ptr<las::DiscretOperator<T,las::BilinearFieldPhi<T> > > 
a(const field::ElementFieldBase<T,ElementDerived> & f, phiFunction test)
{
    return boost::shared_ptr<las::DiscretOperator<T,las::BilinearFieldPhi<T> > >
            (
                new las::BilinearFieldPhi<T>(f)
            );
}

template<typename T, typename ET_Derived>
boost::shared_ptr<las::DiscretOperator<T,SourceDiscreteOperator<T> > > 
a(const array::ET_Array<T,ET_Derived> & f, phiFunction test)
{
    return boost::shared_ptr<las::DiscretOperator<T,SourceDiscreteOperator<T> > >
    (
        new SourceDiscreteOperator<T>(f)
    );
}


/// BILINEAR FORM OF (IMPLICIT) CONTINIOUS FILED AND TEST FUNCTION
template<typename T>
boost::shared_ptr<las::DiscretOperator<T,las::BilinearImplFieldPhi<T> > > 
a(field::GeometricField<T> & f, phiFunction test)
{
    return boost::shared_ptr<las::DiscretOperator<T,las::BilinearImplFieldPhi<T> > >
    (
        new las::BilinearImplFieldPhi<T>(f,1.)
    );
}

/// BILINEAR FORM OF (IMPLICIT) CONTINIOUS FILED AND TEST FUNCTION
template<typename T>
boost::shared_ptr<las::DiscretOperator<T,las::BilinearImplFieldPhi<T> > > 
a(Scalar mulVal, field::GeometricField<T> & f, phiFunction test)
{
    return boost::shared_ptr<las::DiscretOperator<T,las::BilinearImplFieldPhi<T> > >
    (
        new las::BilinearImplFieldPhi<T>(f,mulVal)
    );
}


/// BILINEAR FORM OF SCALAR FIELD AND GRADIENT OF TEST FUNCTION
template<typename ElementDerived>
boost::shared_ptr<las::DiscretOperator<Vector,WeakGradVDiscretOperator> >
a(const field::ElementFieldBase<Scalar,ElementDerived> & f, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Vector,WeakGradVDiscretOperator> >
    (
        new WeakGradVDiscretOperator(f)
    );
}

inline boost::shared_ptr<las::DiscretOperator<Vector,WeakGradVDiscretOperator> >
a(const field::GeometricField<Scalar> & f, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Vector,WeakGradVDiscretOperator> >
    (
        new WeakGradVDiscretOperator(field::GeometricFieldElementProxy<Scalar>(f))
    );
}

/// BILINEAR FORM OF VECTOR FIELD AND GRADIENT OF TEST FUNCTION
template<typename ElementDerived>
boost::shared_ptr<las::DiscretOperator<Scalar,WeakDivWDiscretOperator> >
a(const field::ElementFieldBase<Vector,ElementDerived> & f, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,WeakDivWDiscretOperator> >
    (
        new WeakDivWDiscretOperator(f)
    );
}

boost::shared_ptr<las::DiscretOperator<Scalar,WeakDivWDiscretOperator> >
inline a(const field::GeometricField<Vector> & f, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,WeakDivWDiscretOperator> >
    (
        new WeakDivWDiscretOperator(field::GeometricFieldElementProxy<Vector>(f))
    );
}



/// SCALAR EQUATION FROM BILINEAR FORM OF FIELD GRADIENT AND TEST FUNCTION GRADIENT
boost::shared_ptr<las::DiscretOperator<Scalar,LaplaceDiscreteOperator<Scalar,Scalar> > >
inline a(field::GradientOperator grad_field, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,LaplaceDiscreteOperator<Scalar,Scalar> > >
    (
        new LaplaceDiscreteOperator<Scalar,Scalar>(*grad_field.field(),1.)
    );
}

boost::shared_ptr<las::DiscretOperator<Scalar,LaplaceDiscreteOperator<Scalar,Scalar> > >
inline a(field::ElementFieldScalarMulLeftOperand<Vector,field::GradientOperator> mul_grad_field, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,LaplaceDiscreteOperator<Scalar,Scalar> > >
    (
        new LaplaceDiscreteOperator<Scalar,Scalar>(*mul_grad_field.field().field(),mul_grad_field.singleValue())
    );
}


/// VECTOR EQUATION FROM BILINEAR FORM OF FIELD GRADIENT AND TEST FUNCTION GRADIENT
boost::shared_ptr<las::DiscretOperator<Vector,LaplaceDiscreteOperator<Vector,Scalar> > >
inline a(field::VectorGradientOperator grad_field, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Vector,LaplaceDiscreteOperator<Vector,Scalar> > >
    (
        new LaplaceDiscreteOperator<Vector,Scalar>(*grad_field.field(),1.)
    );
}


boost::shared_ptr<las::DiscretOperator<Vector,LaplaceDiscreteOperator<Vector,Scalar> > >
inline a(field::ElementFieldScalarMulLeftOperand<Tensor,field::VectorGradientOperator> mul_grad_field, grad_phi test)
{
    return boost::shared_ptr<las::DiscretOperator<Vector,LaplaceDiscreteOperator<Vector,Scalar> > >
    (
        new LaplaceDiscreteOperator<Vector,Scalar>(*mul_grad_field.field().field(),mul_grad_field.singleValue())
    );
}


/// BILINEAR FORM FROM CONVECTIVE DERIVATIVE AND TEST FUNCTION
template<typename T>
inline boost::shared_ptr<las::DiscretOperator<T,field::ConvectiveDerivativeOperator<T> > >
a(const field::ConvectiveDerivative<T> & cDeriv, weak::phiFunction testFun)
{
    return boost::shared_ptr<las::DiscretOperator<T,field::ConvectiveDerivativeOperator<T> > >
    (
        new field::ConvectiveDerivativeOperator<T>(cDeriv.velocity(),cDeriv.field())
    );
}

template<typename T>
inline boost::shared_ptr<las::DiscretOperator<T,field::ConvectiveDerivativeOperator<T> > >
a(field::ElementFieldScalarMulLeftOperand<T,field::ConvectiveDerivative<T> > mul_conv, phiFunction test)
{
    return boost::shared_ptr<las::DiscretOperator<T,field::ConvectiveDerivativeOperator<T> > >
    (
        new field::ConvectiveDerivativeOperator<T>(mul_conv.field().velocity(), mul_conv.field().field(),mul_conv.singleValue())
    );
}


/// BOUNDARY INTEGRAL OF FIELD NORMAL GRADIENT MULTIPLIED BY TEST FUNCTION
template<typename T>
boost::shared_ptr<las::DiscretOperator<T,BoundaryIntegralDDNGeometricFieldPhi<T> > >
bInt(Scalar coeff, ddnForm<T> derivForm, phiFunction testFun)
{
    return boost::shared_ptr<las::DiscretOperator<T,BoundaryIntegralDDNGeometricFieldPhi<T> > >
    (
        new BoundaryIntegralDDNGeometricFieldPhi<T>(derivForm.field(),coeff)
    );
}


}//weak

}//SEM











#endif // WEAKFORMDEFINITION_H

