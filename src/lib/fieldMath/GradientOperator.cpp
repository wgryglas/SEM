#include "GradientOperator.h"

#include "fields/GeometricField.h"

namespace SEM { namespace field {
 
    
//-------------------- discontinous field grad ---------------------------//    
numArray< Vector > GradientOfDiscontinousField::element(size_t e) const 
{
    numArray<Vector> result;
    mesh()[e].grad(m_field.element(e),result);
    return result;
}    

//-------------------- continous field grad ---------------------------//    
GradientOperator::GradientOperator(SEM::field::GeometricField< SEM::Scalar >& f) 
: baseType(f.mesh()), m_field(f)
{
}

numArray< Vector > GradientOperator::element(size_t e) const 
{
    numArray<Vector> result;
    mesh()[e].grad(m_field.element(e),result);
    return result;
}

GeometricField< Scalar > * GradientOperator::field() const{return &m_field; }

GradientOperator grad(GeometricField< Scalar >& f) 
{
    GradientOperator g(f);
    return g;
}

//-------------------- vector field grad ---------------------------//
VectorGradientOperator::VectorGradientOperator(GeometricField< Vector >& f) 
:baseType(f.mesh()), m_field(f)
{
}

GeometricField< Vector > * VectorGradientOperator::field() const
{
    return &m_field;
}

VectorGradientOperator grad(GeometricField< Vector >& field) 
{
    return VectorGradientOperator(field);
}




}//field
}//SEM