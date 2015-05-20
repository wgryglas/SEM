#include "DivergenceOperator.h"

#include "fields/GeometricField.h"

namespace SEM { namespace field {
    
DivergenceOperator::DivergenceOperator(const GeometricField< Vector >& field) 
: baseType(field.mesh()), m_field(field)
{
}

numArray< Scalar > DivergenceOperator::element(size_t e) const 
{
    numArray<Scalar> result;
    mesh()[e].div(m_field.element(e), result);
    return result;
}
    
const GeometricField< Vector >& DivergenceOperator::field() const 
{
    return m_field;
}


DivergenceOperator div(const GeometricField< Vector >& field) 
{
    return DivergenceOperator(field);
}



}//field
}//SEM