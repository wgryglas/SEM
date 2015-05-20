
#ifndef ROTOPERATOR_H
#define ROTOPERATOR_H

#include "fields/ElementFieldBase.h"
#include "fields/GeometricField.h"

namespace SEM { namespace field {

class RotOperator : public ElementFieldBase<Scalar,RotOperator >
{
    const GeometricField<Vector> & m_field;
    typedef ElementFieldBase<Scalar,RotOperator> baseType;
public :
    using baseType::mesh;
    explicit RotOperator(const GeometricField<Vector> & f)
    : baseType(f.mesh()), m_field(f)
    {
    }
    
    numArray<Scalar> element(size_t e) const
    {
        return mesh()[e].rot(m_field.element(e));
    }
};
}


inline field::RotOperator rot(const field::GeometricField<Vector> & f)
{
    return field::RotOperator(f);
}



}

#endif // ROTOPERATOR_H
