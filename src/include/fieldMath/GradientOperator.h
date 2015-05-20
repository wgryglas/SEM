#ifndef GRADIENTOPERATOR_H
#define GRADIENTOPERATOR_H

#include "fields/ElementFieldBase.h"
#include "components/Basic.h"
#include "fields/DiscontinousField.h"
#include "iomanagment/InfoStream.h"

namespace SEM { namespace field {

template<typename T>
class GeometricField;
    
class GradientOfDiscontinousField : public ElementFieldBase<Vector,GradientOfDiscontinousField>
{
    DiscontinousField<Scalar> m_field;
    typedef ElementFieldBase<Vector,GradientOfDiscontinousField> baseType;
public:
    template<typename ElDerived>
    GradientOfDiscontinousField(const ElementFieldBase<Scalar,ElDerived> & f)
    : baseType(f.mesh()), m_field(f)
    {
    }
    
    numArray<Vector> element(size_t e) const;
};

template<typename ElDerived>
GradientOfDiscontinousField grad(const ElementFieldBase<Scalar,ElDerived> & f)
{
    return GradientOfDiscontinousField(f);
}

class GradientOperator : public ElementFieldBase<Vector,GradientOperator>
{
    GeometricField<Scalar> & m_field;
    typedef ElementFieldBase<Vector,GradientOperator> baseType;
public:
    explicit GradientOperator(GeometricField<Scalar> &field);
    
    GeometricField<Scalar>* field() const;
    
    numArray<Vector> element(size_t e) const;
};

GradientOperator grad(GeometricField<Scalar> & field);

class VectorGradientOperator : public ElementFieldBase<Tensor,VectorGradientOperator>
{
    GeometricField<Vector> & m_field;
    typedef ElementFieldBase<Tensor,VectorGradientOperator> baseType;
public:
    explicit VectorGradientOperator(GeometricField<Vector> &f);
    
    GeometricField<Vector> * field() const;
    
    numArray<Tensor> element(size_t e) const 
    { 
        ErrorInFunction <<" explicit Vector gradient calculation not implemented yet" <<iomanagment::endProgram;
        return numArray<Tensor>(0);
    }
    
};

VectorGradientOperator grad(GeometricField<Vector> & field);

}//field
}//SEM





#endif // GRADIENTOPERATOR_H
