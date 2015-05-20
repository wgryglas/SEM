#ifndef DIVERGENCEOPERATOR_H
#define DIVERGENCEOPERATOR_H

#include "fields/ElementFieldBase.h"
#include "components/Basic.h"

namespace SEM { namespace field {

template<typename T>
class GeometricField;
    
    
class DivergenceOperator : public ElementFieldBase<Scalar,DivergenceOperator>
{
    const GeometricField<Vector> &m_field;
    typedef ElementFieldBase<Scalar,DivergenceOperator> baseType;
public :
    explicit DivergenceOperator(const GeometricField<Vector> & field);
    
    numArray<Scalar> element(size_t e) const;
    
    const GeometricField<Vector> & field() const;
};

DivergenceOperator div(const GeometricField<Vector> & field);


}//field
}//SEM

#endif // DIVERGENCEOPERATOR_H
