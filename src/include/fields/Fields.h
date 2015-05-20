#ifndef FIELDS_H
#define FIELDS_H

#include "GeometricField.h"
#include "components/Basic.h"
#include "iomanagment/RegistryObject.h"
#include "iomanagment/case.h"

namespace SEM { namespace field {

typedef  GeometricField<Scalar> ScalarField;
typedef  GeometricField<Vector> VectorField;


}//field
}//SEM





#endif // FIELDS_H
