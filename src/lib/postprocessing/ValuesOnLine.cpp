#include "ValuesOnLine.h"

#include "utilities/ImplementationFactory.h"

namespace SEM 
{
    typedef ValuesOnLine<Scalar> ScalarVeluesOnLine;
    typedef ValuesOnLine<Vector> VectorVeluesOnLine;
    
    REGISTER_IMPLEMENATION(PostObject,ScalarVeluesOnLine,"scalarValuesOnLine")
    REGISTER_IMPLEMENATION(PostObject,VectorVeluesOnLine,"vectorValuesOnLine")
}