#ifndef POSTOBJECT_H
#define POSTOBJECT_H

#include "components/Basic.h"
#include "utilities/ImplementationFactory.h"
#include "mesh/Mesh.h"
#include "iomanagment/Register.h"
#include "iomanagment/Dictionary2.h"
#include "time/Time.h"

namespace SEM {

/// \class PostObject 
/// base interface for objects that would evaluate
/// value when time would change
class PostObject 
{
    DECLARE_IMPLEMENTATION_FACTORY(PostObject, const iomanagment::Dictionary&, const mesh::Mesh&, Time&)
    
public:
    virtual void calculate()=0;
    
    virtual ~PostObject()
    {
    }
    
};
}
#endif // POSTOBJECT_H
