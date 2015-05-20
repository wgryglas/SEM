/*
#ifndef EXPLICITEOPERATOR_H
#define EXPLICITEOPERATOR_H


#include "mesh/Mesh.h"
#include "EquationMatrix.h"
#include "EquationAssignment.h"

namespace SEM { namespace las {

template<typename T, typename Derived>
class ExpliciteOperator
{
    REFERENCE_TYPE(ExpliciteOperator<T COMMA Derived>)
    
    ExpliciteOperator(const ExpliciteOperator<T,Derived>& other);
    ExpliciteOperator<T,Derived>& operator=(const ExpliciteOperator<T,Derived>& other);
public:
    ExpliciteOperator() {}
    virtual ~ExpliciteOperator(){}

    template<typename DerivedAssigner>
    void build(const mesh::Mesh &mesh, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const
    {
        static_cast<const Derived*>(this)->build(mesh, rhsVector, assigner);
    }
};


}//las 
}//SEM

#endif // EXPLICITEOPERATOR_H

*/
