// #ifndef IMPLICITEOPERATOR_H
// #define IMPLICITEOPERATOR_H
// 
// #include "mesh/Mesh.h"
// #include "EquationMatrix.h"
// #include "EquationAssignment.h"
// #include "fields/GeometricField.h"
// #include "utilities/Reference.h"
// 
// namespace SEM{ namespace las {
// 
// template<typename T, typename Derived>
// class ImpliciteOperator
// {
//     REFERENCE_TYPE(ImpliciteOperator<T COMMA Derived>)
//     
//     ImpliciteOperator(const ImpliciteOperator<T,Derived>& other);
//     ImpliciteOperator<T, Derived>& operator=(const ImpliciteOperator& other);
// public:
//     ImpliciteOperator()
//     {
//     }
//     
//     virtual ~ImpliciteOperator(){}
// 
//     field::GeometricField<T> & solvingField()const { return static_cast<const Derived*>(this)->solvingField();}
//     
//     template<typename DerivedAssigner>
//     void build(const mesh::Mesh &mesh, SEMMatrix & matrix, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
//     {
//         static_cast<const Derived*>(this)->build(mesh,matrix,rhsVector,assigner);
//     }
//     
// };
// 
// }//las
// }//SEM
// 
// #endif // IMPLICITEOPERATOR_H
