// #ifndef EQUATION_PART_H
// #define EQUATION_PART_H
// 
// #include "iomanagment/InfoStream.h"
// #include "fields/GeometricField.h"
// #include "EquationAssignment.h"
// #include "DiscretOperator.h"
// #include "CompoundDiscretOperator.h"
// 
// namespace SEM { namespace las {
// 
//     /** ***********************************************
//      *  \class EquationPair 
//      *  Place holder containing pair of equation and
//      *  field as varibale which shall be solve with 
//      *  this equation
//      * ************************************************/
//     template <typename T, typename DerivedEq>
//     class EquationPart
//     {
//         field::GeometricField<T> & m_field;
//         typename DiscretOperator<T,DerivedEq>::ref m_discreteOperator;
//        
//         /// disallow assigment 
//         EquationPart& operator=(const EquationPart<T,DerivedEq> &other);
//         
//     public:
//         EquationPart(field::GeometricField<T> &f, const typename DiscretOperator<T,DerivedEq>::ref &eBuilder)
//         : m_field(f),m_discreteOperator(eBuilder)
//         {
//         }
//         
//         EquationPart(const EquationPart<T,DerivedEq> &other)
//         : m_field(other.m_field), m_discreteOperator(other.m_discreteOperator)
//         {
//         }
//         
//         field::GeometricField<T> & field() const {return m_field; }
//         typename DiscretOperator<T,DerivedEq>::ref discreteOperator() const {return m_discreteOperator;}
//     };
//     
//     
//     #define MAKE_COMPOUNDEQUATION_PART_OPERATOR(OPERATOR, ASSINGER_NAME)                                \
//     template<typename T, typename Lhs, typename Rhs>                                                    \
//     EquationPart<T,CompoundDiscretOperator<T, Lhs, Rhs, ASSINGER_NAME> >                                \
//     operator OPERATOR                                                                                   \
//     (const EquationPart<T,Lhs > &lhsPair, const EquationPart<T,Rhs> &rhsPair)                           \
//     {                                                                                                   \
//         using namespace iomanagment;                                                                    \
//         if( &lhsPair.field() != &rhsPair.field() )                                                      \
//         {                                                                                               \
//             ErrorInFunction<<"Can't create equation where parts of it are "                             \
//             <<"defined for difrent fields"<<endProgram;                                                 \
//         }                                                                                               \
//         typedef CompoundDiscretOperator<T, Lhs, Rhs, ASSINGER_NAME> CompoundType;                       \
//         typename CompoundType::ref compound(new CompoundType(lhsPair.discreteOperator(),rhsPair.discreteOperator() ) );\
//                                                                                                         \
//         return EquationPart<T,CompoundType>( lhsPair.field(), compound);                                \
//     }                                                                                                   \
//     
//     MAKE_COMPOUNDEQUATION_PART_OPERATOR(+,AddEqAssigment)
//     MAKE_COMPOUNDEQUATION_PART_OPERATOR(-,SubstractEqAssigment)
//     MAKE_COMPOUNDEQUATION_PART_OPERATOR(*,MultiplyEqAssigment)
//     MAKE_COMPOUNDEQUATION_PART_OPERATOR(/,DivideEqAssigment)
//     
//     
// #define MAKE_EQUATION_PART_TYPEDEF(VALUE_NAME, CLASS_TYPE)              \
// public:                                                                 \
//     typedef SEM::las::EquationPart<VALUE_NAME,CLASS_TYPE > equationPart;\
// private:                                                                \
//     
// 
// }//las
// }//SEM
//     
//     #endif //EQUATION_PART_H