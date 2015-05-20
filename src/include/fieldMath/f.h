
#ifndef _f_H_
#define _f_H_

#include "solver/CompoundDiscretOperator.h"
#include "solver/EquationPart.h"
#include "solver/EquationAssignment.h"
#include "SourceDiscreteOperator.h"
#include "solver/DiscretEquation.h"
namespace SEM{


//----------------------- rhs discret operators functions --------------------------------------------------------
    
/** **************************************************
 * \file f
 *  Set of functions and operators which would 
 *  put into equation explicite operator - SourceDiscreteOperator
 *  --> operators which aim is to fill rhs vector with 
 *  values(or integral of those values) of some array 
 *  expression
 * ****************************************************/


/** **************************************************
 * \brief f - function creating default impl. for explicite operator
 *  from expression template
 * ****************************************************/
template<typename T, typename RhsArray>
inline typename SourceDiscreteOperator<T>::ref f(const array::ET_Array<T,RhsArray>& rhsValues)
 {
     return typename SourceDiscreteOperator<T>::ref
     (
         new SourceDiscreteOperator<T>(rhsValues)
     );
 }

 /** **************************************************
  * \brief operator == - operator which converts expression template
  * into SourceDiscreteOperator and assigns this operator into full
  * discret equation containing lhs implicit part and rhs explicit part
  * ****************************************************/ 
 template<typename T, typename ExprArray, typename ImplicitOp>
 las::DiscretEquation<T,ImplicitOp,SourceDiscreteOperator<T> >
 operator ==
 (boost::shared_ptr<las::DiscretOperator<T,ImplicitOp> > implOperator,const array::ET_Array<T,ExprArray> &rhsValues)
 {
     typename SourceDiscreteOperator<T>::ref explOp( new SourceDiscreteOperator<T>(rhsValues) );
     return las::DiscretEquation<T,ImplicitOp,SourceDiscreteOperator<T> >(implOperator, explOp);
 }

#define SEM_ALLOW_OPERATION_BETWEEN_ETARRAY_AND_DISCRET_OPERATOR(OP,ASSIGNER)                                               \
template<typename T, typename DerivedOp, typename ExprArray>                                                                \
boost::shared_ptr                                                                                                           \
<SEM::las::DiscretOperator                                                                                                  \
<T,SEM::las::CompoundDiscretOperator                                                                                        \
<T,DerivedOp,SourceDiscreteOperator<T>,ASSIGNER> > >                                                                        \
operator OP                                                                                                                  \
(boost::shared_ptr<las::DiscretOperator<T,DerivedOp> > otherOper,                                                           \
const array::ET_Array<T,ExprArray> &arrayValues)                                                                            \
 {                                                                                                                          \
    typename SourceDiscreteOperator<T>::ref arrayOper(new SourceDiscreteOperator<T>(arrayValues));                          \
    typedef SEM::las::CompoundDiscretOperator<T,DerivedOp,SourceDiscreteOperator<T>,ASSIGNER> compoundType; \
    return  boost::shared_ptr<SEM::las::DiscretOperator<T,compoundType> >(new compoundType(otherOper,arrayOper));           \
}                                                                                                                           \
template<typename T, typename DerivedOp, typename ExprArray>                                                                \
boost::shared_ptr                                                                                                           \
<SEM::las::DiscretOperator                                                                                                  \
<T,SEM::las::CompoundDiscretOperator                                                                                        \
<T,DerivedOp,SourceDiscreteOperator<T>,ASSIGNER> > >                                                        \
operator OP                                                                                                                  \
(const array::ET_Array<T,ExprArray> &arrayValues,                                                                           \
boost::shared_ptr<las::DiscretOperator<T,DerivedOp> > otherOper)                                                           \
{                                                                                                                           \
    typename SourceDiscreteOperator<T>::ref arrayOper(new SourceDiscreteOperator<T>(arrayValues));                          \
    typedef SEM::las::CompoundDiscretOperator<T,DerivedOp,SourceDiscreteOperator<T>,ASSIGNER> compoundType; \
    return  boost::shared_ptr<SEM::las::DiscretOperator<T,compoundType> >(new compoundType(otherOper,arrayOper));           \
}                                                                                                                           \

SEM_ALLOW_OPERATION_BETWEEN_ETARRAY_AND_DISCRET_OPERATOR(+,SEM::las::AddEqAssigment)
SEM_ALLOW_OPERATION_BETWEEN_ETARRAY_AND_DISCRET_OPERATOR(-,SEM::las::SubstractEqAssigment)
 
 
 
//  //------------------------------- operator == handling -------------------------------------------------------
//  /** ****************************************************************
//   *  \brief operator==
//   *  Operator for making CompoundDiscretOperator from some DiscretOperator
//   *  lhs and array of known values on rhs<--- operator conversion to
//   *  source part of equation.
//   *  Usage of operator "==" was the only way to do this (operation are 
//   *  done on pointers due to problems with rvalue) because we couldn't
//   *  overload "=" on DiscretOperator (should not return shared_ptr
//   *  from this pointer---> in this situation it may be possible if we would
//   *  be sure that this of DiscretOperator is allocated on heap, what in 
//   *  current impl. should happen, but it's too far dengerous in future
//   *  to use so stiff approach).
//   * *****************************************************************/
//  template<typename T, typename Lhs, typename RhsArray>
//  inline typename las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment>::ref
//  operator ==
//  (const boost::shared_ptr<las::DiscretOperator<T, Lhs> > &lhs, const array::ET_Array<T,RhsArray> & rhsValues )
//  {
//      typename SourceDiscreteOperator<T>::ref rhsRef(new SourceDiscreteOperator<T>(rhsValues));
//      
//      typedef las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment> compoundType;
//      
//      return typename compoundType::ref( new compoundType(lhs,rhsRef) );
//  }
// 
//  /** ***********************************************************
//   * \brief operator ==
//   *  Alow assigment from result of f(some_array) to compound 
//   *  operator 
//   * ************************************************************/
//  template<typename T, typename Lhs>
//  inline typename las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment >::ref
//  operator ==
//  (const boost::shared_ptr<las::DiscretOperator<T, Lhs> > &lhs, const boost::shared_ptr<las::DiscretOperator<T,SourceDiscreteOperator<T> > > &rhs)
//  {
//      typedef las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment> compoundType;
//      return typename compoundType::ref(new compoundType(lhs,rhs));
//  }
//  
//  
//  /** ***********************************************************
//   * \brief operator ==
//   *  Alow assigment from IntegratedSourceOperator into compound 
//   *  operator
//   * ************************************************************/ 
//   template<typename T, typename Lhs>
//   inline typename las::CompoundDiscretOperator<T,Lhs,IntegratedSourceOperator<T>, las::AddEqAssigment>::ref
//   operator ==
//   (boost::shared_ptr<las::DiscretOperator<T,Lhs> > lhs, boost::shared_ptr<las::DiscretOperator<T,IntegratedSourceOperator<T> > > rhs)
//   {
//       typedef las::CompoundDiscretOperator<T,Lhs,IntegratedSourceOperator<T>,las::AddEqAssigment> compoundType;
//       return typename compoundType::ref(new compoundType(lhs,rhs));
//   }
//  
//  /** ***********************************************************
//   * \brief operator ==
//   *  Alow assigment of array values to equation part
//   * ************************************************************/
//  template<typename T, typename Lhs, typename RhsArray>
//  inline las::EquationPart<T,las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment > >
//  operator ==
//  ( const las::EquationPart<T,Lhs> &lhsPart, const array::ET_Array<T,RhsArray> & rhsValues)
//  {
//      typename SourceDiscreteOperator<T>::ref rhsRef(new SourceDiscreteOperator<T>(rhsValues));
//      
//      typedef las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment> compoundType; 
//      
//      return las::EquationPart<T,compoundType>
//      (lhsPart.field(), typename compoundType::ref( new compoundType(lhsPart.discreteOperator(),rhsRef) ) );
//  }
// 
//  /** ***********************************************************
//   * \brief operator ==
//   *  Alow assigment of from result of f(...)
//   * ************************************************************/
//  template<typename T, typename Lhs>
//  las::EquationPart<T,las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment > >
//  operator ==
//  (las::EquationPart<T,Lhs> lhsPart, const boost::shared_ptr<las::DiscretOperator<T,SourceDiscreteOperator<T> > > &rhs)
//  {
//      typedef las::CompoundDiscretOperator<T,Lhs,SourceDiscreteOperator<T>,las::AddEqAssigment> compoundType; 
//      return las::EquationPart<T,compoundType>
//      ( lhsPart.field(), typename compoundType::ref( new compoundType(lhsPart.discreteOperator(), rhs) ) );
//  }
//  
//  /** ***********************************************************
//   * \brief operator ==
//   *  Alow assigment to lhs equation the rhs IntegratedSourceOperator
//   * ************************************************************/
//  template<typename T, typename Lhs>
//  las::EquationPart<T,las::CompoundDiscretOperator<T,Lhs,IntegratedSourceOperator<T>,las::AddEqAssigment > >
//  operator ==
//  (las::EquationPart<T,Lhs> lhsPart, const boost::shared_ptr<las::DiscretOperator<T,IntegratedSourceOperator<T> > > &rhs)
//  {
//      typedef las::CompoundDiscretOperator<T,Lhs,IntegratedSourceOperator<T>,las::AddEqAssigment> compoundType; 
//      return las::EquationPart<T,compoundType>
//      ( lhsPart.field(), typename compoundType::ref( new compoundType(lhsPart.discreteOperator(), rhs) ) );
//  }
 
}//SEM












//#include <Eigen\Dense>
//#include "ScalarField.h"
//#include "SolutionControl.h"
//#include "Equation.h"

//using namespace SEM::FIELDS;

//namespace SEM
//{
//namespace FIELDMATH
//{
//	Eigen::VectorXd& f(ScalarField T);

//}
//}

//#include "solver/Equation.h"
//#include "fields/GenericField.h"

//namespace SEM {

//class SourceBuilder : public las::EquationBuilder<Scalar,SourceBuilder >
//{
//    typedef las::EquationBuilder<Scalar,SourceBuilder > baseType;
//public:
//    using baseType::field;

//    LaplaceBuilder(const Coeff &coefficient, field::GeometricField<Scalar> & field)
//        : baseType(field),m_coeff(coefficient)
//    {
//    }

//    template<typename Assigner>
//    void build(las::SEMMatrix &matrix, numArray<Scalar> &rhsVector, const las::AssigmentBase<Assigner> & assigner);

//private:
//    const Coeff & m_coeff;
//};



//}//SEM


#endif //_f_H_
