#ifndef COMPOUND_EQUATION_BUILDER

#include <boost/shared_ptr.hpp>

#include "utilities/Reference.h"
#include "utilities/TypeDefs.h"
#include "EquationMatrix.h"
#include "EquationAssignment.h"
#include "DiscretOperator.h"
#include "iomanagment/InfoStream.h"
#include <boost/fusion/adapted/std_pair.hpp>

class A;
namespace SEM { namespace las {
/** ****************************************************************************
 * \class CompoundDiscretOperator
 * \brief Place holder for parts of eqution.
 * when build mehod is called on it, then it
 * delegates job to rhs and lhs.
 * This class allows to lazy evaluation of whole equation,
 * while each part of equation evaluates it's job on the
 * same matrix and vector, so there will be no temporary
 * memory allocation.
 * When operator eg. + is found between 2 EquationBuilders
 * implementation, then + assigner is passed to rhs part of equation
 * build method.This tells to rhs EquationBuilder that all elements build by this shall be
 * added to final matrix. Rhs equation builder is evalauated
 * acording to that how someone called build method in below class.
 * -------------------------------------------------------------------------------
 * Final equation is a mixture of CompoundEquationBuilders and direct implementation
 * of field discretization. Final equation matrix and vector are filled from
 * the moste left discretizer to most right,
 * because as you can se in this class build implementation
 * first is evaluated build on lhs builder, and then evaluated build on rhs.
 * To start building procedure zero matirx and vector and EqAssigment as assigner
 * shall be pased to this class method. EqAssigment simply assignes values evaluated
 * by it's discretizer.
 * --------------------------------------------------------------------------------
 * \note Making compound equation is performed by ...::ref {aka boost::sher_ptr<...>
 * due to allow using rvalue while building equation. If it would operate on 
 * const references then there would be problem with scope in which temporary object
 * -rvalue created, would be valid.
 ********************************************************************************/

template<typename T, typename Lhs, typename Rhs, typename RhsAssigner>
class CompoundDiscretOperator : public DiscretOperator<T,CompoundDiscretOperator<T,Lhs,Rhs,RhsAssigner> >
{
    typedef DiscretOperator<T,CompoundDiscretOperator<T,Lhs,Rhs,RhsAssigner> > baseType;
    
    typename DiscretOperator<T, Lhs>::ref m_lhs;
    typename DiscretOperator<T, Rhs>::ref m_rhs;
    RhsAssigner m_rhsAssigner;
    
public:
    CompoundDiscretOperator(const typename DiscretOperator<T, Lhs>::ref &lhs, const typename DiscretOperator<T, Rhs>::ref &rhs)
    : m_lhs(lhs), m_rhs(rhs)
    {
    }
    
    SEM::field::GeometricField<T>* field()
    {
        if(m_lhs->field() != m_rhs->field() )
        {
            using namespace SEM::iomanagment;
            ErrorInFunction <<"Can't create compound implicit operator from 2 operators acting on diffrent fields"<<endProgram;
        }
        
        return m_lhs->field();
    }
    
    bool isMatrixDiagonal() const
    {
        return m_lhs->isMatrixDiagonal() && m_rhs->isMatrixDiagonal();
    }
    
    template<typename Assigner>
    void buildMatrixAsDiagonal(const mesh::Mesh& mesh, SEMVector &diagonal, SEM::las::SEMVector & rhsVector, const AssigmentBase<Assigner> & lhsAssigner) const
    {
        m_lhs->buildMatrixAsDiagonal(mesh, diagonal, rhsVector, lhsAssigner);
	m_rhs->buildMatrixAsDiagonal(mesh, diagonal, rhsVector, m_rhsAssigner);
    }
    
    
    template<typename Assigner>
    void buildImplicit(const mesh::Mesh& mesh, SEMMatrix &matrix, SEMVector &rhsVector, const AssigmentBase<Assigner> & lhsAssigner) const
    {
        m_lhs->buildImplicit(mesh,matrix,rhsVector,lhsAssigner);
        m_rhs->buildImplicit(mesh,matrix,rhsVector,m_rhsAssigner);
    }
    
    template<typename Assigner>
    void buildExplicit(const mesh::Mesh& mesh, SEMVector &rhsVector, const AssigmentBase<Assigner> & lhsAssigner) const
    {
        m_lhs->buildExplicit(mesh,rhsVector,lhsAssigner);
        m_rhs->buildExplicit(mesh,rhsVector,m_rhsAssigner);
    }
};

/** ****************************************************************************
 * \brief preparing operator functions which takes as argument
 * 2 template EquationBuilders (discretizer or some compound type)
 * and returns new CompoundDiscretOperator--> builder representing
 * operation between 2 parts of equation. Needed for lazy matrix and rhs vector
 * filling with only one memory allocation.
 * ---------------------------------------------------------------------------
 * \note input arg. are of type DiscretOperator<...>::ref-shared_ptr, 
 * and return type is DiscretOperator<...>::ref type as well(it's derived from, but
 * shared_ptr can handle such a upcasting) --> used shared_ptr due to fact 
 * described in class CompoundDiscretOperator description.
 ********************************************************************************/
template<typename T, typename Lhs, typename Rhs>
boost::shared_ptr< DiscretOperator< T, CompoundDiscretOperator<T,Lhs,Rhs,SEM::las::AddEqAssigment > > >
operator +
(const boost::shared_ptr<DiscretOperator<T, Lhs> > &lhs, const boost::shared_ptr<DiscretOperator<T, Rhs> > &rhs)
{
    typedef CompoundDiscretOperator<T,Lhs,Rhs,SEM::las::AddEqAssigment> compoundType;
    return boost::shared_ptr<DiscretOperator<T,compoundType> >( new compoundType(lhs,rhs) );
}


template<typename T, typename Lhs, typename Rhs>
boost::shared_ptr< DiscretOperator< T,CompoundDiscretOperator<T,Lhs,Rhs,SEM::las::SubstractEqAssigment > > >
operator -
(const boost::shared_ptr<DiscretOperator<T, Lhs> > &lhs, const boost::shared_ptr<DiscretOperator<T, Rhs> > &rhs)
{
    typedef CompoundDiscretOperator<T,Lhs,Rhs,SEM::las::SubstractEqAssigment> compoundType;
    return boost::shared_ptr<DiscretOperator<T,compoundType> >( new compoundType(lhs,rhs) );
}



//     template<typename T, typename Lhs, typename Rhs, typename RhsAssigner>
//     class ImpliciteCompoundOperator : public ImpliciteOperator<T, ImpliciteCompoundOperator<T,Lhs,Rhs,RhsAssigner> >
//     {
//         typedef ImpliciteCompoundOperator<T,Lhs,Rhs,RhsAssigner> type;
//         
//         ImpliciteCompoundOperator(const type &other);
//         type & operator =(const type &other);
//         
//         typename ImpliciteOperator<T,Lhs>::ref m_lhs;
//         typename ImpliciteOperator<T,Rhs>::ref m_rhs;
//         AssigmentBase<RhsAssigner> m_rhsAssigner;
//         
//     public:
//         ImpliciteCompoundOperator(typename ImpliciteOperator<T,Lhs>::ref lhs, typename ImpliciteOperator<T,Rhs>::ref rhs, AssigmentBase<RhsAssigner> rhsAssigner):
//             m_lhs(lhs), m_rhs(rhs), m_rhsAssigner(rhsAssigner)
//         {
//         }
//         
//         const field::GeometricField<T> & solvingField()const
//         {
//             return m_lhs->solvingField();
//         }
//         
//         template<typename DerivedAssigner>
//         void build(const mesh::Mesh &mesh, SEMMatrix & matrix, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
//         {
//             m_lhs->build(mesh,matrix,rhsVector,assigner);
//             m_rhs->build(mesh,matrix,rhsVector,m_rhsAssigner);
//         }        
//     };
//     
//     template<typename T,typename Lhs, typename Rhs, typename RhsAssigner>
//     class ExpliciteCompundOperator : public ExpliciteOperator<T,ExpliciteCompundOperator<T,Lhs,Rhs,RhsAssigner> >
//     {
//         typedef ExpliciteCompundOperator<T,Lhs,Rhs,RhsAssigner> type;
//         
//         ExpliciteCompundOperator(const type &other);
//         type & operator =(const type &other);
//         
//         typename ExpliciteOperator<T,Lhs>::ref m_lhs;
//         typename ExpliciteOperator<T,Rhs>::ref m_rhs;
//         AssigmentBase<RhsAssigner> m_rhsAssigner;
//         
//     public:
//         ExpliciteCompundOperator(typename ExpliciteOperator<T,Lhs>::ref lhs,typename ExpliciteOperator<T,Rhs>::ref rhs, AssigmentBase<RhsAssigner> rhsAssigner):
//         m_lhs(lhs), m_rhs(rhs), m_rhsAssigner(rhsAssigner)
//         {
//         }
//         
//         template<typename DerivedAssigner>
//         void build(const mesh::Mesh &mesh, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
//         {
//             m_lhs->build(mesh,rhsVector,assigner);
//             m_rhs->build(mesh,rhsVector,m_rhsAssigner);
//         }
//     };
//     
//     template<typename T, typename Lhs, typename Rhs, typename RhsAssigner>
//     class DiscretCompundOperator : public DiscretOperator<T, DiscretCompundOperator<T,Lhs,Rhs,RhsAssigner> >
//     {
//         typedef DiscretCompundOperator<T,Lhs,Rhs,RhsAssigner> type;
//         
//         DiscretCompundOperator(const type &other);
//         type & operator =(const type &other);
//         
//         typename DiscretOperator<T,Lhs>::ref m_lhs;
//         typename DiscretOperator<T,Rhs>::ref m_rhs;
//         AssigmentBase<RhsAssigner> m_rhsAssigner;
//         
//     public:
//         DiscretCompundOperator(typename DiscretOperator<T,Lhs>::ref lhs,typename DiscretOperator<T,Rhs>::ref rhs):
//         m_lhs(lhs), m_rhs(rhs)
//         {
//         }
//         
//         const field::GeometricField<T> & solvingField()const
//         {
//             return m_lhs->solvingField();
//         }
//         
//         template<typename DerivedAssigner>
//         void build(const mesh::Mesh &mesh, SEMMatrix & matrix, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
//         {
//             m_lhs->build(mesh,matrix,rhsVector,assigner);
//             m_rhs->build(mesh,matrix,rhsVector,m_rhsAssigner);
//         }
//         
//         template<typename DerivedAssigner>
//         void build(const mesh::Mesh &mesh, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
//         {
//             m_lhs->build(mesh,rhsVector,assigner);
//             m_rhs->build(mesh,rhsVector,m_rhsAssigner);
//         }
//     };
//     
//     
//     /** ************************************************************************************
//      * \brief macro MAKE_COMPOUND_IMPLICIT_OPERATORS_MATH_SUPPORT
//      * Generator for math operators wich would handle operation between 2 implicit type
//      * operators. Result of such a operator is compond_implicit_operator
//      * *************************************************************************************/
// #define MAKE_COMPOUND_IMPLICIT_OPERATORS_MATH_SUPPORT(OPERATOR, ASSIGNER)                               \
//     template<typename T, typename Lhs, typename Rhs>                                                    \
//     typename ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>::ref                                         \
//     operator OPERATOR                                                                                   \
//     (boost::shared_ptr<ImpliciteOperator<T,Lhs> > lhs, boost::shared_ptr<ImpliciteOperator<T,Lhs> > rhs)\
//     {                                                                                                   \
//         if(lhs->solvingField()!=rhs->solvingField())                                                    \
//         {                                                                                               \
//             using namespace iomanagment;                                                                \
//             ErrorInFunction<<"Can't build equation when lhs operators \n"                               \
//                            <<"acts on diffrent fileds ("                                                \
//                            <<lhs->solvingField().name()<<", "                                           \
//                            <<rhs->solvingField().name()<<")"                                            \
//                            <<endProgram;                                                                \
//         }                                                                                               \
//         return boost::make_shared<ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER> >( lhs,rhs,ASSIGNER() );\
//     }                                                                                                   \
//     
//     /** ************************************************************************************
//      * \brief macro MAKE_COMPOUND_EXPLICIT_OPERATORS_MATH_SUPPORT
//      * Generator for math operators wich would handle operation between 2 explicit type
//      * operators. Result of such a operator is compond_explicit_operator
//      * *************************************************************************************/
// #define MAKE_COMPOUND_EXPLICIT_OPERATORS_MATH_SUPPORT(OPERATOR, ASSIGNER)                               \
//     template<typename T, typename Lhs, typename Rhs>                                                    \
//     typename ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                                          \
//     operator OPERATOR                                                                                   \
//     (boost::shared_ptr<ExpliciteOperator<T,Lhs> > lhs, boost::shared_ptr<ExpliciteOperator<T,Rhs> > rhs)\
//     {                                                                                                   \
//         return boost::make_shared<ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER> >(lhs,rhs,ASSIGNER());   \
//     }                                                                                                   \
//     
//     /** ************************************************************************************
//      * \brief macro MAKE_COMPOUND_DISCRET_OPERATORS_MATH_SUPPORT
//      * Generator for math operators wich would handle operation between 2 dual type 
//      * operators(called disctet_operator). Result of such a operator is compond_discret_operator
//      * *************************************************************************************/
// #define MAKE_COMPOUND_DISCRET_OPERATORS_MATH_SUPPORT(OPERATOR, ASSIGNER)                    \
//     template<typename T,typename Lhs, typename Rhs>                                         \
//     typename DiscretCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                                \
//     operator OPERATOR                                                                       \
//     (typename DiscretOperator<T,Lhs>::ref lhs, typename DiscretOperator<T,Rhs>::ref rhs)    \
//     {                                                                                       \
//         return typename DiscretCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                     \
//         (                                                                                   \
//             new DiscretCompundOperator<T,Lhs,Rhs,ASSIGNER>(lhs,rhs)                         \
//         );                                                                                  \
//     }                                                                                       \
//     
//     /** ************************************************************************************
//      * \brief macro MAKE_COMPOUND_DISCRET_TO_EXPLICIT_MATH_OPERATOR_CONVERSION
//      * Generator for operators wich would handle mathematic between discret(dual type operator)
//      * and implicit operator. As result returned is compound_implicit operator
//      * This handles situation when imlicit operator meet dual type operatore, and choose
//      * to convert dual into implicit
//      * *************************************************************************************/
// #define MAKE_COMPOUND_DISCRET_TO_IMPICIT_MATH_OPERATOR_CONVERSION(OPERATOR, ASSIGNER)       \
//     template<typename T,typename Lhs, typename Rhs>                                         \
//     typename ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>::ref                             \
//     operator OPERATOR                                                                       \
//     (typename ImpliciteOperator<T,Lhs>::ref lhs, typename DiscretOperator<T,Rhs>::ref rhs)  \
//     {                                                                                       \
//         if(lhs->solvingField()!=rhs->solvingField())                                        \
//         {                                                                                   \
//             using namespace iomanagment;                                                    \
//             ErrorInFunction<<"Can't build equation when lhs operators \n"                   \
//             <<"acts on diffrent fileds ("                                                   \
//             <<lhs->solvingField().name()<<", "                                              \
//             <<rhs->solvingField().name()<<")"                                               \
//             <<endProgram;                                                                   \
//         }                                                                                   \
//         return typename ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>::ref                  \
//         (                                                                                   \
//         new ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>(lhs,rhs)                          \
//         );                                                                                  \
//     }                                                                                       \
//     template<typename T,typename Lhs, typename Rhs>                                         \
//     typename ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>::ref                             \
//     operator OPERATOR                                                                       \
//     (typename DiscretOperator<T,Lhs>::ref lhs, typename ImpliciteOperator<T,Rhs>::ref rhs)  \
//     {                                                                                       \
//         if(lhs->solvingField()!=rhs->solvingField())                                        \
//         {                                                                                   \
//             using namespace iomanagment;                                                    \
//             ErrorInFunction<<"Can't build equation when lhs operators \n"                   \
//             <<"acts on diffrent fileds ("                                                   \
//             <<lhs->solvingField().name()<<", "                                              \
//             <<rhs->solvingField().name()<<")"                                               \
//             <<endProgram;                                                                   \
//         }                                                                                   \
//         return typename ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>::ref                  \
//         (                                                                                   \
//             new ImpliciteCompoundOperator<T,Lhs,Rhs,ASSIGNER>(lhs,rhs)                      \
//         );                                                                                  \
//     }                                                                                       \
//     
//     /** ************************************************************************************
//      * \brief macro MAKE_COMPOUND_DISCRET_TO_EXPLICIT_MATH_OPERATOR_CONVERSION
//      * Generator for operators wich would handle mathematic between discret(dual type operator)
//      * and explicit operator. As result returned is compound_explicit operator
//      * This handles situation when explicit operator meet dual type operatore, and choose
//      * to convert dual into explicit
//      * *************************************************************************************/
// #define MAKE_COMPOUND_DISCRET_TO_EXPLICIT_MATH_OPERATOR_CONVERSION(OPERATOR, ASSIGNER)      \
//     template<typename T,typename Lhs, typename Rhs>                                         \
//     typename ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                              \
//     operator OPERATOR                                                                       \
//     (typename ExpliciteOperator<T,Lhs>::ref lhs, typename DiscretOperator<T,Rhs>::ref rhs)  \
//     {                                                                                       \
//         return typename ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                   \
//         (                                                                                   \
//             new ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>(lhs,rhs)                       \
//         );                                                                                  \
//     }                                                                                       \
//     template<typename T,typename Lhs, typename Rhs>                                         \
//     typename ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                              \
//     operator OPERATOR                                                                       \
//     (typename DiscretOperator<T,Lhs>::ref lhs, typename ExpliciteOperator<T,Rhs>::ref rhs)  \
//     {                                                                                       \
//         return typename ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>::ref                   \
//         (                                                                                   \
//             new ExpliciteCompundOperator<T,Lhs,Rhs,ASSIGNER>(lhs,rhs)                       \
//         );                                                                                  \
//     }                                                                                       \
//     
//     
//     
//     template<typename T,typename Lhs, typename Rhs>                                         
//     typename DiscretCompundOperator<T,Lhs,Rhs,AddEqAssigment>::ref                                
//     operator +
//     (boost::shared_ptr<DiscretOperator<T,Lhs> > lhs, boost::shared_ptr<DiscretOperator<T,Rhs> > rhs)    
//     {                                                                                       
//         return typename DiscretCompundOperator<T,Lhs,Rhs,AddEqAssigment>::ref                    
//         (                                                                                  
//             new DiscretCompundOperator<T,Lhs,Rhs,AddEqAssigment>(lhs,rhs)                        
//         );                                                                                 
//     }                                                                                       
//     template<typename T,typename Lhs, typename Rhs>                                         
//     typename DiscretCompundOperator<T,Lhs,Rhs,SubstractEqAssigment>::ref                                
//     operator -
//     (boost::shared_ptr<DiscretOperator<T,Lhs> > lhs, boost::shared_ptr<DiscretOperator<T,Rhs> > rhs)    
//     {                                                                                       
//         return typename DiscretCompundOperator<T,Lhs,Rhs,SubstractEqAssigment>::ref                    
//         (                                                                                  
//             new DiscretCompundOperator<T,Lhs,Rhs,SubstractEqAssigment>(lhs,rhs)                        
//         );                                                                                 
//     }
//     
//     
//     
//     /** ************************************************************************************
//     * \brief macro MAKE_COMPOUND_OPERATORS_SUPPORT collection of all necessery functions
//     * to fully support discret operators math - explicit, implicit and dual(called simply 
//     * discrete).
//     * *************************************************************************************/
//     #define MAKE_COMPOUND_OPERATORS_SUPPORT(OPERATOR,ASSIGNER)                              \
//     MAKE_COMPOUND_IMPLICIT_OPERATORS_MATH_SUPPORT(OPERATOR,ASSIGNER)                        \
// 
//     //MAKE_COMPOUND_EXPLICIT_OPERATORS_MATH_SUPPORT(OPERATOR,ASSIGNER)                        \
//     //MAKE_COMPOUND_DISCRET_TO_IMPICIT_MATH_OPERATOR_CONVERSION(OPERATOR,ASSIGNER)            \
//     //MAKE_COMPOUND_DISCRET_TO_EXPLICIT_MATH_OPERATOR_CONVERSION(OPERATOR,ASSIGNER)           \
//     //MAKE_COMPOUND_DISCRET_OPERATORS_MATH_SUPPORT(OPERATOR,ASSIGNER)                         \
//     
//     
//     MAKE_COMPOUND_OPERATORS_SUPPORT(+,AddEqAssigment)                          
//     MAKE_COMPOUND_OPERATORS_SUPPORT(-,SubstractEqAssigment)                    
        
    
}//las
}//SEM

#endif //COMPOUND_EQUATION_BUILDER