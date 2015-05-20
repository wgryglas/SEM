#ifndef _Equation_H_
#define _Equation_H_

#include <vector>

#include <Eigen/Sparse>

#include "iomanagment/InfoStream.h"
#include "utilities/Reference.h"
#include "utilities/TypeDefs.h"
#include "components/Basic.h"
#include "fields/GeometricField.h"
#include "EquationMatrix.h"
#include "EquationAssignment.h"

namespace SEM { namespace las{
  
  

    /** ******************************************************
     * \class DiscretOperator
     * Base class for all discretizers. It's also
     * implemented by CompoundEquationBuilder which
     * controls the way how elements of matrix shall
     * be assigned to final equation matrix. It's full
     * representation of operations between matrix derived
     * from each parto of equation.
     *********************************************************/
    template<typename T, typename Derived>
    class DiscretOperator
    {
        REFERENCE_TYPE(SINGLE_ARG(DiscretOperator<T,Derived>))
        
    public:
        DiscretOperator(){}

        /// \brief field - due to problems with shared_ptr and template functions
        /// couldn't handle 2 types of discret operators - impicit and explicit. 
        /// Because of this fact we introduce here one generic operator, and below
        /// method shall be implemented only by implicit(or both types) operator. Explicit
        /// operator type shall here throw error, because explicit operators do not
        /// need to be linked with field, and shall never be asked for that.
        /// It's not best solution, but right now there is no other solution to 
        /// serve some operators which would be able to behave once as explicit 
        /// and once as implicit depending on equation == operator side. 
        /// Operators need to be pased as pointers due to inheritance of bese operator class, 
        /// and its life time must be dynamicly controlled(oprator need to be pased between functions).
        /// This allows easy shared_ptr. Shared_ptr has problems with auto upcasting in template functions, so 
        /// each shared_ptr<DiscretOperator<...> > need to be directly constructed
        /// as DiscretOperator<..> but not as some derived type. If we would split
        /// operators into explicit and implicit, then we wont have possibility
        /// to serve explicit and implicit at the same time, becasue there wont be way to auto-cast into 
        /// implicit/explicit base type in template fuction argument--> 
        /// in that situation it would have to be done directly but we want to do 
        /// this implicitly, according to equation side. So this solution we don't
        /// use "multi-level" inheritance, but only 1-level and some generic 
        /// DiscretOperator which handle both cases.
        SEM::field::GeometricField<T> *field() 
        {
            return static_cast<Derived*>(this)->field();
        }
        
       bool isMatrixDiagonal() const
       {
	   static_cast<const Derived*>(this)->isMatrixDiagonal();
       }
        
        template<typename Assigner>
        void buildMatrixAsDiagonal(const mesh::Mesh& mesh, SEM::las::SEMVector & diagonal, SEM::las::SEMVector & rhsVector, const SEM::las::AssigmentBase<Assigner> & assigner) const
        {
	    static_cast<const Derived*>(this)->buildMatrixAsDiagonal(mesh, diagonal, rhsVector, assigner);
	}
        
        /// \brief interface which allows to build eqution directly into matrix and rhsVector.
        /// The way how equation shall be bild is delegated to Derived class, which must implement
        /// the same method declaration as below.
        /// \param matrix - eqution matrix where shall be placed discretization result-->it's collection of elemental matrices
        /// \param rhsVector - rhs vector of equation discretization, where result shall be placed.
        /// \param assigner - funcator which takes 2, any type classes, and and assigns by special operator first arg. to second.
        ///                   assigner is replacement of operation: matrix[...]=.... and rhsVector[...]=... because it provides
        ///                   apropriate operator describing wheter value shall be added, substracted or assigned to matrix.
        ///                   User, who implement below interface, shall writa instead matrix[...]=.... and rhsVector[...]=...
        ///                   somthing like this: assigner(...., matrix[...]) and assigner(...,rhsVector[...])
        template<typename Assigner>
        void buildImplicit(const mesh::Mesh& mesh, SEM::las::SEMMatrix &matrix, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assigner) const
        {
            static_cast<const Derived*>(this)->buildImplicit(mesh,matrix,rhsVector, assigner);
        }
        
        /// \brief interface for evaluating explicite operator. Result is assigned only to rhsVector
        /// \param mesh - mesh associated with solving domain
        /// \param rhsVector - rhsVectors for each field entity component
        /// \param assigner - funcator designed to assigne values into rhsVector -->this allows to assigne directly computded 
        ///                   values into rhsVector, without temprary vector
        template<typename Assigner>
        void buildExplicit(const mesh::Mesh &mesh, SEM::las::SEMVector &rhsVector, const SEM::las::AssigmentBase<Assigner> & assigner) const
        {
            static_cast<const Derived*>(this)->buildExplicit(mesh,rhsVector,assigner);
        }
    };

#define DECLARE_AND_BLOCK_DIAGONAL_TYPE(OPERATOR_NAME)                                      \
    bool isMatrixDiagonal() const {return false;}                                           \
    template<typename Assigner>                                                             \
    inline void buildMatrixAsDiagonal(const mesh::Mesh&, SEM::las::SEMVector & diagonal,    \
					SEM::las::SEMVector & rhsVector,                    \
                                  const SEM::las::AssigmentBase<Assigner> & assigner) const \
    {                                                                                       \
        using namespace SEM::iomanagment;                                                   \
        ErrorInFunction << "Trying to obtain diagonal matrix version, while \n"             \
                        << "OPERATOR_NAME operator is defined as non diagonal"              \
                        << endProgram;                                                      \
    }
    
    
#define DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(OPERATOR_NAME,T_NAME)                 \
    SEM::field::GeometricField<T_NAME>* field()                                             \
    {                                                                                       \
        using namespace SEM::iomanagment;                                                   \
        ErrorInFunction << "Trying to use explicit operator buildas implicit.\n"            \
                        << "OPERATOR_NAME shall not appear on lhs of equation \n"           \
                        << "Wrong discret equation definition \n"                           \
                        << endProgram;                                                      \
        return nullptr;                                                                     \
    }                                                                                       \
    template<typename Assigner>                                                             \
    inline void buildImplicit(const SEM::mesh::Mesh &elements,SEM::las::SEMMatrix &matrix,  \
                              SEM::las::SEMVector &rhsVector,                               \
                              const SEM::las::AssigmentBase<Assigner> & assign) const       \
    {                                                                                       \
        using namespace SEM::iomanagment;                                                   \
        ErrorInFunction << "Trying to use explicit operator as implicit.\n"                 \
        << "OPERATOR_NAME shall not appear on lhs of equation \n"                           \
        << "Wrong discret equation definition \n"                                           \
        << endProgram;                                                                      \
    }                                                                                       \
    DECLARE_AND_BLOCK_DIAGONAL_TYPE(OPERATOR_NAME)
    
  
    // OLD STUFF USING IMPLICIT AND EXPLICIT OPERATORS SEPARETE - BETTER WAY, BUT DON'T RESOLVE SHARED_PTR PROBLEM
    
    //       template<typename T, typename Derived>
    //       class DiscretOperator : public ImpliciteOperator<T,DiscretOperator<T,Derived> >, public ExpliciteOperator<T,DiscretOperator<T,Derived> >
    //       {
    //           REFERENCE_TYPE(DiscretOperator<T COMMA Derived>)
    //           
    //           typedef DiscretOperator<T,Derived> type;
    //           DiscretOperator(const type& other);
    //           type& operator=(const type& other);
    //       public:
    //           DiscretOperator(){}
    //           virtual ~DiscretOperator(){}
    //           
    //           //Implicite/Explicite operators interface
    //           field::GeometricField<T> & solvingField() const
    //           {
    //               static_cast<const Derived*>(this)->solvingField();
    //           }
    //           
    //           template<typename DerivedAssigner>
    //           void build(const mesh::Mesh &mesh, SEMMatrix & matrix, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
    //           {
    //               static_cast<const Derived*>(this)->build(mesh,matrix,rhsVector,assigner);
    //           }
    //           template<typename DerivedAssigner>
    //           void build(const mesh::Mesh &mesh, SEMVector &rhsVector, const AssigmentBase<DerivedAssigner> &assigner) const 
    //           {
    //               static_cast<const Derived*>(this)->build(mesh,rhsVector,assigner);
    //           }
    //       };
    //     
    
    
}//las
}//SEM


#endif //_Equation_H_
struct stat;