#ifndef RHS_BUILDER_H
#define RHS_BUILDER_H

#include "utilities/Reference.h"
#include "utilities/ET_Array.h"
#include "mesh/Mesh.h"
#include "fields/GeometricField.h"
#include "solver/DiscretOperator.h"
#include "components/Basic.h"
#include "components/CmpTraits.h"
#include "utilities/ArrayFunctions.h"
#include "iomanagment/InfoStream.h"

namespace SEM {
    
    /** ************************************************************
     * \class SourceDiscreteOperator
     * Class for creating rhs vector by integrating array of values
     * --------------------------------------------------------------
     * Wrapper for any ET_Array(simple array or expression) implementation.
     * As result of build method it integrates values of given vector
     * (assumed that are nodal values) and places it into rhs-vector
     * of matrix-vector linear equation to solve.
     * ---------------------------------------------------------------
     * Purpose of this class is to allow placing in easy manner some 
     * field of values into rhs vector of SEM equation, eg like below:
     * ddt(T) - laplacia(T) == X*Y-sin(Y),
     * where T - field to solve, and X,Y are some instances of numArrayBase
     * --> so its just an array of values.
     * *************************************************************/
    template<typename T>
    class SourceDiscreteOperator : public SEM::las::DiscretOperator<T,SourceDiscreteOperator<T> >
    {
        typedef las::DiscretOperator<T, SourceDiscreteOperator<T> > baseType;
        numArray<T> m_sourceValues;
        
    public:
        template<typename ExpressionDerived>
        explicit SourceDiscreteOperator(const array::ET_Array<T,ExpressionDerived>& rhsExpr)
        : m_sourceValues(rhsExpr)
        {
        }
        
        DECLARE_AND_BLOCK_IMPLICIT_OPERATOR_FUNCTIONS(SourceDiscreteOperator<T>,T)
        
        template<typename Assigner>
        inline void buildExplicit(const mesh::Mesh &elements, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
        {
            unsigned int dimSize = CmpTraits<T>::dim();
            
            //iterate over field dimmension
            for(unsigned int dim=0; dim<dimSize; ++dim)
            {
                typename CmpTraits<T>::const_CmpArrayRef cmpValues = CmpTraits<T>::cmpArray(m_sourceValues,dim);
                
                //iterate over mesh elements
                for(size_t e=0; e<elements.size(); ++e)
                {
                    // mass matrix non zero components - diagonal components
                    typename numArray2D<Scalar>::const_mappedArray localMassVector = elements[e].massMatrix().sliceArray(elements[e].localNodesInMatrix());
                    
                    typename numArray<Scalar>::indexMapped localRhs = rhsVector[dim].slice(elements[e].indexVectorMask());
                    
                    auto cmpLocalValues = cmpValues.slice(elements[e].indexVectorMask());
                    
                    assign(localMassVector*cmpLocalValues, localRhs);
                }
                
            }
            
        }
    
        
        
        
    };
    
}//SEM


#endif //RHS_BUILDER_H 