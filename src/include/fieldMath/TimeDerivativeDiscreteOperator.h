
#ifndef TIMEDERIVATIVEDISCRETEOPERATOR_H
#define TIMEDERIVATIVEDISCRETEOPERATOR_H

#include <vector>

#include "utilities/Reference.h"
#include "fields/GeometricField.h"
#include "elements/RealElement.h"
#include "mesh/Mesh.h"
#include "iomanagment/InfoStream.h"
#include "solver/DiscretOperator.h"
#include "solver/EquationPart.h"
#include "components/Basic.h"
#include "components/CmpTraits.h"
#include "solver/EquationAssignment.h"

namespace SEM {
/** ****************************************
 *  \class TimeFirstDerivativeBuilder
 *  First time derivative discretizator. 
 *  This class is supposed to place coefficients
 *  into equation matrix/rhsVector for some 
 *  provided GeometricField in method build
 *  (Build method is called externaly, when 
 *  program decide to create equation matrix 
 *  and rhsVector).
 *  -------------------------------------
 *  Implements CRTP-DiscretOperator.
 * *****************************************/
template<typename T>
class TimeDerivativeDiscreteOperator : public las::DiscretOperator<T,TimeDerivativeDiscreteOperator<T> >
{
    std::vector<Scalar> m_coeffs;
    Scalar m_dt;
    SEM::field::GeometricField<T> &m_field;
public:
    
    TimeDerivativeDiscreteOperator(SEM::field::GeometricField<T>&field, Scalar dt, const std::vector<Scalar> &timeCoeffs)
    : m_coeffs(timeCoeffs),m_dt(dt),m_field(field)
    {
    }
    
    field::GeometricField<T> * field() { return &m_field;}
    
    template<typename DerivedAssigner>
    void buildImplicit(const mesh::Mesh &mesh, las::SEMMatrix & matrix, las::SEMVector &rhsVector, const las::AssigmentBase<DerivedAssigner> &assigner) const;
    
    template<typename DerivedAssigner>
    void buildExplicit(const mesh::Mesh &mesh, las::SEMVector &rhsVector, const las::AssigmentBase<DerivedAssigner> &assigner) const;
  
// private:
//     template<typename Assigner>
//     void buildEuler( las::SEMMatrix &matrix, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const;
//     
//     template<typename Assigner>
//     void buildEulerExplicit(las::SEMVector& rhsVector, const las::AssigmentBase<Assigner>& assign) const;
//     
//     template<typename Assigner>
//     void buildWithSelected(las::SEMMatrix &matrix, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const;
//     
//     template<typename Assigner>
//     void buildWithSelectedExplicit(las::SEMVector& rhsVector, const las::AssigmentBase<Assigner>& assign) const;
};

///////////////////////////////////////////////////////////////////////////
//      TimeFirstDerivativeBuilder DEFINITIONS
template<typename T>
template<typename DerivedAssigner>
void TimeDerivativeDiscreteOperator<T>::buildImplicit(const mesh::Mesh &mesh, las::SEMMatrix & matrix, las::SEMVector &rhsVector, const las::AssigmentBase<DerivedAssigner> &assign) const 
{
    using namespace iomanagment;
    
    //Steady state, do nothing
    if(m_coeffs.size()==0)
        return;
    
    //Check if there are at laste 2 coefficients
    if(m_coeffs.size()<2)
        ErrorInFunction<<"Can't discretize first time derivative of fild with less\n"
        <<"then 2 coefficients for time steps"<<endProgram;
    
    //Select scheme to discretize according to number of already known(cached) filed solutions
// //     if(m_field.oldFieldsNumber() >= m_coeffs.size()-1)
// //         buildWithSelected(matrix,rhsVector,assign);
// //     else
// //         buildEuler(matrix,rhsVector, assign);
//     
//     buildWithSelected(matrix,rhsVector, assign);
    
    
    const mesh::Mesh & elements = m_field.mesh();
    
    size_t dimSize = CmpTraits<T>::dim();
    
    //Iterate over elements
    for(size_t e=0; e<elements.size(); ++e)
    {
        numArray2D<Scalar>::const_mappedArray massVector(elements[e].massMatrix(), elements[e].localNodesInMatrix());
        
        //iterate over filed dimmensions
        for(size_t dim =0; dim < dimSize; ++dim)
        {
            numArray2D<Scalar>::diagonalType matrixDiag = matrix[dim][e].diagonal();
            
            //Assign discretization elements to local matrix
            Scalar coeff = m_coeffs[0]/m_dt;
            assign( massVector*coeff, matrixDiag );
            
            //Get element rhs vector
            numArray<Scalar>::indexMapped localRhs(rhsVector[dim], elements[e].indexVectorMask());
            
            //Assign as many previous times as discretization scheme demands to
            for(int i=1; i<m_coeffs.size(); ++i)
            {
                //Assign discretization elements to part of rhs-vector
                coeff = -m_coeffs[i]/m_dt;
                const std::vector<int>& vecMask = elements[e].indexVectorMask();
                assign( massVector*CmpTraits<T>::cmpArray(m_field.oldField(i-1),dim).slice(vecMask)*coeff, localRhs );
                
            }
        }
        
    }
    
}

template<typename T>
template<typename DerivedAssigner>
void TimeDerivativeDiscreteOperator<T>::buildExplicit(const mesh::Mesh &mesh, las::SEMVector &rhsVector, const las::AssigmentBase<DerivedAssigner> &assign) const
{
    using namespace iomanagment;
    
    //Steady state, do nothing
    if(m_coeffs.size()==0)
        return;
    
    //Check if there are at laste 2 coefficients
    if(m_coeffs.size()<2)
        ErrorInFunction<<"Can't discretize first time derivative of fild with less\n"
        <<"then 2 coefficients for time steps"<<endProgram;
    
    //Select scheme to discretize according to number of already known(cached) filed solutions
//     if(m_field.oldFieldsNumber() >= m_coeffs.size()-1)
//         buildWithSelectedExplicit(rhsVector,assign);
//     else
//         buildEulerExplicit(rhsVector, assign);
    
    const mesh::Mesh & elements = m_field.mesh();
    
    size_t dimSize = CmpTraits<T>::dim();
    
    //Iterate over elements
    for(size_t dim =0; dim < dimSize; ++dim)
    {
        numArray<Scalar> timeDeriv =m_coeffs[0]*CmpTraits<T>::cmpArray(m_field,dim);
        for(size_t c=1; c<m_coeffs.size(); ++c)
        {
            timeDeriv += m_coeffs[c]*CmpTraits<T>::cmpArray(m_field.cachedField(c-1),dim);
        }
        
        for(size_t e=0; e<elements.size(); ++e)
        {
            numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
            //Assign discretization elements to rhs vector
            numArray<Scalar>::indexMapped localRhs(rhsVector[dim], elements[e].indexVectorMask());
            assign(massVector*timeDeriv.slice(elements[e].indexVectorMask()), localRhs);
        }
    }
    
}


// template<typename T>
// template<typename Assigner>
// void TimeDerivativeDiscreteOperator<T>::buildEuler(las::SEMMatrix &matrix, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
// {
//     const mesh::Mesh & elements = m_field.mesh();
//     
//     size_t dimSize = CmpTraits<T>::dim();
//     
//     //Iterate over elements
//     for(size_t e=0; e<elements.size(); ++e)
//     {
//         numArray2D<Scalar>::const_mappedArray massVector(elements[e].massMatrix(), elements[e].localNodesInMatrix());
//         
//         for(size_t dim =0; dim < dimSize; ++dim)
//         {
//             numArray2D<Scalar>::diagonalType matrixDiag = matrix[dim][e].diagonal();
//             
//             //Assign discretization elements to local matrix
//             assign( massVector*(1./m_dt), matrixDiag );
//             
//             //Assign discretization elements to rhs vector
//             const std::vector<int>& vecMask = elements[e].indexVectorMask();
//             
//             auto localRhs = rhsVector[dim].slice(vecMask);
//             auto fieldCmp = CmpTraits<T>::cmpArray(m_field.oldField(0),dim);
//             auto localFieldCmp = fieldCmp.slice(vecMask);
//             assign(massVector*localFieldCmp*(1./m_dt), localRhs);
//         }
//     }
//     
// }

// template<typename T>
// template<typename Assigner>
// void TimeDerivativeDiscreteOperator<T>::buildEulerExplicit(las::SEMVector& rhsVector, const las::AssigmentBase<Assigner>& assign) const
// {
//     const mesh::Mesh & elements = m_field.mesh();
//     
//     size_t dimSize = CmpTraits<T>::dim();
//     
//     //Iterate over elements
//     for(size_t dim =0; dim < dimSize; ++dim)
//     {
//         numArray<Scalar> timeDeriv = ( CmpTraits<T>::cmpArray(m_field,dim) - CmpTraits<T>::cmpArray(m_field.oldField(0),dim) )/m_dt;
//         
//         for(size_t e=0; e<elements.size(); ++e)
//         {
//             numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
//             //Assign discretization elements to rhs vector
//             numArray<Scalar>::indexMapped localRhs(rhsVector[dim], elements[e].indexVectorMask());
//             assign(massVector*timeDeriv.slice(elements[e].indexVectorMask()), localRhs);
//         }
//     }
// }

// template<typename T>
// template<typename Assigner>
// void TimeDerivativeDiscreteOperator<T>::buildWithSelected(las::SEMMatrix &matrix, las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> & assign) const
// {
//     const mesh::Mesh & elements = m_field.mesh();
//     
//     size_t dimSize = CmpTraits<T>::dim();
//     
//     //Iterate over elements
//     for(size_t e=0; e<elements.size(); ++e)
//     {
//         numArray2D<Scalar>::const_mappedArray massVector(elements[e].massMatrix(), elements[e].localNodesInMatrix());
//         
//         //iterate over filed dimmensions
//         for(size_t dim =0; dim < dimSize; ++dim)
//         {
//             numArray2D<Scalar>::diagonalType matrixDiag = matrix[dim][e].diagonal();
//             
//             //Assign discretization elements to local matrix
//             Scalar coeff = m_coeffs[0]/m_dt;
//             assign( massVector*coeff, matrixDiag );
//             
//             //Get element rhs vector
//             numArray<Scalar>::indexMapped localRhs(rhsVector[dim], elements[e].indexVectorMask());
//             
//             //Assign as many previous times as discretization scheme demands to
//             for(int i=1; i<m_coeffs.size(); ++i)
//             {
//                 //Assign discretization elements to part of rhs-vector
//                 coeff = -m_coeffs[i]/m_dt;
//                 const std::vector<int>& vecMask = elements[e].indexVectorMask();
//                 assign( massVector*CmpTraits<T>::cmpArray(m_field.oldField(i-1),dim).slice(vecMask)*coeff, localRhs );
//                 
//             }
//         }
//         
//     }
// }

// template<typename T>
// template<typename Assigner>
// void TimeDerivativeDiscreteOperator<T>::buildWithSelectedExplicit(las::SEMVector &rhsVector, const las::AssigmentBase<Assigner> &assign) const
// {
//     const mesh::Mesh & elements = m_field.mesh();
//     
//     size_t dimSize = CmpTraits<T>::dim();
//     
//     //Iterate over elements
//     for(size_t dim =0; dim < dimSize; ++dim)
//     {
//         numArray<Scalar> timeDeriv =m_coeffs[0]*CmpTraits<T>::cmpArray(m_field,dim);
//         for(size_t c=1; c<m_coeffs.size(); ++c)
//         {
//             timeDeriv += m_coeffs[c]*CmpTraits<T>::cmpArray(m_field.cachedField(c-1),dim);
//         }
//         
//         for(size_t e=0; e<elements.size(); ++e)
//         {
//             numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
//             //Assign discretization elements to rhs vector
//             numArray<Scalar>::indexMapped localRhs(rhsVector[dim], elements[e].indexVectorMask());
//             assign(massVector*timeDeriv.slice(elements[e].indexVectorMask()), localRhs);
//         }
//     }
// }


}//SEM



#endif // TIMEDERIVATIVEDISCRETEOPERATOR_H
