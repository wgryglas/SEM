
#ifndef _Solver_H_
#define _Solver_H_

#include "Eigen/Dense"
#include "Eigen/Sparse"

#include "time/Timmer.h"
#include "utilities/Utilities.h"
#include "utilities/VectorUtils.h"
#include "utilities/numArray.h"
#include "utilities/ArrayFunctions.h"
#include "mesh/Mesh.h"

#include "DiscretOperator.h"
#include "EquationPart.h"
#include "EquationMatrix.h"
#include "numArrayEigen.h"
#include "DiscretEquation.h"

#include "components/Basic.h"
#include "components/CmpTraits.h"

#include "iomanagment/InfoStream.h"
#include "iomanagment/case.h"
#include "iomanagment/SolutionControl.h"

#include "fieldMath/WeakFormDefinition.h"

namespace SEM{ namespace las {

    /** ********************************************************************
     * \class Solver - base interface for runtime chose solver
     * *********************************************************************/
    template<typename T>
    struct Solver
    {
        virtual void solveEquation(field::GeometricField<T> &field, SEMMatrix & matrix, SEMVector & rhsVector) const= 0;
    };

    
    /** ********************************************************************
     * \class EigenLinearSolver - base interface for runtime chose one of
     *                            eigen linear solvers
     * *********************************************************************/
    struct EigenLinearSolver 
    {
        typedef Eigen::SparseMatrix<Scalar> EigenMatrix;
        typedef numArrayEigen<Scalar> EigenArray;
        typedef typename EigenArray::EigenData EigenVector;
        

        virtual void solve(EigenMatrix &matrix, EigenVector &rhsVector, EigenVector &solution) const =0;
    };
    
    
    // forward decalaration of used below function (definition under below class)
    template<typename T>
    void applyDirichletBCInMatrixVector(field::GeometricField<T> &solField, SEMMatrix &matrix, SEMVector &rhsVector);
    
    
    /** ********************************************************************
     * \class EigenSolver - EigenLinearSolver wrapper allowing to use it with
     *                      SEM fields and mesh interface
     * *********************************************************************/
    template<typename T>
    struct EigenSolver : public Solver<T>
    {
      typedef Eigen::Triplet<Scalar> Coeff;
      typedef typename EigenLinearSolver::EigenMatrix EigenMatrix;
      typedef typename EigenLinearSolver::EigenVector EigenVector;
      typedef typename EigenLinearSolver::EigenArray EigenArray;

      EigenSolver(const EigenLinearSolver* linearSolver) : m_linearSolver(linearSolver) {}
      ~EigenSolver() { delete m_linearSolver; }
      
      
      void applyDirichletValsAndCreateDirichletInfo
      (
            field::GeometricField<T> & solField,
            numArray<bool> & dirMask,
            numArray<int> & vecToNonDirVecMap,
            numArray<int> & nonDirVecToVecMap
      ) const;
      
      template<typename CmpDerived>
      void assemblyReducedMatrixVector
      (
            const numArray<bool> & dirMask,
            const numArray<int> & vecToNonDirVecMap,
            const numArray<int> & nonDirVecToVecMap,
            const mesh::Mesh &elements,
            const numArrayBase<CmpDerived> &cmpField,
            ElMatrix & matrix,
            ElVector & rhsVector,
            EigenMatrix &e_mat,
            EigenArray &e_rhs
      ) const;
      
      
      void assemblyMatrixVector
      (
              const mesh::Mesh & mesh,
              ElMatrix & matrix,
              ElVector & rhsVector,
              EigenMatrix &e_mat,
              EigenArray &e_rhs
      ) const;
      
      
      void assemblyMatrixVectorApplyDirichletBC
      (
            const field::GeometricField<T> & solField,
            size_t dim,
            const mesh::Mesh & mesh,
            ElMatrix & matrix,
            ElVector & rhsVector,
            EigenMatrix &e_mat,
            EigenArray &e_rhs
      )const;
      
      
      void solveEquation(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &rhsVector) const;

      
    private:
        const EigenLinearSolver* m_linearSolver;
    };
    
    template<typename T>
    void EigenSolver<T>::applyDirichletValsAndCreateDirichletInfo
    (
        field::GeometricField<T> & solField,
        numArray<bool> & dirMask,
        numArray<int> & vecToNonDirVecMap,
        numArray<int> & nonDirVecToVecMap
    ) const
    {
        //Make dirichlet mask and apply dirichlet values to current solution
        dirMask.resize(solField.size());
        dirMask = false;
        for(typename field::GeometricField<T>::Patch* patch : solField.boundaryFields() )
        {
            //check if it's exacly dirichlet BC
            if(patch->typeFamily()==field::DIRICHLET)
            {
                const mesh::Boundary & boundary = patch->boundaryEdges();
                
                //iterate over egdes
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const mesh::BoundaryEdge & edge = boundary[e];
                    solField.slice(edge.gNodesIds()) = (*patch)[e];;
                    dirMask.slice(edge.gNodesIds()) = true;
                }
            }
        }
        
        size_t dirNumber=0;
        for(const bool & flag : dirMask)
        {
            if(flag) ++dirNumber;
        }
        
        //create information about mapping from/to new linear equation system
        nonDirVecToVecMap.resize(solField.size()-dirNumber);
        vecToNonDirVecMap.resize(solField.size()); 
        size_t curNonDir=0;
        for(size_t n=0; n<dirMask.size(); ++n)
        {
            if(!dirMask[n])
            {
                vecToNonDirVecMap[n] = curNonDir;
                nonDirVecToVecMap[curNonDir] = n;
                ++curNonDir;
            }
        }
    }
    
    template<typename T>
    template<typename CmpDerived>
    void EigenSolver<T>::assemblyReducedMatrixVector
    (
        const numArray<bool> & dirMask,
        const numArray<int> & vecToNonDirVecMap,
        const numArray<int> & nonDirVecToVecMap,
        const mesh::Mesh &elements,
        const numArrayBase<CmpDerived> &cmpField,
        ElMatrix & matrix,
        ElVector & rhsVector,
        EigenMatrix &e_mat,
        EigenArray &e_rhs
    ) const
    {
        using namespace iomanagment;
        InfoArrow<<"Assembling local equations into global matrix"<<std::endl;
        
        e_rhs = rhsVector.slice(nonDirVecToVecMap);
        
        //estimate number of sparse matrix entries
        size_t tripleSize=0;
        for(size_t i=0; i<matrix.size(); ++i)
            tripleSize+=matrix[i].expSize(); //expSize is equal to the number of elements in array2D
        
        //reserve space for sparse matrix entries
        std::vector<Coeff> triples;
        triples.reserve(tripleSize);            
         
        
        //Build matrix (special treatment with row/column corresponding to dirichlet node)   
        for(size_t e=0; e<matrix.size(); ++e)
        {
            const std::vector<int>& elMap = elements[e].indexVectorMask();
            for(size_t i=0; i<matrix[e].size(); ++i)
            {
                size_t iNew= vecToNonDirVecMap[elMap[i]];
                for(size_t j=0; j<matrix[e][i].size();++j)
                {
                    if(!dirMask[elMap[i]] && !dirMask[elMap[j]]) // put coeff if this is not dirichlet row/column
                    {
                        size_t jNew= vecToNonDirVecMap[elMap[j]];
                        triples.push_back( Coeff(iNew, jNew, matrix[e][i][j] ) );
                    }
                    else if(i!=j && dirMask[elMap[j]]) //move column*solution to rhs (if row is not dirichlet)
                    {
                        e_rhs[iNew] -= matrix[e][i][j]*cmpField[elMap[j]];
                    }
                }
            }
        }
        
        //setup matrix from triple list
        e_mat.setFromTriplets(triples.begin(), triples.end());
        e_mat.makeCompressed();
        
    }
    
    template<typename T>
    void EigenSolver<T>::assemblyMatrixVector
    (
        const mesh::Mesh & mesh,
        ElMatrix & matrix,
        ElVector & rhsVector,
        EigenMatrix & e_matrix,
        EigenArray & e_rhsVector   
    ) const 
    {
        using namespace iomanagment;
        
        InfoArrow<<"Assembling local equations into global matrix"<<std::endl;
        
        //Fill vector:
        e_rhsVector.resize(rhsVector.size());
        e_rhsVector=rhsVector;
        
        //Fill triples vector -- vector of sparse matrix coefficients
        size_t tripleSize=0;
        for(size_t i=0; i<matrix.size(); ++i)
            tripleSize+=matrix[i].expSize(); //expSize is equal to the number of elements in array2D
            
        std::vector<Coeff> triples;
        triples.reserve(tripleSize);
        
        for(size_t e=0; e<matrix.size(); ++e)
        {
            const std::vector<int>& mask = mesh[e].indexVectorMask();
            
            for(size_t i=0; i<matrix[e].size(); ++i)
                for(size_t j=0; j<matrix[e][i].size();++j)
                    triples.push_back( Coeff(mask[i], mask[j], matrix[e][i][j] ) );
        }
        
        //setup matrix from triple vector
        e_matrix.setFromTriplets(triples.begin(), triples.end());
        e_matrix.makeCompressed();
    }

    template<typename T>
    void EigenSolver<T>::assemblyMatrixVectorApplyDirichletBC
    (
        const field::GeometricField<T> & solField,
        size_t dim,
        const mesh::Mesh & mesh,
        ElMatrix & matrix,
        ElVector & rhsVector,
        EigenMatrix &e_matrix,
        EigenArray &e_rhsVector
    )const
    {
        using namespace iomanagment;
        
        InfoArrow<<"Assembling local equations into global matrix"<<std::endl;
        
        //Fill vector:
        e_rhsVector=rhsVector;
        
        //Fill triples vector -- vector of sparse matrix coefficients
        size_t tripleSize=0;
        for(size_t i=0; i<matrix.size(); ++i)
            tripleSize+=matrix[i].expSize(); //expSize is equal to the number of elements in array2D
            
        std::vector<Coeff> triples;
        triples.reserve(tripleSize);
        
        //Make dirichlet nodes mask:
        numArray<bool> mask(rhsVector.size(), false);
        numArray<Scalar> dirichletValues(rhsVector.size(),0.);
        
        for(typename field::GeometricField<T>::Patch* patch : solField.boundaryFields() )
        {
            //check if it's exacly dirichlet BC
            if(patch->typeFamily()==field::DIRICHLET)
            {
                const mesh::Boundary & boundary = patch->boundaryEdges();
                
                //iterate over egdes
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const mesh::BoundaryEdge & edge = boundary[e];
                    auto patchValues = CmpTraits<T>::cmpArray((*patch)[e],dim);
                    dirichletValues.slice(edge.gNodesIds()) = patchValues;
                    mask.slice(edge.gNodesIds()) = true;
                }
            }
        }
        
        //Build matrix (special treatment with row/column corresponding to dirichlet node)   
        for(size_t e=0; e<matrix.size(); ++e)
        {
            const std::vector<int>& map = mesh[e].indexVectorMask();
            for(size_t i=0; i<matrix[e].size(); ++i)
            {
                size_t iM = map[i];
                for(size_t j=0; j<matrix[e][i].size();++j)
                {
                    size_t jM = map[j];
                    
                    if(!mask[iM] && !mask[jM]) // put coeff if this is not dirichlet row/column
                        triples.push_back( Coeff(iM, jM, matrix[e][i][j] ) );
                    else if(i!=j && mask[jM]) //move column*solution to rhs (if row is not dirichlet)
                        e_rhsVector[iM] -= matrix[e][i][j]*dirichletValues[jM];
                }
            }
        }
        
        //apply 1. to diagonal and set value to rhs
        for(size_t node=0; node < mask.size(); ++node)
        {
            if(mask[node])
            {
                triples.push_back(Coeff(node,node,1.));
                e_rhsVector[node] = dirichletValues[node];
            }
        }
        
        
        //setup matrix from triple vector
        e_matrix.setFromTriplets(triples.begin(), triples.end());
        e_matrix.makeCompressed();
    }
    
//OLD SOLUTION---> APPLYING 1. INTO DIRICHLET (ROW,COL) IN MATRIX + MOVING COLUMN*DIR_VAL INTO RHS
    template<typename T>
    void EigenSolver<T>::solveEquation(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &rhsVector) const
    {
        using namespace SEM::iomanagment;
        
        //applyDirichletBCInMatrixVector(f,matrix,rhsVector);
        
        size_t dimSize = CmpTraits<T>::dim();
        size_t matrixSize = f.mesh().nodesNumber();
        
        for(size_t dim=0; dim < dimSize; ++dim)
        {
            InfoArrow << "solving field "<<dim<<" component ... " <<std::endl;
            
            EigenMatrix e_mat(matrixSize,matrixSize);
            EigenArray e_rhs(matrixSize);
            EigenArray solution = CmpTraits<T>::cmpArray(f,dim);
            
            // assembly eigen matrix and vector
            //assemblyMatrixVector(f.mesh(), matrix[dim], rhsVector[dim], e_mat, e_rhs);
            assemblyMatrixVectorApplyDirichletBC(f, dim, f.mesh(), matrix[dim], rhsVector[dim], e_mat, e_rhs);
            // solve with eigen solver
            m_linearSolver->solve(e_mat,e_rhs.eigenData(), solution.eigenData());
            
            //assign solution
            CmpTraits<T>::cmpArray(f,dim) = solution;
        }
    }
    
    //NEW APPROACH - REDUCE MATRIX SIZE BY NUMBER OF DIRICHLET NODES
//     template<typename T>
//     void EigenSolver<T>::solveEquation(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &rhsVector) const
//     {
//         using namespace SEM::iomanagment;
//         
//         numArray<bool> mask;
//         numArray<int> vecToNonDirVecMap, nonDirVecToVecMap;
//         
//         // get information how to modify equation system and apply dirichlet boundary conditions
//         // into resulting filed
//         applyDirichletValsAndCreateDirichletInfo(f,mask, vecToNonDirVecMap, nonDirVecToVecMap);
//         
//         
//         size_t dimSize = CmpTraits<T>::dim();
//         size_t matrixSize = nonDirVecToVecMap.size();
//         
//         for(size_t dim=0; dim < dimSize; ++dim)
//         {
//             InfoArrow << "solving field "<<dim<<" component ... " <<std::endl;
//             
//             EigenMatrix e_mat(matrixSize,matrixSize);
//             EigenArray e_rhs(matrixSize);
//             
//             //set initial gues
//             EigenArray solution = CmpTraits<T>::cmpArray(f,dim).slice(nonDirVecToVecMap);
//             
//             // assembly eigen matrix and vector
//             assemblyReducedMatrixVector
//             (
//                 mask,vecToNonDirVecMap,nonDirVecToVecMap,f.mesh(),
//                 CmpTraits<T>::cmpArray(f,dim),matrix[dim],rhsVector[dim],e_mat,e_rhs
//             );
//             // solve with eigen solver
//             m_linearSolver->solve(e_mat,e_rhs.eigenData(), solution.eigenData());
//             
//             //assign solution
//             CmpTraits<T>::cmpArray(f,dim).slice(nonDirVecToVecMap) = solution;
//         }
//     }
    
    
    /** ***********************************************************************
     * \class EigenConjugateGradinet - eigen solver implementation 
     *                                 using conjugate gradinet 
     * ************************************************************************/
    struct EigenConjugateGradient : public EigenLinearSolver
    {
        EigenConjugateGradient(const size_t &maxIter=100, Scalar tolerance=1e-6 );

    protected:
        void solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const;
    private:
        size_t m_maxIter;
        Scalar m_tolerance;
    };

    struct EigenBiCGSTAB : public EigenLinearSolver
    {
            EigenBiCGSTAB(const size_t &maxIter=100, Scalar tolerance=1e-6);
        protected:
            void solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const;
        private:
            size_t m_maxIter;
            Scalar m_tolerance;
    };

    struct EigenSimplicialCholesky : public EigenLinearSolver
    {
    protected:
        void solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const;
    };

    struct EigenSimplicialLLT : public EigenLinearSolver
    {
    protected:
        void solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const;
    };

    struct EigenSimplicialLDLT : public EigenLinearSolver
    {
    protected:
        void solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const;
    };

    
    
    /** ***********************************************************************
     * \class SEM_CG - own ConjugateGradient implementation
     * ************************************************************************/
    template<typename T>
    class SEM_CG : public Solver<T>
    {
        Scalar m_tolerance;
        size_t m_maxIter;
        
        template<typename DerivedArray>
        numArray<Scalar> matrixMul(const mesh::Mesh& mesh, const ElMatrix& matrix, const numArrayBase<DerivedArray>& vec) const
        {
            numArray<Scalar> result(vec.size(),0.);
            for(size_t e=0; e < mesh.size(); ++e)
            {
                const std::vector<int>& elMask = mesh[e].indexVectorMask();
                for(size_t row=0; row < matrix[e].size(); ++row)
                {
                    size_t index = elMask[row];
                    numArrayIndexMapped<const DerivedArray > localVec = vec.slice(elMask);
                    result[index] +=  array::sum( matrix[e][row]* localVec );
                }
            }
            return result;
        }
        
        Scalar arrayDotProd(const numArray<Scalar>& a1, const numArray<Scalar> &a2) const
        {
            return array::sum(a1*a2);
        }
        
        void makeMaskAndApplyDirichletValues(field::GeometricField<T> &solField, std::vector<int> &maskedIndexes) const
        {
            //Make dirichlet nodes mask:
            numArray<bool> mask(solField.size(),false);
            
            size_t dirNodesEstim=0;
            for(typename field::GeometricField<T>::Patch* patch: solField.boundaryFields() )
            {
                //check if it's exacly dirichlet BC
                if(patch->typeFamily()==field::DIRICHLET)
                {
                    const mesh::Boundary & boundary = patch->boundaryEdges();
                    
                    //iterate over egdes
                    for(size_t e=0; e<boundary.size(); ++e)
                    {
                        const mesh::BoundaryEdge & edge = boundary[e];
                        mask.slice(edge.gNodesIds()) = true;
                        solField.slice(edge.gNodesIds()) = (*patch)[e];
                        dirNodesEstim+=edge.gNodesIds().size();
                    }
                }
            }
            
            maskedIndexes.clear();
            maskedIndexes.reserve(dirNodesEstim);
            for(size_t i=0; i< mask.size(); ++i)
            {
                if(mask[i])
                {
                    maskedIndexes.push_back(i);
                }
            }
        }
        
        template<typename Derived>
        void solve(numArrayBase<Derived> &fieldCmp, const mesh::Mesh& mesh, const ElMatrix& matrix, const ElVector& rhsVector, const std::vector<int> &mask) const
        {
            using namespace SEM::iomanagment;
            
            numArray<Scalar> r(rhsVector.size(),0.);
            
            bool normalize = mask.size()==0;
            
            if(normalize)
            {
                fieldCmp -= array::avg(fieldCmp);
                r = rhsVector -array::avg(rhsVector) - matrixMul(mesh, matrix,fieldCmp);            
            }
            else
            {
                r = rhsVector - matrixMul(mesh, matrix,fieldCmp);  
            }
            
            r.slice(mask) = 0.;
            numArray<Scalar> p = r;
            Scalar rsold = arrayDotProd(r,r);
            
            for(size_t iter=0; iter < m_maxIter; ++iter)
            {
                numArray<Scalar> Ap = matrixMul(mesh,matrix,p);
                Scalar alpha = rsold / (arrayDotProd(p,Ap));
                fieldCmp += alpha*p;
                
                if(normalize) 
                    fieldCmp -= array::avg(fieldCmp);
                
                r -= alpha*Ap;
                r.slice(mask) =0.;
                Scalar rsnew = arrayDotProd(r,r);
                
                if(std::sqrt(rsnew) < m_tolerance)
                {
                    InfoArrow << "SEM_CG solution achieved in "<<iter<<" iteration with error:: "<<std::sqrt(rsnew) <<std::endl;
                    return;
                }
                
                p = r + rsnew/rsold * p;
                rsold = rsnew;
            }
            
            InfoArrow << "SEM_CG solution achieved maximum numbero of iteration"<< std::endl;
            
        }
        
    public:
        SEM_CG(Scalar tolerance, size_t maxIter) : m_tolerance(tolerance), m_maxIter(maxIter)
        {
        }
        
        void solveEquation(field::GeometricField<T> &field, SEMMatrix & matrix, SEMVector & rhsVector) const
        {
            std::vector<int> mask;
            
            makeMaskAndApplyDirichletValues(field,mask);
            
            size_t dimSize = CmpTraits<T>::dim();
            for(size_t dim =0; dim < dimSize; ++dim)
            {
                typename CmpTraits<T>::CmpArrayRef fieldCmp = CmpTraits<T>::cmpArray(field,dim);
                solve(fieldCmp, field.mesh(), matrix[dim], rhsVector[dim], mask);
            }
        }
    };
    
    
    template<typename T>
    Solver<T>* getSolver(SolverType solverType)
    {
        EigenLinearSolver * linearSolver;
        switch(solverType)
        {
            case ConjugateGradient:
                linearSolver =  new EigenConjugateGradient(Case::solutionControl().maxIteration(), Case::solutionControl().tolerance());
                break;
            case SimplicialCholesky:
                linearSolver = new EigenSimplicialCholesky();
                break;
            case SimplicialLLT:
                linearSolver= new EigenSimplicialLLT();
                break;
            case SimplicialLDLT:
                linearSolver =  new EigenSimplicialLDLT();
                break;
            case SEM_ConjugateGradient:
                return new SEM_CG<T>(Case::solutionControl().tolerance(), Case::solutionControl().maxIteration());
                
            default:
                linearSolver = new EigenBiCGSTAB(Case::solutionControl().maxIteration(), Case::solutionControl().tolerance());
                break;
        }
        
        return new EigenSolver<T>(linearSolver);   
    }
    
    
    template<typename T>
    Solver<T>* takeChoosenSolver(field::GeometricField<T> &f)
    {
        return getSolver<T>(Case::solutionControl().solverType());
    }
    
    /// \brief applyDirichletBC
    /// Function which removes rows and columns(coressponding to dirichlet-boundary node)
    /// from local matrices, and column*dirichletValue product places in
    /// rhs-vector.In place wher row==column(for dirichlet node) value is set to 1,
    /// to force solution at dirichlet node.
    ///
    /// Important: this method shall be used only when global matrix is build from
    /// local matrices, because only in some local matrices at diagonal value is set
    /// to "1", but in other matrices(which are build from the sam node)
    /// is "0"-->securing duplication of "1" value.
    /// Rhs vector has no such a feature, because rhs vector is treated globaly all the
    /// time, due to fact that rhs-vector in nature is not sparse. Nevertheless rhs-vector
    /// is affected by changes made in this function.
    ///
    /// Important issue is (described widely insied implementation)
    /// -special care when clearing columns/rows. Whole process is subdivided in
    /// 3 loops over patche nodes. Firstly we remove columns/rows and product of column-value,
    /// and we make mask which flags those nodes which are already "dirichletized" on patches
    /// that are setup to be dirichlet.
    /// Due to fact, that boundary is described in manner of edges, so there is possibility
    /// that some of dirichlet nodes appears in non-dirichlet boundary-edges, and local matrix
    /// coressponding to element with that edge can have not been cleaned yet. So we make 2nd
    /// loop over patches, but this time we clear column/row for node only if earlier this
    /// node was masked as dirichlet.
    ///
    /// Finally,we clean mask from flags(we flag this time nodes which already have value 1 setup),
    /// then (when matrix is cleaned from dirichlet node influence) we iterate
    /// again over dirichlet patches, and set value 1 to diagonal of some nodes of some matrices.
    void applyDirichletBC(field::GeometricField<Scalar> & field, ElMatrix &matrix, ElVector &rhsVector);

    /// \brief applyDirichletBCInMatrixVector
    /// Function which removes rows and columns(coressponding to dirichlet-boundary node)
    /// from local matrices, and column*dirichletValue product places in
    /// rhs-vector.In place wher row==column(for dirichlet node) value is set to 1,
    /// to force solution at dirichlet node.
    ///
    /// Important: this method shall be used only when global matrix is build from
    /// local matrices, because only in some local matrices at diagonal value is set
    /// to "1", but in other matrices(which are build from the sam node)
    /// is "0"-->securing duplication of "1" value.
    /// Rhs vector has no such a feature, because rhs vector is treated globaly all the
    /// time, due to fact that rhs-vector in nature is not sparse. Nevertheless rhs-vector
    /// is affected by changes made in this function.
    ///
    /// Important issue is (described widely insied implementation)
    /// -special care when clearing columns/rows. Whole process is subdivided in
    /// 3 loops over patche nodes. Firstly we remove columns/rows and product of column-value,
    /// and we make mask which flags those nodes which are already "dirichletized" on patches
    /// that are setup to be dirichlet.
    /// Due to fact, that boundary is described in manner of edges, so there is possibility
    /// that some of dirichlet nodes appears in non-dirichlet boundary-edges, and local matrix
    /// coressponding to element with that edge can have not been cleaned yet. So we make 2nd
    /// loop over patches, but this time we clear column/row for node only if earlier this
    /// node was masked as dirichlet.
    ///
    /// Finally,we clean mask from flags(we flag this time nodes which already have value 1 setup),
    /// then (when matrix is cleaned from dirichlet node influence) we iterate
    /// again over dirichlet patches, and set value 1 to diagonal of some nodes of some matrices.
    template<typename T>
    void applyDirichletBCInMatrixVector(field::GeometricField<T> &solField, SEMMatrix &matrix, SEMVector &rhsVector)
    {
        using namespace mesh;
        using namespace field;
        
        const Mesh & mesh = solField.mesh();
        
        size_t dimSize = CmpTraits<T>::dim();
        
        numArray<bool> mask(solField.size(),false);
        
        // CLEAR ROWS AND COLUMNS CORRESPONDING TO DIRICHLET BC,
        // MOVE COLUMN *PATCH_VALUE TO RHAS_VECTOR
        // iterate over patches with dirichlet BC
        for(typename GeometricField<T>::Patch* patch: solField.boundaryFields() )
        {
            //check if it's exacly dirichlet BC
            if(patch->typeFamily()==field::DIRICHLET)
            {
                const Boundary & boundary = patch->boundaryEdges();
                
                //iterate over egdes
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const BoundaryEdge & edge = boundary[e];
                    
                    auto patchValues = (*patch)[e];
                    
                    //iterate over nodes in boundary
                    for(const int & matrixNodeIndex : edge.lNodesIdsInElementMatrix())
                    {
                        //iterate over dimensions
                        for(size_t dim=0; dim < dimSize; ++dim)
                        {
                            //get node value
                            auto patchCmpValues =CmpTraits<T>::cmpArray(patchValues,dim);
                            
                            // clear row coresponding to node local index in local matrix
                            matrix[dim][edge.elementId()].row(matrixNodeIndex) = 0.;
                            
                            // Get local matrix column at node local index
                            numArray2D<Scalar>::columnType column = matrix[dim][edge.elementId()].column(matrixNodeIndex);
                            
                            //(note that values for rhsVector[node] won't be substracted, because
                            // we have already set 0 in column[node] while clearing row)
                            rhsVector[dim].slice(edge.element(mesh).indexVectorMask()) -= column*(patchCmpValues);
                            
                            // clear column in local matrix
                            column = 0.;
                        }
                    }
                    
                    //sign nodes which are dirichlet
                    mask.slice(edge.gNodesIds())=true;
                }
            }
        }
        
        // iterate over patches with non dirichlet BC,
        // because in interfacial node between 2 difrent type patches
        // can appear non-cleard* row/column in elemental matrix.
        // BC are defined by edges, and equations are prepared in
        // terms of nodal values. So there is posible situation, when
        // one node is used by 2 edges of 2 difrent type. Further
        // those edges can belong to 2 elements. Finaly, 2 elements
        // have 2 own elemental matrices, but row/column in node
        // place would be cleaned* only in this matrix, where edge was
        // defined as dirichlet BC. To force apropriate cleaning* row/column in
        // elements that "don't know" that have dirichlet BC on some
        // node we do iteration again, but only over those patches
        // which differs from dirichlet patch(those are setup above).
        // Because 1 value in local matrix (node,node) location is set(setup above),
        // so in duplicated nodal entries we now can just simply set
        // apropriate column/row to 0, and substract column values multiplied by
        // dirichlet BC value corresponding to unfortunate node from rhsVector.
        // note- rhsVector for dirichlet node is already corretly setup, so
        // in this place we shall ommit substraction.
        //
        // * cleaned - row/culumn set to 0, diagonal element set to 1,
        //            and rhs vector substracted by column multiplied
        //            by referencing nodal value.
        for (typename GeometricField<T>::Patch* patch : solField.boundaryFields() )
        {
            if(patch->typeFamily()==field::DIRICHLET)
                continue;
            
            const Boundary & boundary = patch->boundaryEdges();
            
            //iterate over egdes
            for(size_t e=0; e<boundary.size(); ++e)
            {
                const BoundaryEdge & edge = boundary[e];
                // get sub vector from dirichlet mask, that would refer to nodes in edge
                numArray<bool>::indexMapped subMask = mask.slice(edge.gNodesIds());
                
                // iterate over nodes in edge
                for(int edgeNodeIndex=0; edgeNodeIndex<subMask.size(); ++edgeNodeIndex)
                {
                    // if node is dirichlet type then in local matrix row and column for this node
                    // must be cleaned
                    if(subMask[edgeNodeIndex])
                    {
                        int matrixNodeIndex=edge.lNodesIdsInElementMatrix()[edgeNodeIndex];
                        
                        //iterate over dimensions in field
                        for(size_t dim=0; dim < dimSize; ++dim)
                        {
                            //clear row
                            matrix[dim][edge.elementId()].row(matrixNodeIndex)=0.;
                            
                            //take column corresponding to node
                            numArray2D<Scalar>::columnType column = matrix[dim][edge.elementId()].column(matrixNodeIndex);
                            
                            // assign column(corresponding to node) multiplied by knonw values(dirichlt
                            // BC value for this node is already stored in rhsVector<--applied in first patch loop).
                            // Note that we first cleaned row, so value in matrix(node,node)
                            // is equal to 0, so applying it to rhsVector will not affect already
                            // correct number.
                            rhsVector[dim].slice(edge.element(mesh).indexVectorMask()) -= column*rhsVector[dim][edge.gNodesIds()[edgeNodeIndex]];
                            
                            //now we can clear values in column(value for matrixNodeIndex again cleaned)
                            column = 0.;
                        }
                    }
                }
            }
        }
        
        
        //clear mask,because it will be again used to flag nodes which already have
        // setup value 1 on diagonal in one local matrix coressponding to dirichlet node.
        // (note, more then 1 local matrix can be linked with the same dirichlet node)
        mask=false;
        
        // SET ON DIAGONAL 1 VALUE (OMMI THE SAME NODE IN OTHER ELMENTS)
        // SETUP DIRICHLET NODES MASK OF PREPARED NODES, APPLY DIRICHLET VALUES TO RHS_VECTOR
        for (typename GeometricField<T>::Patch* patch : solField.boundaryFields() )
        {
            //check if it's exacly dirichlet BC
            if(patch->typeFamily()==field::DIRICHLET)
            {
                const Boundary & boundary = patch->boundaryEdges();
                
                //iterate over egdes
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const BoundaryEdge & edge = boundary[e];
                    auto patchValues = (*patch)[e];
                    
                    //iterate over nodes in boundary
                    for(const int & matrixNodeIndex : edge.lNodesIdsInElementMatrix())
                    {
                        // set coefficient value in diagonal place to 1
                        //-forces exact solution for this node
                        // Place it only when this node haven't been
                        // set as dirichlet yet in mask vector - "1" value
                        // won't be duplicated at asempling process.
                        if(!mask[ edge.element(mesh).indexVectorMask()[matrixNodeIndex] ] )
                        {
                            //iterate over field dimensions
                            for(size_t dim=0; dim < dimSize; ++dim)
                            {
                                matrix[dim][edge.elementId()][matrixNodeIndex][matrixNodeIndex]  = 1.;
                            }
                            mask[ edge.element(mesh).indexVectorMask()[matrixNodeIndex] ]=true;
                        }
                    }
                    // replace elements of rhs vector(known values vector), at indexes of boundary nodes,
                    // to exact solution value
                    //iterate over field dimensions
                    for(size_t dim=0; dim < dimSize; ++dim)
                    {
                        rhsVector[dim].slice(edge.gNodesIds()) = CmpTraits<T>::cmpArray(patchValues,dim);
                    }
                }
            }
        }
        
        // Check if all boundary are neuman type - if yes, then we need to fix solution in one node
        bool allNeuman = true;
        
        for (typename GeometricField<T>::Patch* patch : solField.boundaryFields())
        {
            if(patch->typeFamily()==field::DIRICHLET)
                allNeuman = false;
        }
        
        if(allNeuman)
        {
            typename GeometricField<T>::Patch* patch = *(solField.boundaryFields().begin()+1);
            
            const BoundaryEdge &edge = patch->boundaryEdges()[0];
            
            int lNodeId = edge.lNodesIdsInElementMatrix()[1]; // <-- take 2nd node, what would ensure that it belongs only to only one element
            int gNodeId = edge.gNodesIds()[1];
            int element = edge.elementId();
            
            
            for(size_t dim=0; dim< dimSize; ++dim)
            {
                matrix[dim][element].row(lNodeId) = 0.;
                matrix[dim][element].column(lNodeId) = 0.;
                matrix[dim][element][lNodeId][lNodeId] = 1.;
                rhsVector[dim][gNodeId] = 0.;
            }
        }
    }
    
    
    template<typename T>
    void applyDirichletBCToSolution(field::GeometricField<T> &solField)
    {
        using namespace field;
        
        for(typename GeometricField<T>::Patch* patch: solField.boundaryFields() )
        {
            //check if it's exacly dirichlet BC
            if(patch->typeFamily()==field::DIRICHLET)
            {
                const mesh::Boundary & boundary = patch->boundaryEdges();
                
                //iterate over egdes
                for(size_t e=0; e<boundary.size(); ++e)
                {
                    const mesh::BoundaryEdge & edge = boundary[e];
                    
                    // assign edge values to boundary
                    solField.slice(edge.gNodesIds()) = (*patch)[e];
                }
            }
        }
 
    }
    
    
    template<typename T>
    void solve(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &vector, Solver<T>* solver)
    {
        using namespace SEM::iomanagment;
        
        time::Timer timer;
        solver->solveEquation(f,matrix,vector);
        InfoArrow<<"solution time :"<<timer.elapsed()<<endl;
        EndInfo;
        
        T minValue = SEM::array::min(f);
        T maxValue = SEM::array::max(f);
        
        iomanagment::write(minValue,InfoArrow<<"min")<<endl;
        iomanagment::write(maxValue,InfoArrow<<"max")<<endl;
        
        
        if(!CmpTraits<T>::isValid(minValue) || !CmpTraits<T>::isValid(maxValue) )
            ErrorInFunction<<"solution reusult is not valid"<<endProgram;
        
        delete solver;
    }
    
    template<typename T>
    void solve(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &vector, SolverType solverType)
    {
        solve(f,matrix,vector,getSolver<T>(solverType));
    }
    
    template<typename T>
    void solve(field::GeometricField<T> &f, SEMMatrix &matrix, SEMVector &vector)
    {
        solve(f,matrix,vector,takeChoosenSolver(f));
    }
    
    //Solve explicit equation - asssign build vector into field component fields
    template<typename T>
    void solve(field::GeometricField<T> &f, SEMVector & vector)
    {
           for(size_t d=0; d<CmpTraits<T>::dim(); ++d)
           {
               CmpTraits<T>::cmpArray(f,d) = vector[d];
           }
    }
    
    
    /// \brief solve - function which takes equations builders(objects which can
    /// do discretization of it's part of equation and apply those values to
    /// preallocated matrix and rhs vector), fills matrix and vectors with those
    /// builders, then applys neumen BC (clear columns and rows from matrix, and
    /// adds known values to rhsVector), finaly chose solver which is specified
    /// to use for this equation and performs computation with assigning result
    /// to field.
    template<typename T, typename Implicit, typename Explicite>
    void solve(DiscretEquation<T,Implicit, Explicite> equation, Solver<T>* solver)
    {
        using namespace iomanagment;

        
        SEMVector rhsVector;
        SEMMatrix matrix;
        
        time::Timer timer;
        Info<<"Building equations ..."<<std::endl;
        
        if( equation.isDiagonal() )
        {
            Info<<"Selected explicit solver ..."<<std::endl;
            applyDirichletBCToSolution(equation.field());
            equation.buildExplicitEquation(rhsVector);
        }
        else
        {
            equation.buildEquation(matrix,rhsVector);
        }
        
        InfoArrow<<"Equation system build in: "<<timer.elapsed()<<std::endl;

        if( equation.isDiagonal() )
        {
            solve(equation.field(), rhsVector);
        }
        else
        {
            solve(equation.field(), matrix, rhsVector, solver);
            InfoArrow<<"Equation system solved in: "<<timer.elapsed()<<std::endl;
        }
        
    }
    
    template<typename T, typename Implicit, typename Explicite>
    void solve(DiscretEquation<T,Implicit,Explicite> equation)
    {
        solve(equation, takeChoosenSolver( equation.field() ) );
    }
    
    template<typename T, typename Implicit, typename Explicite>
    void solve(DiscretEquation<T,Implicit,Explicite> equation, SolverType solverType)
    {
        solve(equation,getSolver<T>(solverType) );
    }

    template<typename T, typename Implicit>
    void solve(boost::shared_ptr<DiscretOperator<T,Implicit> > implOperator)
    {
        typedef DummyExplicitOperator<T> Dummy;
        typename Dummy::ref explOperator(new Dummy);
        DiscretEquation<T,Implicit,Dummy> equation(implOperator,explOperator);
        solve(equation, takeChoosenSolver(equation.field()) );
    }
    
    template<typename T, typename Implicit>
    void solve(boost::shared_ptr<DiscretOperator<T,Implicit> > implOperator, SolverType solverType)
    {
        typedef DummyExplicitOperator<T> Dummy;
        typename Dummy::ref explOperator(new Dummy);
        DiscretEquation<T,Implicit,Dummy> equation(implOperator,explOperator);
        solve( equation, getSolver<T>(solverType) );
    }
    

    //-------------------------------------------------------------
    
    template<typename T, typename Implicit, typename Explicite>
    void makeEquationSystem(DiscretEquation<T,Implicit, Explicite> equation, SEMMatrix & matrix, SEMVector & vector)
    {
        time::Timer timer;
        Info<<"Building equations ..."<<std::endl;
        equation.buildEquation(matrix,vector);
        InfoArrow<<"Equation system build in: "<<timer.elapsed()<<std::endl;
    }
    
    template<typename T, typename Implicit, typename Explicite>
    void makeEquationSystem(boost::shared_ptr<DiscretEquation<T,Implicit, Explicite> > equation, SEMMatrix & matrix, SEMVector & vector)
    {
        makeEquationSystem(*equation, matrix, vector);
    }
    
    template<typename T, typename Implicit>
    void makeEquationSystem(boost::shared_ptr<DiscretOperator<T,Implicit> > implOperator, SEMMatrix & matrix, SEMVector & vector)
    {
        typedef DummyExplicitOperator<T> Dummy;
        typename Dummy::ref explOperator(new Dummy);
        DiscretEquation<T,Implicit,Dummy> equation(implOperator,explOperator);
        makeEquationSystem(equation, matrix, vector);
    }
    
}//las
}//SEM








#endif //_Solver_H_
