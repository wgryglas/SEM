
#include <iostream>

#include"Solver.h"

#include "time/Timmer.h"

namespace SEM { namespace las {


EigenConjugateGradient::EigenConjugateGradient(const size_t &maxIter, Scalar tolerance)
    : m_maxIter(maxIter), m_tolerance(tolerance)
{
}

void EigenConjugateGradient::solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const
{
    using namespace iomanagment;
    using namespace std;

    InfoArrow<<"Selceted ConjugateGradient linear solver..."<<endl;

    Eigen::ConjugateGradient<EigenMatrix> solver(matrix);
    solver.setMaxIterations(m_maxIter);
    solver.setTolerance(m_tolerance);

    solution = solver.solveWithGuess(rhsVector, solution);

    InfoArrow<<"solver iteration: #"<<solver.iterations()<<endl;
    InfoArrow<<"solutution estimated error: "<<solver.error()<<endl;
}

EigenBiCGSTAB::EigenBiCGSTAB(const size_t &maxIter, Scalar tolerance)
    : m_maxIter(maxIter), m_tolerance(tolerance)
{
}

void EigenBiCGSTAB::solve(EigenMatrix &matrix, EigenVector &rhsVector,EigenVector &solution) const
{
    using namespace iomanagment;
    using namespace std;

    InfoArrow<<"Selceted BiCGSTAB linear solver..."<<endl;

    Eigen::BiCGSTAB<EigenMatrix, Eigen::IncompleteLUT<Scalar> > solver(matrix);
    solver.setMaxIterations(m_maxIter);
    solver.setTolerance(m_tolerance);

    solution = solver.solveWithGuess(rhsVector, solution);

    InfoArrow<<"solver iteration: #"<<solver.iterations()<<endl;
    InfoArrow<<"solutution estimated error: "<<solver.error()<<endl;
}

void EigenSimplicialCholesky::solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const
{
    using namespace iomanagment;
    using namespace std;

    InfoArrow<<"Selceted SimplicialCholesky linear solver..."<<endl;

    Eigen::SimplicialCholesky<EigenMatrix> solver(matrix);
    solution = solver.solve(rhsVector);
}

void EigenSimplicialLLT::solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const
{
    using namespace iomanagment;
    using namespace std;

    InfoArrow<<"Selceted SimplicialLLT linear solver..."<<endl;

    Eigen::SimplicialLLT<EigenMatrix> solver(matrix);
    solution = solver.solve(rhsVector);

}

void EigenSimplicialLDLT::solve(EigenMatrix &matrix,EigenVector &rhsVector,EigenVector &solution) const
{
    using namespace iomanagment;
    using namespace std;

    InfoArrow<<"Selceted SimplicialLDLT linear solver..."<<endl;
    Eigen::SimplicialLDLT<EigenMatrix> solver(matrix);
    solution = solver.solve(rhsVector);
}

void applyDirichletBC(field::GeometricField<Scalar> &field, ElMatrix &matrix, ElVector &rhsVector)
{
    using namespace mesh;
    using namespace field;

    const Mesh & mesh = field.mesh();

    numArray<bool> mask(field.size(),false);

    // CLEAR ROWS AND COLUMNS CORRESPONDING TO DIRICHLET BC,
    // MOVE COLUMN *PATCH_VALUE TO RHAS_VECTOR
    // iterate over patches with dirichlet BC
    for(typename GeometricField<Scalar>::Patch* patch : field.boundaryFields() )
    {
        //check if it's exacly dirichlet BC
        if(patch->typeFamily()==field::DIRICHLET)
        {
            const Boundary & boundary = patch->boundaryEdges();

            //iterate over egdes
            for(size_t e=0; e<patch->size(); ++e)
            {
                const BoundaryEdge & edge = boundary[e];
                auto patchValues = (*patch)[e];

                //iterate over nodes in boundary
                for(const int & matrixNodeIndex: edge.lNodesIdsInElementMatrix())
                {
                    // clear row coresponding to node local index in local matrix
                    matrix[edge.elementId()].row(matrixNodeIndex) = 0.;

                    // Get local matrix column at node local index
                    numArray2D<Scalar>::columnType column = matrix[edge.elementId()].column(matrixNodeIndex);

                    //(note that values for rhsVector[node] won't be substracted, because
                    // we have already set 0 in column[node] while clearing row)
                    rhsVector.slice(edge.element(mesh).indexVectorMask()) -= column*(patchValues);

                    // clear column in local matrix
                    column = 0.;
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
    for (typename GeometricField<Scalar>::Patch* patch : field.boundaryFields() )
    {
        if(patch->typeFamily()==field::DIRICHLET)
            continue;

        const Boundary & boundary = patch->boundaryEdges();

        //iterate over egdes
        for(size_t e=0; e<patch->size(); ++e)
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

                    //clear row
                    matrix[edge.elementId()].row(matrixNodeIndex)=0.;

                    //take column corresponding to node
                    numArray2D<Scalar>::columnType column = matrix[edge.elementId()].column(matrixNodeIndex);

                    // assign column(corresponding to node) multiplied by knonw values(dirichlt
                    // BC value for this node is already stored in rhsVector<--applied in first patch loop).
                    // Note that we first cleaned row, so value in matrix(node,node)
                    // is equal to 0, so applying it to rhsVector will not affect already
                    // correct number.
                    rhsVector.slice(edge.element(mesh).indexVectorMask()) -= column*rhsVector[edge.gNodesIds()[edgeNodeIndex]];
                    
                    //now we can clear values in column(value for matrixNodeIndex again cleaned)
                    column = 0.;
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
    for (typename GeometricField<Scalar>::Patch* patch : field.boundaryFields() )
    {
        //check if it's exacly dirichlet BC
        if(patch->typeFamily()==field::DIRICHLET)
        {
            const Boundary & boundary = patch->boundaryEdges();

            //iterate over egdes
            for(size_t e=0; e<patch->size(); ++e)
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
                        matrix[edge.elementId()][matrixNodeIndex][matrixNodeIndex]  = 1.;
                        mask[ edge.element(mesh).indexVectorMask()[matrixNodeIndex] ]=true;
                    }
                }
                // replace elements of rhs vector(known values vector), at indexes of boundary nodes,
                // to exact solution value
                rhsVector.slice(edge.gNodesIds()) = patchValues;
            }
        }
    }


}


}//las
}//SEM
