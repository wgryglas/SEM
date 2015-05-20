#ifndef _SE_TLGL_H
#define _SE_TLGL_H

#include <vector>

#include "boost/array.hpp"

#include "SpectralElement.h"
#include "components/Basic.h"
#include "utilities/Reference.h"

namespace SEM{ namespace mesh
{
        /** **********************************************
         *  \class SE_TLGL
         *  Implementation of SpectralElement interface
         *  which allows to discretize triangular elements
         *************************************************/
	class SE_TLGL: public SpectralElement
	{
        REFERENCE_TYPE(SE_TLGL)
		public:
            SE_TLGL();

            SE_TLGL(unsigned int NSpec);
            
            unsigned int numberOfRealNodes() const;
            

            void generateSpectralNodesAndWeights();

            void generateSpectralNodesRealCoordinates
            (
                Matrix<Vector>::type & rSCoords,
                const std::vector<Vector> & rCoords
            );

            void generateDerivativeMatrix();
            
            Vector mapToLocal(const Vector &point, const std::vector<Vector> & nodes ) const;
            
            numArray<Scalar> computeInterpolationCoefficients(Scalar ksi, Scalar etha) const;

            /// \brief getEdgeLocalNodesIdInMatrix
            ///         Instead of getEdgeLocalNodesId it returns edge ids in vector.
            ///         Eeach local ID in this vector referes to matrix equation Row/Colum
            ///         in local numeration. "getEdgeLocalNodesId" evaluates mapping in the
            ///         same manner as nodes are geometricaly stored, and this method in manner
            ///         as they are represented in local equation matrix.
            /// \param edgeId -one of edges, numerated from 0-bottom counterclockwise
            /// \return vector with mapping local numeration of all nodes.
            std::vector<int> getEdgeLocalNodesIdInMatrix(int edgeId) const;

            std::vector<int> getEdgeGloblaNodesId(int edgeId, const Matrix<int>::type & mask) const;

            std::vector<boost::array<int,2> > getEdgeLocalNodesId(int edgeId) const;

            void defineMaskArraySize(Matrix<int>::type & mask);

            void allocateMaskArray_VertexNode
            (
                int vertexNr, 
                int nodeGlobalAdress,
                Matrix<int>::type & addressArray
            );

            void allocateMaskArray_NodesInEdge
            (
                int edgeNr,
                std::vector<int> edgeGlobalAdress,
                Matrix<int>::type & addressArray
            );

            void allocateMaskArray_InteriorNodes
            (
                int startAddress,
                Matrix<int>::type & addressArray
            );

            /// \brief transformMatrixMaskToVectorMask
            /// \param mask - local matrix mask of global indexes
            /// \param vecMask(return value) - local vectoral mask of global indexes,
            ///        just like columns/rows are assembled from elemental matrix to global matrix
            void transformMatrixMaskToVectorMask(const Matrix<int>::type &mask, std::vector<int> & vecMask) const;

            Scalar x_ksi(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar y_ksi(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar x_etha(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar y_etha(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar jacobianCoeff(const std::vector<Vector>& rC) const;
            
            Scalar jacobian(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar ksi_x(int i, int j,  const std::vector<Vector> &rC) const;
            
            Scalar ksi_y(int i, int j,  const std::vector<Vector> &rC) const;
            
            Scalar etha_x(int i, int j, const std::vector<Vector> &rC) const;
            
            Scalar etha_y(int i, int j, const std::vector<Vector> &rC) const;
            
            Vector normal(unsigned int edgeId, const std::vector<Vector> &rC) const;
            
            void generateH_Matrix
                (
                    HMatrix & H_Matrix,
                    const std::vector<Vector> & realCoords
                );

            void generateM_Matrix
                (
                numArray2D<double> & M_Matrix,
                const std::vector<Vector> & realCoords
                ) const;

            void generateStiff_Matrix
                (
                Matrix4<double>::type & Stiff_Matrix,
                double scalarValue,
                const HMatrix & H_Matrix
                ) const;

            numArray2D<Scalar> generateStiff_Matrix
            (
                const double& scalarValue,
                const HMatrix & H_Matrix
            ) const;

            //Calculate rotation of vector filed in node ksi,etha
            Scalar rotU( const numArray<Vector> &u, size_t ksi, size_t etha, const std::vector<Vector> &rC) const;
            
            numArray<Scalar> rot(const numArray<Vector> & U, const std::vector<Vector> &rCoords ) const;
            
            //Calculate transofromation jacobian of 1D integral over specified edge
            Scalar gamma(size_t node, unsigned int edge, const std::vector<Vector> & rC) const;
            
            void generateWeekVddx
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const;
            
            void generateWeekVddy     
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const;
            
            void generate_ddx
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const;
            
            void generate_ddy
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const;
            
            numArray<Scalar> boundInt_NxRotRotU
            (
                const numArray<Vector> & U, 
                const std::vector<Vector> &rCoords, 
                const size_t &edge
            ) const;
            
            numArray2D< Scalar > convMatrix
            (
                const numArray< Vector >& U, 
                const std::vector< Vector >& rC
            ) const;
            
            // void generateStiff_Matrix
            //(
            //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
            //	array<double,4>& tensorValue,
            //	vector<vector<array<double,6>>>& H_Matrix
            //)=0;


            //void generateStiff_Matrix
            //(
            //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
            //	vector<vector<double>>& scalarNodalValues,
            //	vector<vector<array<double,6>>>& H_Matrix
            //);

            //void generateStiff_Matrix
            //(
            //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
            //	vector<vector<array<double,4>>>& tensorNodalValues,
            //	vector<vector<array<double,6>>>& H_Matrix
            //);
            
            void generateEdgeNeumanBCIntegralCoefficient
            (
                int edgeNumber,
                const std::vector<Vector> & elementCoordinates,
                std::vector<double>& resultNodalVal
            );


        protected:
            int generateNumberOfInteriorNodes();

            void generateLocalIndexesInMatrix(VectorToMatrixMap &map) const;

            std::vector<int> generateNumberOfInteriorInEdgesNodes();

            virtual ~SE_TLGL();

        private:

            void qdqLN(double x, int N, double*q, double*dq, double*LN);

            double WJ(int j);
        };

        std::ostream& operator<<(std::ostream& o, SE_TLGL& e);
} //mesh
} //SEM

#endif //_SE_TLGL_H
