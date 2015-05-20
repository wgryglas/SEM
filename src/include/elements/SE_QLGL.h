#ifndef _SE_QLGL_H
#define _SE_QLGL_H

#include <iostream>
#include <vector>
#include <list>

#include "boost/array.hpp"

#include "Eigen/Dense"

#include "SpectralElement.h"
#include "components/Basic.h"
#include "utilities/Reference.h"

namespace SEM{ namespace mesh{

    class SE_QLGL: public SpectralElement
    {
        REFERENCE_TYPE(SE_QLGL)
    public:
        SE_QLGL();

        SE_QLGL(unsigned int NSpec);

        unsigned int numberOfRealNodes() const;
        
        void generateSpectralNodesAndWeights();

        void generateSpectralNodesRealCoordinates
            (
            Matrix<Vector>::type & rSCoords,
            const std::vector<Vector> & rCoords
            );

        void generateDerivativeMatrix();
        
        numArray<Scalar> computeInterpolationCoefficients(Scalar ksi, Scalar etha) const;
        
        Vector mapToLocal(const Vector &point, const std::vector<Vector> & nodes ) const;
        
        double jacobian(int ksi, int etha, const std::vector<Vector>& rC) const;
        double jacobian(double xksi, double xetha, double yksi, double yetha) const;
        double x_ksi(int ksi, int etha, const std::vector<Vector>& rC) const;
        double x_etha(int ksi, int etha, const std::vector<Vector>& rC) const;
        double y_ksi(int ksi, int etha, const std::vector<Vector>& rC) const;
        double y_etha(int ksi, int etha, const std::vector<Vector>& rC) const;
        
        double ksi_x(int ksi, int etha, const std::vector<Vector>& rC, const Scalar &jacobian) const
        {
            return y_etha(ksi,etha,rC)/jacobian;
        }
        double ksi_y(int ksi, int etha, const std::vector<Vector>& rC, const Scalar &jacobian) const
        {
            return -x_etha(ksi,etha,rC)/jacobian;
        }
        double etha_x(int ksi, int etha, const std::vector<Vector>& rC, const Scalar &jacobian) const
        {
            return -y_ksi(ksi,etha,rC)/jacobian;
        }
        double etha_y(int ksi, int etha, const std::vector<Vector>& rC, const Scalar &jacobian) const
        {
            return x_ksi(ksi,etha,rC)/jacobian;
        }
        
        Vector normal(unsigned int edge, const std::vector<Vector> &rC) const;
        
    
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

        std::vector< boost::array<int,2> > getEdgeLocalNodesId(int edgeId) const;

        void defineMaskArraySize(Matrix<int>::type & mask);

        void allocateMaskArray_VertexNode(int vertexNr, int nodeGlobalAdress, Matrix<int>::type & addressArray);

        void allocateMaskArray_NodesInEdge
            (
            int edgeNr,
            std::vector<int> edgeGlobalAdress,
            Matrix<int>::type & addressArray
            );

        /// \brief transformMatrixMaskToVectorMask
        /// \param mask - local matrix mask of global indexes
        /// \param vecMask(return value) - local vectoral mask of global indexes,
        ///        just like columns/rows are assembled from elemental matrix to global matrix
        void transformMatrixMaskToVectorMask(const Matrix<int>::type &mask, std::vector<int> & vecMask) const;

        void allocateMaskArray_InteriorNodes(int startAddress, Matrix<int>::type & addressArray);

        void generateH_Matrix
            (
            HMatrix & H_Matrix,
            const std::vector< Vector > & realCoords
            );

        void generateM_Matrix
            (
            numArray2D<Scalar> & M_Matrix,
            const std::vector< Vector > & realCoords
            )const;

        void generateStiff_Matrix
            (
            Matrix4<double>::type & Stiff_Matrix,
            double scalarValue,
            const HMatrix& H_Matrix
            ) const;

        numArray2D<Scalar> generateStiff_Matrix
            (
                const double &scalarValue,
                const HMatrix & H_Matrix
            ) const;

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
            
            
        numArray<Scalar> rot(const numArray<Vector> & U, const std::vector<Vector> &rCoords ) const
        {
            numArray<Scalar> result(NSpec*NSpec,0.);
            for(size_t ksi=0; ksi<NSpec; ++ksi)
            {
                for(size_t etha=0; etha<NSpec; ++etha)
                {
                    result[ksi+etha*NSpec] = rotU(U,ksi,etha,rCoords);
                }
            }
                
            return result;
        }
            
        numArray2D< Scalar > convMatrix
        (
            const numArray< Vector >& U,
            const std::vector< Vector >& rC
        ) const;            
            
    private:            
        /// \function rotU 
        /// calculates rotation of vector filed in node ksi,etha
        Scalar rotU(const numArray< Vector >& U,size_t ksi, size_t etha, const std::vector<Vector> & rCoords) const;
        
        /// \function gamma 
        /// calculates coefficient required for apropiate integral caluclacion (jacobian of 1D integral transofrmation)
        /// it's value is equal to dS/dt(S-edge length, dt-independent variable on edge).
        /// In case of linear mapping value is equal to edge length divided by 2
        Scalar gamma(size_t node, unsigned int edge, const std::vector<Vector> &rCoords) const;
        
    public:
        numArray<Scalar> boundInt_NxRotRotU(const numArray<Vector> & U, const std::vector<Vector> &rCoords, const size_t &edge) const;
            
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
            const std::vector< Vector > & elementCoordinates,
            std::vector<double> & resultNodalVal
        );


    protected:
        int generateNumberOfInteriorNodes();

        void generateLocalIndexesInMatrix(VectorToMatrixMap &map) const;

        std::vector<int> generateNumberOfInteriorInEdgesNodes();

        virtual ~SE_QLGL();

    private:

        void qdqLN(double x, int N, double*q, double*dq, double*LN);

        double WJ(int j);

    };


    std::ostream& operator<<(std::ostream& o, SE_QLGL& e);

}//mesh
}//SEM

#endif //_SE_QLGL_H

