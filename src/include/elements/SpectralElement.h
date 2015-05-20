#ifndef _SpectralElement_H
#define _SpectralElement_H

#include <list>
#include <vector>
#include <iostream>

#include "boost/array.hpp"

#include "utilities/TypeDefs.h"
#include "utilities/Reference.h"
#include "utilities/numArray2D.h"
#include "components/Basic.h"
#include "utilities/ImplementationFactory.h"

namespace SEM{ namespace mesh{

    typedef Matrix<boost::array<double,6> >::type HMatrix;

    class  SpectralElement
    {
        REFERENCE_TYPE(SpectralElement)

        DECLARE_IMPLEMENTATION_FACTORY(SpectralElement,unsigned int)
        
    protected:
        unsigned int NSpec;
        unsigned int numberOfInteriorNodes;
        std::vector<int> numberOfInteriorInEdgesNodes;

        std::string type;

        std::vector<double> nodes;
        std::vector<double> weights;
        Matrix<double>::type D_Matrix;//derivative matrix

        VectorToMatrixMap m_localIndexesInMatrix;
        
        numArray<Scalar> m_barycentricInterpWeights;
        

    public:
        SpectralElement();

        SpectralElement( unsigned int NSpec, const std::string &typeName);

        int getNumberOfInteriorNodes();

        int getNumberOfInteriorInEdgeNodes(int edgeNumber);

        int getSpectralSize() const;

        std::vector<double> getSpectralNodes();

        std::vector<double> getSpectralWeights();

        virtual unsigned int numberOfRealNodes() const=0;
        
        /// Return element node indexes in local matrix row/column
        inline const VectorToMatrixMap & localNodesInMatrix() const
        {
            return m_localIndexesInMatrix;
        }

        /// Return matrix with coordinates in real space. Coordinates need to be placed in local
        /// numeration of indexes (they are translated later by mask)
        virtual void generateSpectralNodesRealCoordinates
        (
            Matrix<Vector>::type & rSCoords,
            const std::vector<Vector> & rCoords
        )=0;

        /// THIS METHOD SHOULD FILL UP NODES,WEIGHTS FIELDS, ACORDINGLY TO THIS ELEMENT TYPE DEFINITION
        virtual void generateSpectralNodesAndWeights ()=0;

        /// THIS METHOD SHOULD FILL UP DERIVATVE MATIRX
        /// (derivative marix-> d(l(ksi_i))/d(ksi_j), where l- shape function, ksi_i - spectral node at i index))
        virtual void generateDerivativeMatrix()=0;

        
        virtual Vector mapToLocal(const Vector &point, const std::vector<Vector> & nodes ) const =0;
        
        
        /// function to calculate interpolation coefficients. Dot product of returned
        /// vector and nodal solution in this element would produce value of interpolation
        /// local ksi,etha coordinates.
        virtual numArray<Scalar> computeInterpolationCoefficients(Scalar ksi, Scalar etha) const =0;

        /// \brief getEdgeLocalNodesIdInMatrix
        ///         Instead of getEdgeLocalNodesId it returns edge ids in vector.
        ///         Eeach local ID in this vector referes to matrix equation Row/Colum
        ///         in local numeration. "getEdgeLocalNodesId" evaluates mapping in the
        ///         same manner as nodes are geometricaly stored, and this method in manner
        ///         as they are represented in local equation matrix.
        /// \param edgeId -one of edges, numerated from 0-bottom counterclockwise
        /// \return vector with mapping local numeration of all nodes.
        virtual std::vector<int> getEdgeLocalNodesIdInMatrix(int edgeId) const =0;

        /// THIS METHOD NEED TO RETURN ELEMENTS FROM "MASK" THAT CORESPONDS TO "EDGE_ID'S" EDGE
        virtual std::vector<int> getEdgeGloblaNodesId(int edgeId, const Matrix<int>::type & mask) const =0;

        /// THIS METHOD NEED TO RETURN ELEMENTS THAT LOCALY CORESPONDS TO "EDGE_ID'S" EDGE
        virtual std::vector<boost::array<int,2> > getEdgeLocalNodesId(int edgeId) const =0;

        /// DEFINE THE SIZE OF ELEMENT MASK FPR ALL SPECTRAL NODES ASIGNED TO THIS (THE ARRAY)
        /// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
        virtual void defineMaskArraySize(Matrix<int>::type & mask)=0;

        /// ALLOCATE IN PROPER PLACE IN MASK ARRAY THE ELEMENT VERTEX NODES GIVEN GLOBAL ADRESS(GIVEN FROM MESH, NOT GENERATED AS SPECTRAL)
        /// 0-always left bottom, next vertex counterclockwise
        /// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
        virtual void allocateMaskArray_VertexNode
            (
            int vertexNr,
            int nodeGlobalAdress,
            Matrix<int>::type & addressArray
            )=0;

        /// ALLOCATE INTO GLOBAL ADRESS MASK ARRAY SPECTRAL NODES ON EDGE(WITHOUT VERTEX NODE)
        /// 0-bottom edge, and next counterclockwise. Nodes are given in this direction, means alwasy next node is this, which is counterclockwise
        /// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
        virtual void allocateMaskArray_NodesInEdge
            (
            int edgeNr,
            std::vector<int> edgeGlobalAdress,
            Matrix<int>::type & addressArray
            )=0;

        /// ALLOCATE INTO GLOBAL ADRESS MASK ARRAY SPECTRAL NODES INSIDE NODES
        /// NODES SHOULD BE INDEXED FROM START ADDRESS TO SOME APROPRIATE VALUE
        /// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
        virtual void allocateMaskArray_InteriorNodes
            (
            int startAddress,
            Matrix<int>::type & addressArray
            )=0;

        /// \brief transformMatrixMaskToVectorMask
        /// \param mask - local matrix mask of global indexes
        /// \param vecMask(return value) - local vectoral mask of global indexes,
        ///        just like columns/rows are assembled from elemental matrix to global matrix
        virtual void transformMatrixMaskToVectorMask(const Matrix<int>::type &mask, std::vector<int> & vecMask) const =0;

        /// Fill up the "H-Matrix"---> this is the evaluation of this parts of equation in laplacian discret week form:
        /// H_Matrix[i][j][0] = (ksix*ksix)*J*W[i]*W[j];          H_Matrix[i][j][3] = (ksiy*ksiy)*J*weights[i]*weights[j];
        /// H_Matrix[i][j][1] = (ksix*etax)*J*W[i]*W[j];			 H_Matrix[i][j][4] = (ksiy*etay)*J*weights[i]*weights[j];
        /// H_Matrix[i][j][2] = (etax*etax)*J*W[i]*W[j];			 H_Matrix[i][j][5] = (etay*etay)*J*weights[i]*weights[j];
        /// i,j - element spectral node index
        /// W - weights
        /// J - jacobian at i,j
        virtual void generateH_Matrix
            (
            HMatrix & H_Matrix,
            const std::vector< Vector > & realCoords
            )=0;

        /// Fill up the "M-Matrix"---> this is part of equation required for building the mass matrix. Its for is:
        /// M_Matrix[i][j] = J*W.V[i]*W.V[j];
        /// i,j - element spectral node index
        /// W - weights
        /// J - jacobian at i,j
        virtual void generateM_Matrix
            (
            numArray2D<double> & M_Matrix,
            const std::vector< Vector > & realCoords
            )const=0;

        /// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
        /// with scalar value-constant for all nodes
        /// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
        ///  i,j - coresponds to local i,j node)
        /// Returns matrix in local nodal form
        virtual void generateStiff_Matrix
            (
            Matrix4<double>::type & Stiff_Matrix,
            double scalarValue,
            const HMatrix & H_Matrix
            )const=0;

        /// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
        /// with scalar value-constant for all nodes
        /// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
        ///  i,j - coresponds to local i,j node)
        /// Returns matrix in elemental matrix form
            virtual numArray2D<Scalar> generateStiff_Matrix
            (
                const double& scalarValue,
                const HMatrix & H_Matrix
            ) const=0;

            
        /// \function normal
        /// function which shall calculate edge normal vector
        virtual Vector normal(unsigned int edgeId, const std::vector<Vector> &rC) const=0;
            
       
        /// \brief generateWeekVddx 
        /// function which evaluates expression v[s] = integrate( u * d/dx(phi_qq) )
        /// where p,q indexes coressponding to local nodal values, and s-it's 
        /// representation in element vector.
       virtual void generateWeekVddx
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const =0;
            
        /// \brief generateWeekVddy 
        /// function which evaluates expression v[s] = integrate( u * d/dx(phi_pq) )
        /// where p,q indexes coressponding to local nodal values, and s-it's 
        /// representation in element vector.            
       virtual void generateWeekVddy     
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const =0;
       
        /// \brief generate_ddx 
        /// function which evaluates derivative d(u)/dx in spectral nodes
        virtual void generate_ddx
            (
                const numArray<Scalar>& u, 
                const std::vector<Vector> &rCoords, 
                numArray<Scalar>& result
            ) const =0;
        
        /// \brief generate_ddy 
        /// function which evaluates derivative d(u)/dx in spectral nodes
        virtual void generate_ddy
        (
            const numArray<Scalar>& u, 
            const std::vector<Vector> &rCoords, 
            numArray<Scalar>& result
        ) const =0;
       
        virtual numArray<Scalar> rot(const numArray<Vector> & U, const std::vector<Vector> &rCoords ) const=0;
        
        
        //Calculate boundary integral of N x Rot(Rot(U))*Phi, where N-normal, U-vector field, Phi-test function
       virtual numArray<Scalar> boundInt_NxRotRotU(const numArray<Vector> & U, const std::vector<Vector> &rCoords, const size_t &edge) const=0;
            

       virtual numArray2D<Scalar> convMatrix(const numArray<Vector> &U, const std::vector<Vector> &rC) const =0; 
       
        ///// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
        ///// with tensor coefficient value-constant for all nodes
        ///// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
        /////  i,j - coresponds to local i,j node)
        //virtual void generateStiff_Matrix
        //(
        //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
        //	array<double,4>& tensorValue,
        //	vector<vector<array<double,3>>>& H_Matrix
        //)=0;


        ///// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
        ///// with scalar coefficient field
        ///// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
        /////  i,j - coresponds to local i,j node)
        //virtual void generateStiff_Matrix
        //(
        //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
        //	vector<vector<double>>& scalarNodalValues,
        //	vector<vector<array<double,3>>>& H_Matrix
        //)=0;

        ///// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
        ///// With tensor coefficient filed
        ///// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
        /////  i,j - coresponds to local i,j node)
        //virtual void generateStiff_Matrix
        //(
        //  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
        //	vector<vector<array<double,4>>>& tensorNodalValues,
        //	vector<vector<array<double,3>>>& H_Matrix
        //)=0;


        /// Calcualte coefficient value of neuman BC integral for nodes in apropriate edge
        /// This value is of the form: int_edge(h*Phipq)dEdge = h*sqrt((dx/dksi)^2+(dy/dksi)^2)*wp*w0
        /// This value multiplied by flux value on edge gives influence of Neuman BC to equation(known value)
        virtual void generateEdgeNeumanBCIntegralCoefficient
        (
            int edgeNumber,
            const std::vector<Vector> & elementCoordinates,
            std::vector<double> & resultNodalVal
        )=0;


        virtual ~SpectralElement()
        {
        }

    protected:
        static numArray<Scalar> computeInterpolationBarycentricWeights(const std::vector<Scalar> & nodes);
        
    protected:
        
        /// RETRUN THE NUMBER OF INTERNIOR NODES ( THE COUNT OF NODES, WHERE NONE LIES ON ANY EDGE OR VERTEX )
        virtual int generateNumberOfInteriorNodes()=0;

        virtual void generateLocalIndexesInMatrix(VectorToMatrixMap &map) const =0;

        /// RETRUN NUMBERS(ON EACH EDGE RESPECTIVLY) OF INTERNIOR IN EDGE NODES ( THE COUNT OF NODES, WHERE NONE LIES VERTEX OR INSIDE )
        virtual std::vector<int> generateNumberOfInteriorInEdgesNodes()=0;


        friend std::ostream& operator <<( std::ostream& o, const SpectralElement& se);
        
        
    };


}//mesh
}//SEM
#endif //_SpectralElement_H
