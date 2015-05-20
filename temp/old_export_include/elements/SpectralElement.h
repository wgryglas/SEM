#ifndef _SpectralElement_H
#define _SpectralElement_H

#include <list>
#include <vector>
#include <array>
#include <iostream>

using namespace std;

namespace SEM
{
	namespace MESH
	{

		class  SpectralElement
		{

		protected:
			int NSpec;
			int numberOfInteriorNodes;
			vector<int> numberOfInteriorInEdgesNodes;

			std::string type;

			vector<double> nodes;
			vector<double> weights;
			vector<vector<double> > D_Matrix;//derivative matrix
		public:
			SpectralElement
				(
				);

			SpectralElement
				(
				const int NSpec,
				const std::string type
				);

			int getNumberOfInteriorNodes
				(
				);

			int getNumberOfInteriorInEdgeNodes
				(
				int edgeNumber
				);


			int getSpectralSize();

			vector<double> getSpectralNodes
				(
				);

			vector<double> getSpectralWeights
				(
				);

			// Return matrix with coordinates in real space. Coordinates need to be placed in local
			// numeration of indexes (they are translated later by mask)
			virtual void generateSpectralNodesRealCoordinates
			(
				vector< vector< array<double,2> > > & rSCoords,
				vector< array<double,2> > & rCoords
			)=0;


			// THIS METHOD SHOULD FILL UP NODES,WEIGHTS FIELDS, ACORDINGLY TO THIS ELEMENT TYPE DEFINITION
			virtual void generateSpectralNodesAndWeights
				(
				)=0;

			// THIS METHOD SHOULD FILL UP DERIVATVE MATIRX
			// (derivative marix-> d(l(ksi_i))/d(ksi_j), where l- shape function, ksi_i - spectral node at i index))
			virtual void generateDerivativeMatrix
				(
				)=0;

			// THIS METHOD NEED TO RETURN ELEMENTS FROM "MASK" THAT CORESPONDS TO "EDGE_ID'S" EDGE
			virtual vector<int> getEdgeGloblaNodesId
			(
			int edgeId,
			vector< vector<int> > & mask
			)=0;

			// THIS METHOD NEED TO RETURN ELEMENTS THAT LOCALY CORESPONDS TO "EDGE_ID'S" EDGE
			virtual vector< array<int,2> > getEdgeLocalNodesId
			(
				int edgeId
			)=0;

			//DEFINE THE SIZE OF ELEMENT MASK FPR ALL SPECTRAL NODES ASIGNED TO THIS (THE ARRAY)
			// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
			virtual void defineMaskArraySize
				(
				vector< vector<int> > & mask
				)=0;

			// ALLOCATE IN PROPER PLACE IN MASK ARRAY THE ELEMENT VERTEX NODES GIVEN GLOBAL ADRESS(GIVEN FROM MESH, NOT GENERATED AS SPECTRAL)
			//0-always left bottom, next vertex counterclockwise
			// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)

			virtual void allocateMaskArray_VertexNode
				(
				int vertexNr,
				int nodeGlobalAdress,
				vector< vector<int> > & addressArray
				)=0;

			// ALLOCATE INTO GLOBAL ADRESS MASK ARRAY SPECTRAL NODES ON EDGE(WITHOUT VERTEX NODE)
			// 0-bottom edge, and next counterclockwise. Nodes are given in this direction, means alwasy next node is this, which is counterclockwise
			// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
			virtual void allocateMaskArray_NodesInEdge
				(
				int edgeNr,
				vector<int> edgeGlobalAdress,
				vector< vector<int> > & addressArray
				)=0;

			// ALLOCATE INTO GLOBAL ADRESS MASK ARRAY SPECTRAL NODES INSIDE NODES
			// NODES SHOULD BE INDEXED FROM START ADDRESS TO SOME APROPRIATE VALUE
			// Nodes shell be numbered from left to right,and then up(--> means a[i...N][0]-bottom edge, and so on)
			virtual void allocateMaskArray_InteriorNodes
				(
				int startAddress,
				vector< vector<int> > & addressArray
				)=0;

			//Fill up the "H-Matrix"---> this is the evaluation of this parts of equation in laplacian discret week form:
			// H_Matrix[i][j][0] = (ksix*ksix)*J*W[i]*W[j];          H_Matrix[i][j][3] = (ksiy*ksiy)*J*weights[i]*weights[j];
			// H_Matrix[i][j][1] = (ksix*etax)*J*W[i]*W[j];			 H_Matrix[i][j][4] = (ksiy*etay)*J*weights[i]*weights[j];
			// H_Matrix[i][j][2] = (etax*etax)*J*W[i]*W[j];			 H_Matrix[i][j][5] = (etay*etay)*J*weights[i]*weights[j];
			// i,j - element spectral node index
			// W - weights
			// J - jacobian at i,j
			virtual void generateH_Matrix
				(
				vector< vector< array<double,6> > > & H_Matrix,
				vector< array<double,2> > & realCoords
				)=0;

			// Fill up the "M-Matrix"---> this is part of equation required for building the mass matrix. Its for is:
			// M_Matrix[i][j] = J*W.V[i]*W.V[j];
			// i,j - element spectral node index
			// W - weights
			// J - jacobian at i,j
			virtual void generateM_Matrix
				(
				vector< vector<double> > & M_Matrix,
				vector< array<double,2> > & realCoords
				)=0;

			// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
			// with scalar value-constant for all nodes
			// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
			//  i,j - coresponds to local i,j node)
			virtual void generateStiff_Matrix
				(
				vector< vector< vector< vector<double> > > > & Stiff_Matrix,
				double scalarValue,
				vector< vector< array<double,6> > > & H_Matrix
				)=0;

			//// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
			//// with tensor coefficient value-constant for all nodes
			//// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
			////  i,j - coresponds to local i,j node)
			//virtual void generateStiff_Matrix
			//(
			//  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
			//	array<double,4>& tensorValue,
			//	vector<vector<array<double,3>>>& H_Matrix
			//)=0;


			//// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
			//// with scalar coefficient field
			//// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
			////  i,j - coresponds to local i,j node)
			//virtual void generateStiff_Matrix
			//(
			//  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
			//	vector<vector<double>>& scalarNodalValues,
			//	vector<vector<array<double,3>>>& H_Matrix
			//)=0;

			//// Basing on constant parts of equation (H_Matrix, D_Matrix) generate stiffnes matrix
			//// With tensor coefficient filed
			//// (Stiff_Matrix[alpha][betha][i][j], where alpha,betha-coresponds to equation for local alpha,betha node
			////  i,j - coresponds to local i,j node)
			//virtual void generateStiff_Matrix
			//(
			//  vector<vector<vector<vector<double>>>>& Stiff_Matrix,
			//	vector<vector<array<double,4>>>& tensorNodalValues,
			//	vector<vector<array<double,3>>>& H_Matrix
			//)=0;


			//Calcualte coefficient value of neuman BC integral for nodes in apropriate edge
			// This value is of the form: int_edge(h*Phipq)dEdge = h*sqrt((dx/dksi)^2+(dy/dksi)^2)*wp*w0
			// This value multiplied by flux value on edge gives influence of Neuman BC to equation(known value)
			virtual void generateEdgeNeumanBCIntegralCoefficient
			(
				int edgeNumber,
				vector< array<double,2> > & elementCoordinates,
				vector<double> & resultNodalVal
			)=0;



			virtual ~SpectralElement
				(
				)
			{
			};

		protected:
			// RETRUN THE NUMBER OF INTERNIOR NODES ( THE COUNT OF NODES, WHERE NONE LIES ON ANY EDGE OR VERTEX )
			virtual int generateNumberOfInteriorNodes
				(
				)=0;

			// RETRUN NUMBERS(ON EACH EDGE RESPECTIVLY) OF INTERNIOR IN EDGE NODES ( THE COUNT OF NODES, WHERE NONE LIES VERTEX OR INSIDE )
			virtual vector<int> generateNumberOfInteriorInEdgesNodes
				(
				)=0;

		};
	}
}
#endif _SpectralElement_H
