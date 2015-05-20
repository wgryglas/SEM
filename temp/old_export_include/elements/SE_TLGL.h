
#ifndef _SE_TLGL_H_
#define _SE_TLGL_H___

#include "SpectralElement.h"


namespace SEM
{
namespace MESH
{

	class SE_TLGL: public SpectralElement
	{
		public:
			SE_TLGL
			(
			);

			SE_TLGL
			(
			int NSpec
			);

			void generateSpectralNodesAndWeights
			(
			);

			void generateSpectralNodesRealCoordinates
			(
				vector<vector<array<double,2>>>& rSCoords, 
				vector<array<double,2>>& rCoords
			);


			void generateDerivativeMatrix
			(
			);

			vector<int> getEdgeGloblaNodesId
			(
				int edgeId,
				vector<vector<int>>& mask
			);

			vector<array<int,2>> getEdgeLocalNodesId
			(
				int edgeId
			);

			void defineMaskArraySize
			(
				vector<vector<int>>& mask
			);

			void allocateMaskArray_VertexNode
			(
				int vertexNr, 
				int nodeGlobalAdress,
				vector<vector<int>>& addressArray
			);

			void allocateMaskArray_NodesInEdge
			(
				int edgeNr,
				vector<int> edgeGlobalAdress,
				vector<vector<int>>& addressArray
			);

			void allocateMaskArray_InteriorNodes
			(
				int startAddress,
				vector<vector<int>>& addressArray
			);

			void generateH_Matrix
			(
				vector<vector<array<double,6>>>& H_Matrix,
				vector<array<double,2>>& realCoords 
			);

			void generateM_Matrix
			(
				vector<vector<double>>& M_Matrix,
				vector<array<double,2>>& realCoords
			);

			void generateStiff_Matrix
			(
				vector<vector<vector<vector<double>>>>& Stiff_Matrix,
				double scalarValue,
				vector<vector<array<double,6>>>& H_Matrix
			);

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
				vector<array<double,2>>& elementCoordinates,
				vector<double>& resultNodalVal
			);


		protected:
			int generateNumberOfInteriorNodes
				(
				);

			vector<int> generateNumberOfInteriorInEdgesNodes
				(
				);

			virtual ~SE_TLGL
				(
				);

		private:

			void qdqLN
				(
				double x,
				int N,
				double*q,
				double*dq,
				double*LN
				);

			double WJ
				(
				int j
				);


		};


		std::ostream& operator<<
			(
			std::ostream& o,
			SE_TLGL& e
			);

}
}



#endif _SE_TLGL_H_