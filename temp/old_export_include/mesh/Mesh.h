#ifndef _Mesh_H_
#define _Mesh_H_



#include <vector>
#include <array>
#include <map>

#include <Eigen/Dense>
#include <Eigen/StdVector>

#include "Quad.h"
#include "Triangle.h"
#include "RealElement.h"
#include "SpectralElement.h"
#include "SE_QLGL.h"
#include "SE_TLGL.h"
#include "MeshInfoElement.h"
#include "IOFile.h"
#include "InfoStream.h"

using namespace std;
using namespace SEM::IOStream;

namespace SEM
{
namespace MESH
{
	struct boundaryEdge
	{
		int element;
		int edge;
		vector<int> nodesId;
		vector<double> neumanBCPreValue;
		vector<array<int,2>> nodesLocalId;
	};

	class Mesh
	{
	public:
		std::vector<RealElement*> rEnts;
		std::map<string, vector<boundaryEdge>> boundaryInfo;
		std::vector<array<double,2>> sCoords; //coordinates of all nodes(spcetral and normal)
		int nREnts;
		int nSEnts;
		int nQuad;
		int nTri;
		int nNodes;
		int nSNodes;

		Mesh
		(
			IOStream::IOFile& meshFile,
			IOStream::IOFile& elementDefFile
		);

		~Mesh
		(
		)
		{
		};


		void writeSpectralNodes
		(
		);


	private:
		void generateMeshInfo
		(
			vector<MeshInfoElement*>& elements,
			map<string, vector<boundaryEdge>>& bInfo,
			int nNodes,
			int nEnts,
			std::vector<int>& elementType,
			std::vector<std::vector<int>>& elementNodes
		);

		void makeElementNodeAdressMask
		(
			vector<vector<int>>& mask,
			map<string,vector<boundaryEdge>>& bInfo,
			int* currNodeId,
			MeshInfoElement& eI,
			SpectralElement& sE
		);

		void evaluateSpectralCoordinates
		(
			vector<array<double,2>>& spectralC
		);

		void convertToCoordinates
		(
			vector<array<double,2>>& coodD,
			vector<array<string,2>>& coordS
		);

		void convertToElementNodesAndType
		(
			vector<vector<int>>& eNodesI,
			vector<int>& elementType,
			vector<vector<string>>& eNodesS
		);

		void createSpectralElements
		(
			vector<SpectralElement*>& sEnts,
			IOObject& elements,
			vector<vector<int>>& nodes
		);

		void convertToBoundaryInfo
		(
			map<string, vector<boundaryEdge>>& boundaryI,
			map<string,IOSubObject*> boundaryS
		);

		void evaluateNeumanBCCoefficients
		(
		);

	};


	ostream& operator<<
	(
		ostream& o,
		Mesh& mesh
	);
}
}

#endif _Mesh_H_
