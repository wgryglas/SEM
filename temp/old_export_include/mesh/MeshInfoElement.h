#ifndef _MeshInfoElement_H_
#define _MeshInfoElement_H_

#include <vector>

#include"MeshInfoNode.h"
#include"MeshInfoEdge.h"

#include<map>

using namespace std;

namespace SEM
{
namespace MESH
{
	class MeshInfoElement 
	{
	public :
		int id;
		vector<MeshInfoEdge*> edges;
		vector<MeshInfoNode*> nodes;

		MeshInfoElement
		(
			int id 
		);

		void addNode
		(
			MeshInfoNode* n
		);
		
		void addEdge
		(
			MeshInfoEdge* e
		);

		void swapEdge
		(
			MeshInfoEdge& from, 
			MeshInfoEdge& to
		);

		bool isQuad
		(
		);

		bool isTri
		(
		);

	};

}
}





#endif _MeshInfoElement_H_