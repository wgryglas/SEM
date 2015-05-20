#ifndef _MeshInfoElement_H_
#define _MeshInfoElement_H_

#include <vector>

#include"MeshInfoNode.h"
#include"MeshInfoEdge.h"

#include<map>

using namespace std;

namespace SEM
{
namespace mesh
{
	class MeshInfoElement 
	{
        int m_id;
	public :
		vector<MeshInfoEdge*> edges;
		vector<MeshInfoNode*> nodes;

        MeshInfoElement(int m_id);

        void addNode(MeshInfoNode* n);
		
        void addEdge(MeshInfoEdge* e);

        void swapEdge( MeshInfoEdge& from, MeshInfoEdge& to);

        bool isQuad();

        bool isTri();

        inline int id() const {return m_id;}

	};

}
}





#endif //_MeshInfoElement_H_
