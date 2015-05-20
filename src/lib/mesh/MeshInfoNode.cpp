#include"MeshInfoNode.h"
namespace SEM{ namespace mesh{

    MeshInfoNode::MeshInfoNode(int id) : m_globalWrittenNode(-1), id(id)
	{
    }
	
    void MeshInfoNode::setAsWritten(int globalWrittenNode)
	{
		this->written = true;
		this->m_globalWrittenNode = globalWrittenNode;
	}

} //mesh
} //SEM
