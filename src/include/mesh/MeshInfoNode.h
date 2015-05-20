
#ifndef _MeshInfoNode_H_
#define _MeshInfoNode_H_

#include"MeshInfoSubGeometry.h"

namespace SEM{ namespace mesh{

class MeshInfoNode: public MeshInfoSubGeometry
	{
	private:
        int m_globalWrittenNode;
		
	public:
		int id;

        MeshInfoNode(int id);
		
        void setAsWritten(int m_globalWrittenNode);

        inline int globalWrittenNode() const {return m_globalWrittenNode;}
	};

}//mesh
}//SEM


#endif //_MeshInfoNode_H_
