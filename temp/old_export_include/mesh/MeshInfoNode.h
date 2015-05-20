
#ifndef _MeshInfoNode_H_
#define _MeshInfoNode_H_

#include"MeshInfoSubGeometry.h"

namespace SEM
{
namespace MESH
{

class MeshInfoNode: public MeshInfoSubGeometry
	{
	private:
		int globalWrittenNode;
		
	public:
		int id;

		MeshInfoNode
		(
			int id
		);
		
		void setAsWritten
		(
			int globalWrittenNode
		);

		int getGlobalWrittenNode
		(
		);

	};

}
}


#endif _MeshInfoNode_H_