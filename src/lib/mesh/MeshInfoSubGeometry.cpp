#include"MeshInfoSubGeometry.h"

namespace SEM
{
namespace mesh
{
	MeshInfoSubGeometry::MeshInfoSubGeometry
	()
	:written(false)
	{};
	
	
	/*void MeshInfoSubGeometry::assignElement
	(MeshInfoElement& e)
	{
		elements.push_back(&e);
	};*/
	

	bool MeshInfoSubGeometry::isWritten
	()
	{
		return written;
	}


}
}