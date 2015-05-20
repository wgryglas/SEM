
#ifndef _MeshInfoSubGeometry_H_
#define _MeshInfoSubGeometry_H_

#include<vector>

using namespace std;

namespace SEM
{
namespace  mesh
{
	class MeshInfoElement;

	class MeshInfoSubGeometry
	{
	protected:
		bool written;

	public:
		
		MeshInfoSubGeometry();

		bool isWritten();
	};



}
}

#endif //_MeshInfoSubGeometry_H_
