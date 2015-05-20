
#ifndef _Triangle_H_
#define _Triangle_H_


#include "RealElement.h"
#include "SpectralElement.h"

namespace SEM
{
namespace MESH
{
	class Triangle: public RealElement
	{
	public:

		Triangle(Triangle& q);
		void operator=(const Triangle& a);
		
		Triangle
		(
			vector<array<double,2>> nodes, 
			SpectralElement& sEnt,
			vector<vector<int>>& mask
		);

	};


	std::ostream& operator<<(std::ostream&, Triangle& q);
}
}





#endif _Triangle_H_