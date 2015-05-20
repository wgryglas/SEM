
#ifndef _Quad_H_
#define _Quad_H_

#include "RealElement.h"
#include "SpectralElement.h"

namespace SEM
{
namespace MESH
{
	class Quad: public RealElement
	{
	public:

		Quad(Quad& q);
		void operator=(const Quad& a);
		
		Quad
		(
			vector<array<double,2>> nodes, 
			SpectralElement& sEnt,
			vector<vector<int>>& mask
		);

	};


	std::ostream& operator<<(std::ostream&, Quad& q);
}
}
#endif _Quad_H_