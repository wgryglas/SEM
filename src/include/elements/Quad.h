
#ifndef _Quad_H_
#define _Quad_H_

#include "RealElement.h"
#include "SpectralElement.h"
#include "utilities/Reference.h"

namespace SEM{ namespace mesh {
	class Quad: public RealElement
	{
        REFERENCE_TYPE(Quad)
	public:

		Quad(Quad& q);
        Quad & operator=(const Quad& a);
		
		Quad
		(
            const std::vector<Vector> &m_nodes,
            SpectralElement &m_sEnt,
            const Matrix<int>::type &m_indexMask
        );


        friend std::ostream& operator<<(std::ostream&, Quad& q);
	};


} //mesh
} //SEM
#endif //_Quad_H_
