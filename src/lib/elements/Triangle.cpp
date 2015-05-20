#include "Triangle.h"

namespace SEM{ namespace mesh{

	Triangle::Triangle
	(
        const std::vector<Vector> &nodes,
		SpectralElement& sEnt,
        const Matrix<int>::type &mask
	) : RealElement(nodes, sEnt, mask)
	{
    }
	
    Triangle::Triangle(Triangle& q)
        : RealElement(q.m_nodes, q.m_sEnt, q.m_indexMask)
    {
    }

    Triangle & Triangle::operator=(const Triangle& q)
	{
		m_nodes = q.m_nodes;
		m_sEnt = q.m_sEnt;
        return *this;
	}

    std::ostream& operator<<(std::ostream& o, Triangle& q)
	{
		o<<"============ Triangle element ============"<<std::endl;
		o<<"Coordinates:"<<std::endl;
		for(int i=0;i<3;i++)
			o<<"["<<q.m_nodes[i][0]<<"; "<<q.m_nodes[i][1]<<"]"<<std::endl;
		o<<"======================================="<<std::endl;

		return o;
	}

} //mesh
} //SEM
