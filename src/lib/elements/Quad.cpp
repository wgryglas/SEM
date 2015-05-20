#include "Quad.h"

namespace SEM{ namespace mesh{

    Quad::Quad(Quad& q): RealElement(q.m_nodes, q.m_sEnt, q.m_indexMask)
	{
    }

	Quad::Quad
	(
        const std::vector<Vector> &nodes,
        SpectralElement & sEnt,
        const Matrix<int>::type & mask
	) : RealElement(nodes, sEnt, mask)
	{
    }
	
    Quad & Quad::operator= (const Quad& q)
	{
		m_nodes = q.m_nodes;
		m_sEnt = q.m_sEnt;

        return *this;
	}

    std::ostream& operator<< ( std::ostream& o, Quad& q )
	{
		o<<"============ Quad element ============"<<std::endl;
		o<<"Coordinates:"<<std::endl;
		for(int i=0;i<4;i++)
			o<<"["<<q.m_nodes[i][0]<<"; "<<q.m_nodes[i][1]<<"]"<<std::endl;
		o<<"======================================="<<std::endl;

		return o;
	}

}//mesh
}//SEM
