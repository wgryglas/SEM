
#ifndef _Triangle_H_
#define _Triangle_H_


#include "RealElement.h"
#include "SpectralElement.h"
#include "utilities/Reference.h"

namespace SEM{ namespace mesh{

	class Triangle: public RealElement
	{
        REFERENCE_TYPE(Triangle)
	public:
		Triangle(Triangle& q);
        Triangle & operator=(const Triangle& a);
		
		Triangle
		(
            const std::vector<Vector > &m_nodes,
			SpectralElement& m_sEnt,
            const Matrix<int>::type &m_indexMask
		);


        friend std::ostream& operator<<(std::ostream&, Triangle& q);

	};


}//mesh
}//SEM





#endif //_Triangle_H_
