
#ifndef _EquationMatrix_H_
#define _EquationMatrix_H_

#include<vector>

#include "components/Basic.h"
#include "utilities/numArray2D.h"



/// \namespace las = linear algebra solution
namespace SEM { namespace las {

typedef std::vector<numArray2D<Scalar> > ElMatrix;
typedef std::vector<ElMatrix> SEMMatrix;
typedef SEM::numArray<Scalar> ElVector;
typedef std::vector<ElVector> SEMVector;


inline void print(const SEMMatrix & matrix, std::ostream & o = std::cout)
{
    int dim = 0;
    for (const ElMatrix & dMat: matrix)
    {    
        o <<" ------------ Field dimmension "<<dim<<"-----------"<<std::endl;
        int element = 0;
        for(const numArray2D<Scalar> & eMat : dMat)
        {
            o << " --------------- Element "<< element <<"-------------------"<<std::endl;
            
            for(int i=0; i<eMat.size(); ++i)
            {
                    for(int j=0; j<eMat[i].size(); ++j)
                        o << eMat[i][j] <<"  ";
                    
                    o << std::endl;
            }
            ++element;
        }
        ++dim;
    }
}

}//las
}//SEM


#endif //_EquationMatrix_H_
