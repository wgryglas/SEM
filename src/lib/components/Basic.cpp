#include "Basic.h"

#include <vector>

#include "iomanagment/ReadWriteStream.h"

namespace SEM
{

EigenVector::EigenVector(Scalar x, Scalar y)
    : VectorBase(x,y)
{
}

EigenVector::EigenVector(const EigenVector &v): VectorBase(v)
{
}

EigenVector & EigenVector::operator =(const EigenVector &v)
{
    if(*this!=v)
    {
        for(int i=0;i<SPACE_DIM;i++)
            (*this)(i)=v(i);
    }

    return *this;
}


std::ostream & operator <<(std::ostream & stream, const EigenVector & v)
{
    stream<<'(';
    for(int i=0;i<SPACE_DIM; i++)
    {
        stream<<v(i);

        if(i!=SPACE_DIM-1)
            stream<<' ';
    }
    stream<<')';
}

std::istream & operator >> (std::istream & stream, EigenVector & v)
{
    boost::array<Scalar,SPACE_DIM> input;
    SEM::iomanagment::read<boost::array<Scalar,SPACE_DIM> >(stream,input);

    v(0)=input[0];
    v(1)=input[1];

    return stream;
}


}//SEM
