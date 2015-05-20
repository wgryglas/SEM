#ifndef LIST_H
#define LIST_H

#include <ostream>
#include <vector>
#include <valarray>

#include "Eigen/Dense"
#include "utilities/TypeDefs.h"
#include "iomanagment/ReadWriteStream.h"

namespace SEM {

/////////////////////////////////////////
/// TypeProvider- interface to allow using
/// more complex containers with multiple
/// class templates
/////////////////////////////////////////
template <typename T, template<class> class TypeProvider>
class ListWrap
{
public:
    typedef typename TypeProvider<T>::type Container;

protected:
    Container m_data;

public:
    ListWrap(size_t size=0): m_data(size)
    {
    }

    ListWrap(const Container& data): m_data(data)
    {
    }

    ListWrap<T,TypeProvider>& operator =(const ListWrap<T,TypeProvider>& l)
    {
        m_data = l.m_data;
        return *this;
    }

    ListWrap<T,TypeProvider>& operator =(const Container& data)
    {
        m_data = data;
        return *this;
    }

    Container& data()
    {
        return m_data;
    }

    const Container& data() const
    {
        return m_data;
    }

    size_t size() const
    {
        return m_data.size();
    }

    void resize(size_t size)
    {
        m_data.resize(size);
    }

    T& operator [](size_t loc)
    {
        return m_data[loc];
    }

    const T& operator [](size_t loc) const
    {
        return m_data[loc];
    }

    virtual ~ListWrap()
    {
    }
};


/// Those operators are required due to fact
/// that read<T>/write<T> function resolves
/// read/write algorithm by template parameters.
/// When this class would be subclassed, and new
/// class will have diffrent template arguments,
/// then it won't be possible to read/write this.
template< typename T, template<class> class TypeProvider>
std::istream & operator >> (std::istream & in, ListWrap<T,TypeProvider> & l)
{
    SEM::iomanagment::read<ListWrap<T,TypeProvider> >(in,l);
    return in;
}

template< typename T, template<class> class TypeProvider>
std::ostream & operator << (std::ostream& out,const ListWrap<T,TypeProvider> & l)
{
    SEM::iomanagment::write<ListWrap<T,TypeProvider> >(l,out);
    return out;
}

} //SEM






#endif // LIST_H
