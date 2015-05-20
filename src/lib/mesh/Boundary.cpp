#include "Boundary.h"

#include "Mesh.h"
#include "iomanagment/ReadWriteStream.h"

namespace SEM { namespace mesh {

Boundary::Boundary(const std::string &name, Mesh &mesh, size_t size):
    baseType(size), m_name(name), m_mesh(mesh)
{
}

Boundary::Boundary(const Boundary &other)
    : baseType(other), m_name(other.m_name), m_mesh(other.m_mesh)
{
}

Boundary &Boundary::operator =(const Boundary &other)
{
    baseType::operator =(other);
    m_name = other.m_name;
    m_mesh = other.m_mesh;
}

const RealElement &Boundary::element(size_t edgeInBoundary) const
{
    return at(edgeInBoundary).element(m_mesh);
}

RealElement &Boundary::element(size_t edgeInBoundary)
{
    return at(edgeInBoundary).element(m_mesh);
}

Vector Boundary::normal(size_t edge) const 
{
    return at(edge).element(m_mesh).normal(at(edge).edgeId());
}

SEM::numArrayIndexMapped< const numArray< Vector > > Boundary::nodes(std::size_t edge) const 
{
    return m_mesh.spectralNodes().slice(at(edge).gNodesIds());
}

std::istream &operator >>(std::istream &in, Boundary & out)
{
    iomanagment::read<std::vector<BoundaryEdge> >(in, out);
    return in;
}

std::ostream &operator <<(std::ostream &out, const Boundary & in)
{
    iomanagment::write<std::vector<BoundaryEdge> >(in, out);
    return out;
}



}//mesh
}//SEM
