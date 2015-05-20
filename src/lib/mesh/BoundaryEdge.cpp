#include "BoundaryEdge.h"

#include "Mesh.h"

#include "iomanagment/ReadWriteStream.h"

namespace SEM{ namespace mesh{

BoundaryEdge::BoundaryEdge()
{
}

BoundaryEdge::BoundaryEdge(const BoundaryEdge &other)
    : m_elementId(other.m_elementId),
      m_edgeId(other.m_edgeId),
      m_nodesId(other.m_nodesId),
      m_nodesLocalId(other.m_nodesLocalId),
      m_neumanBCPreValue(other.m_neumanBCPreValue)
{
}

BoundaryEdge &BoundaryEdge::operator =(const BoundaryEdge &other)
{
      m_elementId=other.m_elementId;
      m_edgeId=other.m_edgeId;
      m_nodesId=other.m_nodesId;
      m_nodesLocalId=other.m_nodesLocalId;
      m_neumanBCPreValue=other.m_neumanBCPreValue;
}

RealElement & BoundaryEdge::element(Mesh & mesh) { return mesh[m_elementId]; }

const RealElement & BoundaryEdge::element(const Mesh &mesh) const { return mesh[m_elementId]; }

bool BoundaryEdge::operator ==(const BoundaryEdge &other) const
{
    return (m_nodesId==other.m_nodesId) && (m_elementId==other.m_elementId);
}

bool BoundaryEdge::operator !=(const BoundaryEdge &other) const
{
    return (m_nodesId!=other.m_nodesId) || (m_elementId!=other.m_elementId);
}


std::istream & operator >>(std::istream &stream, BoundaryEdge & e)
{
    boost::array<int,2> v;
    iomanagment::read<boost::array<int,2> >(stream, v);
    e.m_elementId = v[0];
    e.m_edgeId    = v[1];
}

std::ostream & operator <<(std::ostream &stream, const BoundaryEdge & e)
{
    stream<<"Boundary Edge::"<<endl
          <<"Element id="<<e.elementId()<<endl
          <<"Edge id in element="<<e.edgeId()<<endl;
}


}//mesh
}//SEM
