#ifndef BOUNDARYEDGE_H
#define BOUNDARYEDGE_H

#include "elements/RealElement.h"
#include "utilities/Reference.h"
#include "utilities/numArray.h"
#include "iomanagment/InfoStream.h"

namespace SEM { namespace mesh {

class Mesh;
/////////////////////////////////////////////
/// \class BoundaryEdge
/// For storing information about specific
/// boundary edge.
/// ----------------------------------------
/// note: Information stored in this class
/// shall be written by Mesh class instance
/// at the mesh build time. Mesh class is
/// currently fired of BoundaryEdge.
/////////////////////////////////////////////
class BoundaryEdge
{
    REFERENCE_TYPE(BoundaryEdge)

    friend class Mesh;

    int m_elementId;
    int m_edgeId;

    std::vector<int> m_nodesId;
    numArray<double> m_neumanBCPreValue;
    std::vector<boost::array<int,2> > m_nodesLocalId;
    std::vector<int> m_nodesLocalIdInMatrix;

public:
    /// Empty constructor shall not exist here, due to
    /// possible NULL-pointer exception with Mesh member variable.
    /// But for reading perpous (supported read from strem
    /// in this library only for classes with empty const.)
    /// this constructor is here, so user shall keep in mind
    /// that mesh must be set.
    //BoundaryEdge();

    BoundaryEdge();

    BoundaryEdge(const BoundaryEdge& other);

    BoundaryEdge& operator =(const BoundaryEdge& other);

    /// access nodes ids
    inline const std::vector<int> & gNodesIds() const { return m_nodesId; }

    /// access nodes local-element ids
    inline const std::vector<boost::array<int,2> > & lNodesIds() const { return m_nodesLocalId; }

    /// access nodes local-element ids in localElementMatrix
    inline const std::vector<int> & lNodesIdsInElementMatrix() const { return m_nodesLocalIdInMatrix;}

    /// access neumanBC values on edge ---> in the future it must be moved out of here into GeometricField in some general manner for BC
    inline const numArray<double> & neumanBCPreValue() const {return m_neumanBCPreValue; }

    /// access elemet id for this edge
    inline int elementId() const { return m_elementId; }

    /// access edge id relative to element
    inline int edgeId() const    { return m_edgeId; }
    

    /// access element where this edge belongs
    RealElement& element(Mesh & mesh);
    const RealElement& element(const Mesh & mesh) const;

    bool operator ==(const BoundaryEdge& other) const;
    bool operator !=(const BoundaryEdge& other) const;


    /// Operator for easy reading element and edge from list in stream
    friend std::istream& operator >>(std::istream& stream, BoundaryEdge &e);

    /// boundary edge info plotting --> right now shall not be used as feature of program, for debug case
    friend std::ostream& operator <<(std::ostream& stream, const BoundaryEdge &e);
};

//#include "Mesh.h"
//inline RealElement & BoundaryEdge::element() const { return *(m_mesh[m_elementId]); }


}//mesh
}//SEM


#endif // BOUNDARYEDGE_H
