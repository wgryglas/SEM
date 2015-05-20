#ifndef BOUNDARY_H
#define BOUNDARY_H

#include <string>
#include <vector>
#include <iostream>

#include "utilities/ListWrap.h"
#include "utilities/Reference.h"
#include "BoundaryEdge.h"
#include "elements/RealElement.h"

namespace SEM { namespace mesh {

/// forward Mesh decalaration
class Mesh;
////////////////////////////////////////////////
/// \brief The Boundary class
/// -------------------------------------------
/// Class containing mesh edges on specific
/// boundary which is described by "name" variable.
/// Derived from std::vector, so any std like
/// iteration over edges is possible.
/// --------------------------------------------
/// Additional element access by edge id in
/// boundary --> because here is kept mesh ref.,
/// not in each edge object
////////////////////////////////////////////////
class Boundary : public std::vector<BoundaryEdge>
{
    REFERENCE_TYPE(Boundary)
    typedef  std::vector<BoundaryEdge> baseType;

    /// \variable m_name - boundary name
    std::string m_name;

    /// \variable m_mesh - mesh reference
    Mesh &m_mesh;
    
    public:
        using baseType::operator [];


        /// base constructor
        Boundary(const std::string & name, Mesh &mesh, size_t size=0);

        /// copy constructor
        Boundary(const Boundary & other);

        /// assigment operator
        Boundary& operator =(const Boundary & other);

        /// destructor
        virtual ~Boundary(){}

        // Delegate access to edge- pass mesh reference
        /// \return real element refrence to specific in boundary edge
        /// const access element from edge in boundary
        const RealElement& element(size_t edgeInBoundary) const;

        /// \return real element refrence to specific in boundary edge
        /// access edge element
        RealElement &element(size_t edgeInBoundary);

        /// \return name of boundary
        inline std::string name() const { return m_name;}
        
        Vector normal(size_t edge) const;
        
        numArrayIndexMapped<const numArray<Vector> > nodes(size_t edge) const;
        
        
};

/// \param in - stream
/// \param oout - Boundary object
/// input from stream oprator - required because read/write method
/// can operate directly on template container, not on its derived
/// and instantied versions
std::istream& operator >>(std::istream& in, Boundary & out);

/// \param out - stream
/// \param in - Boundary object
/// input from stream oprator - required because read/write method
/// can operate directly on template container, not on its derived
/// and instantied versions
std::ostream& operator <<(std::ostream& out, const Boundary & in);


}//mesh
}//SEM






#endif // BOUNDARY_H
