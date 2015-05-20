#ifndef _Mesh_H_
#define _Mesh_H_

#include <vector>
#include <map>

#include "boost/array.hpp"
#include "boost/filesystem/path.hpp"

#include "elements/Quad.h"
#include "elements/Triangle.h"
#include "elements/RealElement.h"
#include "elements/SpectralElement.h"


#include "MeshInfoElement.h"
#include "BoundaryEdge.h"
#include "Boundary.h"

#include "iomanagment/DictEntry.h"
#include "iomanagment/Dictionary2.h"
#include "iomanagment/InfoStream.h"

#include "utilities/Reference.h"

namespace SEM{ namespace mesh{

    class Mesh
	{
        REFERENCE_TYPE(Mesh)
	public:
        /// \typedef BoundaryMesh - type for holding all boundary information
        typedef std::vector<Boundary> BoundaryMesh;
        /// \variable rEnts - all elements (real elements, spectral are only inner type in real) in mesh
		std::vector<RealElement*> rEnts;

        typedef std::vector<RealElement*>::iterator iterator;
        typedef std::vector<RealElement*>::const_iterator const_iterator;

        /// \variable m_boundaryMesh - collection of Boundary objects in mesh
        BoundaryMesh m_boundaryMesh;

        ///  \variable sCoords - coordinates of all nodes(spcetral and normal)
        numArray<Vector> sCoords;

        //Sizes - just information easy to access
        /// \variable nREnts - number of RealElements in mesh
        size_t nREnts;
        /// \variable nSEnts - number of SpectralElements in mesh
        size_t nSEnts;
        /// \variable nQuad - number of Quad-type elements in mesh
        size_t nQuad;
        /// \variable nTri - number of Tri-type elements in mesh
        size_t nTri;
        /// \variable nNodes - number of nodes in mesh (nodes provided by file)
        size_t nNodes;
        /// \variable nSNodes - number of all nodes in mesh (nodes generated as spcetral for all elements)
        size_t nSNodes;


        /// \brief Mesh constructor, which reads mesh from path provided by case
        Mesh(const boost::filesystem::path &meshPath,const boost::filesystem::path &elementsPath);

        /// \param meshFile - dictionary read from "mesh" file
        /// \param elementDefFile- dictionary read from "elements" file
        /// constructor from 2 neccessery files dictionaries
        Mesh( iomanagment::Dictionary & meshFile, iomanagment::Dictionary & elementDefFile );


        /// destructor
        virtual ~Mesh()
		{
        }

        /// \brief method for writing nodes into file
        /// \param meshDir path to directory where nodes shall be placed
        void writeSpectralNodes(const boost::filesystem::path & filePath) const;

        /// \param index- global element index
        /// \return element
        /// access operator to element
        inline RealElement& operator [](size_t index) { return *(rEnts[index]); }

        /// \brief begin iterator over elements pointers
        /// \return vector<RealElement*>::iterator
        inline iterator begin() { return rEnts.begin();}
        /// \brief begin const_iterator over elements pointers
        /// \return vector<RealElement*>::const_iterator
        inline const_iterator begin() const { return rEnts.begin();}

        /// \brief end iterator over elements pointers
        /// \return vector<RealElement*>::iterator
        inline iterator end() { return rEnts.end();}

        /// \brief end iterator over elements pointers
        /// \return vector<RealElement*>::const_iterator
        inline const_iterator end() const { return rEnts.begin();}

        /// \param index- global element index
        /// \return element
        /// const access operator to element
        inline const RealElement& operator [](size_t index) const { return *(rEnts[index]); }

        /// \return BoundaryMesh obj. ref - reference to list of all Boundaries in mesh
        inline BoundaryMesh & boundaryMesh() { return m_boundaryMesh; }

        inline size_t size() const {return rEnts.size(); }

        /// \param name - boundary name
        /// \return Boundary object ref. assigned for string name
        Boundary & boundary(const std::string & name);

        /// \brief getter for spectral coordinates
        /// \return numerical list of coordinates. List can be easy mapped to elemetnal nodes.
        inline const numArray<Vector> & spectralNodes() const { return sCoords; }
        
        /// \param name - boundary name
        /// \return Boundary object ref. assigned for string name
        /// const method version.
        inline const BoundaryMesh & boundaryMesh() const { return m_boundaryMesh; }

        inline const Boundary & boundaryAt(const std::string& name) const
        {
            for(BoundaryMesh::const_iterator bM = m_boundaryMesh.begin();
                                             bM!= m_boundaryMesh.end();
                                             ++bM)
            {
                if( name == bM->name() )
                {
                   return *bM;
                }
            }

            ErrorInFunction<<"There is no boundary named "<<name
                           <<" in mesh"<<iomanagment::endProgram;

            return m_boundaryMesh[0];
        }

        /// \brief nodesNumber - getter for nodes number
        /// \return number of all nodes in mesh
        inline int nodesNumber() const { return nSNodes; }


	private:
        /// perfom mesh building from 2 dictionaries defining whole mesh
        void buildMesh(iomanagment::Dictionary &meshFile, iomanagment::Dictionary &elementDefFile);


        /// collect mesh infromations
        void generateMeshInfo
        (
            std::vector<MeshInfoElement*>& elements,
            BoundaryMesh & bInfo,
            int nNodes,
            int nEnts,
            std::vector<std::vector<int> > & elementNodes
        );

        /// create node adress mask for each element
        void makeElementNodeAdressMask
        (
            Matrix<int>::type& mask,
            BoundaryMesh & bInfo,
            size_t* currNodeId,
            MeshInfoElement& eI,
            SpectralElement& sE
        );

        /// Calculate internal coordinates for each spectral node
        void evaluateSpectralCoordinates( std::vector<Vector> & spectralC );

        /// assign spectral elements for each real element
        /// spectral element instances can be asigned for many elements
        /// runtime chose spectral element from file specification
        /// !!!!! in fututre it shall be done by RUN_TIME_SELECTION_TABLE
        /// --> object-creator which will return appropriate instance of specific abstract class
        void createSpectralElements
        (
            const iomanagment::Dictionary &elements,
            const Matrix<int>::type& elementNodes,
            std::vector<SpectralElement*> &sEnts
        );


        /// \param boundaryDict - dictionary named "boundary" in mesh file
        /// read boundary information - elementId and edgeId
        void readBoundaryInformation
        (
            const iomanagment::Dictionary& boundaryDict,
            BoundaryMesh & boundaryI
        );

        /// assign coefficient values for neuman boundary condition
        /// those coefficient comes from integration, apropriate for
        /// thise kinde of BC, function over boundary element edge.
        /// Those coeffs are multiplied by neuman value on edge to
        /// apply influence of this bc to solution
        void evaluateNeumanBCCoefficients();

    };

    std::ostream& operator<<( std::ostream& o, Mesh& mesh );


}//mesh
}//SEM

#endif //_Mesh_H_
