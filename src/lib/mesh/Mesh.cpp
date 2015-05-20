#include "Mesh.h"

#include "iomanagment/case.h"
#include "iomanagment/InfoStream.h"
#include "iomanagment/ReadWriteStream.h"

namespace SEM
{
namespace mesh
{
Mesh::Mesh(const boost::filesystem::path &meshPath,const boost::filesystem::path &elementsPath)
: nQuad(0),nTri(0),nSNodes(0)
{
    using namespace iomanagment;
    Dictionary mesh;
    meshPath>>mesh;

    Dictionary elements;
    elementsPath>>elements;

    buildMesh(mesh,elements);
}

Mesh::Mesh
	( 
        SEM::iomanagment::Dictionary &meshFile,
        SEM::iomanagment::Dictionary &elementDefFile
    ): nQuad(0),nTri(0),nSNodes(0)
	{
        buildMesh(meshFile, elementDefFile);
//        using namespace std;
//        using boost::array;
//        using namespace SEM::iomanagment;

//        Info<<"Creating mesh ..."<<endl;

//        Arrow<<"Reading coordinates"<<endl;
//        vector<Vector> coords;
//        meshFile.entry("coordinates")>>coords;
//        nNodes = coords.size();

//        Arrow<<"Reading nodes"<<endl;
//        vector<vector<int> > elementNodes;
//        meshFile.entry("elements")>>elementNodes;
//        nREnts = elementNodes.size();


//        Arrow<<"building spectral elements ..."<<endl;
//		vector<SpectralElement*> sElements(elementNodes.size());
//        createSpectralElements(elementDefFile, elementNodes, sElements);


//        Arrow<<"reading boundary mesh"<<endl;
//        readBoundaryInformation(meshFile.subDictionary("boundary"),m_boundaryMesh);

//        Arrow<<"building mesh topology... "<<endl;
//		vector<MeshInfoElement*> elementsInfo( elementNodes.size());
//        generateMeshInfo( elementsInfo, m_boundaryMesh, nNodes, nREnts, elementNodes );

//        Arrow<<"building internal spectral nodes in mesh... "<<endl;
//        vector<Matrix<int>::type> masks(elementsInfo.size());
		
//		for(int e=0;e< elementsInfo.size();e++)
//		{
//            makeElementNodeAdressMask(masks[e], m_boundaryMesh, &nSNodes,*(elementsInfo[e]),*(sElements[e]));
//		}

//        Arrow<<"building phisical elements ..."<<endl;
//        vector<Vector> eCorrd;

//        rEnts.resize(nREnts);
//		for( int e=0; e<elementNodes.size(); e++ )
//		{
//			switch ( elementNodes[e].size() )
//			{
//			case 3:
//				eCorrd.clear();
//				eCorrd.resize(3);

//                eCorrd[0] = coords[ elementNodes[e][0] ];
//                eCorrd[1] = coords[ elementNodes[e][1] ];
//                eCorrd[2] = coords[ elementNodes[e][2] ];

//                rEnts[e]=new Triangle( eCorrd, *sElements[e], masks[e] );
//				nTri++;

//				break;
//			case 4:

//				eCorrd.clear();
//				eCorrd.resize(4);

//                eCorrd[0] = coords[ elementNodes[e][0] ];
//                eCorrd[1] = coords[ elementNodes[e][1] ];
//                eCorrd[2] = coords[ elementNodes[e][2] ];
//                eCorrd[3] = coords[ elementNodes[e][3] ];

//                rEnts[e]= new Quad( eCorrd, *sElements[e], masks[e] );
//				nQuad++;

//				break;
//			default:
//                Warning<<"Found element with unsupported number of nodes, element will not be used, and program may crash"<<endl;
//				break;
//			}
//		}

//		//Fill spectral nodes coordinates
//		evaluateSpectralCoordinates(sCoords);

//		//Fill last information in boundary- eval. neuman BC coefficients:
//		evaluateNeumanBCCoefficients();
    }

    void Mesh::buildMesh(iomanagment::Dictionary &meshFile, iomanagment::Dictionary &elementDefFile)
    {
        using namespace std;
        using boost::array;
        using namespace SEM::iomanagment;

        Info<<"Creating mesh ..."<<endl;

        InfoArrow<<"Reading coordinates"<<endl;
        vector<Vector> coords;
        meshFile.entry("coordinates")>>coords;
        nNodes = coords.size();

        InfoArrow<<"Reading nodes"<<endl;
        vector<vector<int> > elementNodes;
        meshFile.entry("elements")>>elementNodes;
        nREnts = elementNodes.size();


        InfoArrow<<"building spectral elements ..."<<endl;
        vector<SpectralElement*> sElements(elementNodes.size());
        createSpectralElements(elementDefFile, elementNodes, sElements);


        InfoArrow<<"reading boundary mesh"<<endl;
        readBoundaryInformation(meshFile.subDictionary("boundary"),m_boundaryMesh);

        InfoArrow<<"building mesh topology... "<<endl;
        vector<MeshInfoElement*> elementsInfo( elementNodes.size());
        generateMeshInfo( elementsInfo, m_boundaryMesh, nNodes, nREnts, elementNodes );

        InfoArrow<<"building internal spectral nodes in mesh... "<<endl;
        vector<Matrix<int>::type> masks(elementsInfo.size());

        for(int e=0;e< elementsInfo.size();e++)
        {
            makeElementNodeAdressMask(masks[e], m_boundaryMesh, &nSNodes,*(elementsInfo[e]),*(sElements[e]));
        }

        InfoArrow<<"building phisical elements ..."<<endl;
        vector<Vector> eCorrd;

        rEnts.resize(nREnts);
        for( int e=0; e<elementNodes.size(); e++ )
        {
            switch ( elementNodes[e].size() )
            {
            case 3:
                eCorrd.clear();
                eCorrd.resize(3);

                eCorrd[0] = coords[ elementNodes[e][0] ];
                eCorrd[1] = coords[ elementNodes[e][1] ];
                eCorrd[2] = coords[ elementNodes[e][2] ];

                rEnts[e]=new Triangle( eCorrd, *sElements[e], masks[e] );
                nTri++;

                break;
            case 4:

                eCorrd.clear();
                eCorrd.resize(4);

                eCorrd[0] = coords[ elementNodes[e][0] ];
                eCorrd[1] = coords[ elementNodes[e][1] ];
                eCorrd[2] = coords[ elementNodes[e][2] ];
                eCorrd[3] = coords[ elementNodes[e][3] ];

                rEnts[e]= new Quad( eCorrd, *sElements[e], masks[e] );
                nQuad++;

                break;
            default:
                Warning<<"Found element with unsupported number of nodes, element will not be used, and program may crash"<<endl;
                break;
            }
        }

        //Fill spectral nodes coordinates
        evaluateSpectralCoordinates(sCoords);

        //Fill last information in boundary- eval. neuman BC coefficients:
        evaluateNeumanBCCoefficients();

    }


    void Mesh:: createSpectralElements
    (
        const iomanagment::Dictionary &elements,
        const Matrix<int>::type& elementNodes,
        std::vector<SpectralElement*> &sEnts
    )
    {
        using namespace std;
        using namespace iomanagment;

        int allocated=0;

        for(auto eDictPair : elements.subDictionaries())  
        {
            unsigned int spectralSize = eDictPair.second->entry("spectralSize").getValue<unsigned int>();
            SpectralElement* ptrSE= SpectralElement::Impl(eDictPair.second->entry("object").value())( spectralSize );
            
            if(eDictPair.second->entry("elements") == "allMatch")
            {
                for(int i=0; i<elementNodes.size(); ++i)
                {
                    if(elementNodes[i].size()==ptrSE->numberOfRealNodes())
                    {
                        sEnts[i] = ptrSE;
                        ++allocated;
                    }
                }
            }
            else
            {
                std::vector<int> indexes;
                eDictPair.second->entry("elements")>>indexes;
                for(int&  id : indexes)
                {
                    sEnts[id] = ptrSE;
                    ++allocated;
                }
            }
        }

        if(allocated != sEnts.size())
        {
            ErrorInFunction<<" Wrong definition of spectral elements. Possible not all elements "
                           <<" have subscribed it's own SpectralElement(look in file elements and "
                           <<" into elements field where should be defined for which elements is current SpectralElement"<<endProgram;
        }
    }
    

    void Mesh::readBoundaryInformation
    (
        const iomanagment::Dictionary& boundaryDict,
        BoundaryMesh & boundaryI
    )
    {
        using namespace std;
        using namespace iomanagment;

        const map<string,Dictionary*>& boundaryDicts = boundaryDict.subDictionaries();
        for(map<string,Dictionary*>::const_iterator itr=boundaryDicts.begin(); itr!=boundaryDicts.end();++itr)
        {
            Boundary tmpBoundary(itr->first,*this);

            //Read element and edge into boundary edges list-Boundary class
            itr->second->entry("edges")>>tmpBoundary;
            boundaryI.push_back(tmpBoundary);

    //            //Wrtie read data into BoundaryEdge list(no def. constructro, init with ref to Mesh)
    //            edges.resize( temp.size(), BoundaryEdge(this) );
    //            for(int i=0;i<temp.size(); ++i)
    //            {
    //                edges[i].m_elementId = temp[i][0];
    //                edges[i].m_edgeId = temp[i][1];
    //            }
    //
    //            boundaryI.insert( BoundaryMesh::value_type( itr->first,edges) );
    //            edges.clear();
    //            temp.clear();
        }
    }


    void Mesh::generateMeshInfo
    (
        std::vector<MeshInfoElement*>& elements,
        BoundaryMesh & bInfo,
        int nNodes,
        int nEnts,
        Matrix<int>::type & elementNodes
    )
    {
        using namespace std;

        //init components:
        vector<MeshInfoNode*> nodes(nNodes);
        vector<MeshInfoEdge*> edges;

        //build nodes
        for(int n=0;n<nNodes;n++ )
        {
            nodes[n]=new MeshInfoNode(n);
        }

        //subscribe elements to nodes, build edges and elements 
        int nodeFirst;
        int nodeSecond;
        for(int n=0;n<nEnts;n++)
        {
            elements[n]=new MeshInfoElement(n);

            for(nodeFirst=0,nodeSecond=1; nodeFirst< elementNodes[n].size();nodeFirst++,nodeSecond++)
            {
                elements[n]->addNode(nodes[elementNodes[n][nodeFirst]] );

                //Swith second node index to fill up loop over all edges
                if(nodeSecond==elementNodes[n].size())
                    nodeSecond = 0;

                edges.push_back(new MeshInfoEdge(*nodes[elementNodes[n][nodeFirst]],*nodes[elementNodes[n][nodeSecond]],*(elements[n]) ) );
                elements[n]->addEdge(edges.back()); //Here are duplicated edges in connected elements. They are only temp.
            }
        }

        // assigne boundary info for edges
        // ( !!!notice that at any boundary there will never be duplicated edge):
        BoundaryMesh::iterator boundary;
        vector< BoundaryEdge >::iterator edgeAddress;
        for(boundary = bInfo.begin();boundary!=bInfo.end();++boundary)
        {
            for(edgeAddress = boundary->begin(); edgeAddress!=boundary->end(); ++edgeAddress)
            {
                elements[edgeAddress->elementId()]->edges[edgeAddress->edgeId()]->setAsBoundaryEdge();
            }
        }

        // find similar edges
        vector<MeshInfoEdge*>::iterator itr;
        vector<MeshInfoEdge*>::iterator itr2;
        for(itr=edges.begin(); itr!=edges.end(); ++itr)
        {
            if( (*itr)->isInterior() )
            {
                for( itr2=edges.begin(); itr2!=edges.end(); ++itr2)
                {
                    if( !(*itr2)->isChecked() && (*itr2)->isInterior() )
                    {
                        if( (*itr)->isSimilar(*(*itr2)) )
                        {
                            (*itr)->setAsChecked();
                            (*itr2)->setAsChecked();
                            //replace duplicated edge in second element
                            (*itr2)->getInitialMasterElement().swapEdge(*(*itr2),*(*itr));
                        }
                    }
                }
            }
        }
    }


    void  Mesh::makeElementNodeAdressMask
    (
        Matrix<int>::type & mask,
        BoundaryMesh & bInfo,
        size_t* currNodeId,
        MeshInfoElement& eI,
        SpectralElement& sE
    )
    {
        using namespace std;

        sE.defineMaskArraySize(mask);
        
        //Write nodes
        for(int i=0;i<eI.nodes.size();i++)
        {
            if(eI.nodes[i]->isWritten())
            {
                sE.allocateMaskArray_VertexNode(i,eI.nodes[i]->globalWrittenNode(),mask);
            }
            else
            { //not yet written node
                sE.allocateMaskArray_VertexNode(i,*currNodeId,mask);
                eI.nodes[i]->setAsWritten(*currNodeId);
                (*currNodeId)++;
            }
        }

        //Write interior(in edge) nodes:
        for(int i=0;i<eI.edges.size();i++)
        {
            if(eI.edges[i]->isWritten())
            {
                sE.allocateMaskArray_NodesInEdge(i, eI.edges[i]->getGlobalWrittenNodes(), mask);
            }
            else
            {  //not yet written egde
                vector<int> newAddress(sE.getNumberOfInteriorInEdgeNodes(i));
                for(int e=0;e<sE.getNumberOfInteriorInEdgeNodes(i);e++)
                {
                    newAddress[e]=*currNodeId;
                    (*currNodeId)++;
                }
                sE.allocateMaskArray_NodesInEdge(i,  newAddress, mask);
                eI.edges[i]->setAsWritten(newAddress);
            }
        }

        //Write interior(in element) nodes:
        sE.allocateMaskArray_InteriorNodes(*currNodeId,mask);
        (*currNodeId)+=sE.getNumberOfInteriorNodes();
        

        //Check if element is boundary kind? If yes place its edges nodes adrresses into boundaryInfo structure
    //        map<string, vector<int> > bTemp;
        for(BoundaryMesh::iterator bItr=bInfo.begin();bItr!=bInfo.end(); ++bItr)
        {
            for(int j=0; j<bItr->size(); ++j)
            {
                if( (*bItr)[j].elementId()==eI.id() )
                {
                    (*bItr)[j].m_nodesId = sE.getEdgeGloblaNodesId( (*bItr)[j].edgeId(), mask );
                    (*bItr)[j].m_nodesLocalId = sE.getEdgeLocalNodesId( (*bItr)[j].edgeId() );
                    (*bItr)[j].m_nodesLocalIdInMatrix = sE.getEdgeLocalNodesIdInMatrix((*bItr)[j].edgeId());
                }
            }
        }

    }

    void Mesh::evaluateSpectralCoordinates(std::vector<Vector> & spectralC)
    {
        using namespace std;
        using boost::array;

        spectralC.resize(nSNodes);

        Matrix<Vector>::type temp;

        for(vector<RealElement*>::iterator ele = rEnts.begin(); ele!=rEnts.end();ele++)
        {
            (*ele)->spectralElement().generateSpectralNodesRealCoordinates(temp, (*ele)->nodes());
            
            const Matrix<int>::type &mask = (*ele)->indexMask();

            for(int i=0;i<mask.size();++i)
            {
                for(int j=0;j<mask[i].size();j++)
                {
                    spectralC[ mask[i][j] ] = temp[i][j];
                }
            }
        }

    }

    void Mesh::evaluateNeumanBCCoefficients()
    {
        using namespace std;

        for(BoundaryMesh::iterator bItr = m_boundaryMesh.begin(); bItr!=m_boundaryMesh.end();++bItr)
        {
            for(vector<BoundaryEdge>::iterator edgeItr = bItr->begin(); edgeItr!=bItr->end(); ++edgeItr)
            {
                edgeItr->element(*this).spectralElement().generateEdgeNeumanBCIntegralCoefficient
                                                    (edgeItr->edgeId(), edgeItr->element(*this).nodes(), edgeItr->m_neumanBCPreValue);
    //                rEnts[edge->element]->spectralElement().generateEdgeNeumanBCIntegralCoefficient(edge->edge, rEnts[edge->element]->nodes(), edge->neumanBCPreValue);
            }
        }
    }

    void Mesh::writeSpectralNodes(const boost::filesystem::path &filePath) const 
    {
        using namespace iomanagment;
        
        DictEntry* coordinatesEntry = new DictEntry("coordinates");
        
        *coordinatesEntry<<sCoords;
        
        DictEntry* elementsEntry = new DictEntry("elements");
        
        std::vector< std::vector<int> > elNodes;
        elNodes.reserve(size());
        for(int i=0; i<size(); ++i)
            elNodes.push_back(rEnts[i]->indexVectorMask());
        
        *elementsEntry<<elNodes;
        
        Dictionary nodesDict("spctralNodes");
        nodesDict.add(coordinatesEntry);
        nodesDict.add(elementsEntry);
        
        filePath<<nodesDict;

    }

    Boundary &Mesh::boundary(const string &name)
    {
        BoundaryMesh::iterator itr = m_boundaryMesh.begin();

        for(itr; itr!=m_boundaryMesh.end(); ++itr)
        {
            if(itr->name() == name)
                return *itr;
        }

        ErrorInFunction<<"Mesh boundary access failure, \n"
                        <<"mesh don't contain boundary named "<<name
                        <<iomanagment::endProgram;
    }

    std::ostream& operator<<(std::ostream& o,Mesh& mesh)
    {
        o<<"Mesh object:"<<std::endl;
        o<<"number of elements:\t"<<mesh.nREnts<<std::endl;
        o<<"number of tri elements:\t"<<mesh.nTri<<std::endl;
        o<<"number of quad elements:\t"<<mesh.nQuad<<std::endl;
        o<<"number of nodes:\t"<<mesh.nNodes<<std::endl;
        o<<"number of spectral nodes:\t"<<mesh.nSNodes<<std::endl;

        return o;
    }

}//MESH
}//SEM
