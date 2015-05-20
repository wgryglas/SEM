#include "MeshInfo.h"

namespace SEM
{
namespace mesh
{
	////---------------------------- GEO INFO ELEMENT --------------
	//MeshInfoElement::MeshInfoElement
	//()
	//{
	//}

	//void MeshInfoElement::addNode
	//(MeshInfoNode& n)
	//{
	//	nodes.push_back(n);
	//}

	//void MeshInfoElement::addEdge
	//(MeshInfoEdge& e)
	//{
	//	edges.push_back(e);
	//}

	//void MeshInfoElement::swapEdge
	//(MeshInfoEdge& from, MeshInfoEdge& to)
	//{
	//	edges.insert(find(edges.begin(),edges.end(),from),to);
	//}

	//bool MeshInfoElement::isQuad
	//()
	//{
	//	switch (edges.size())
	//	{
	//	case 4:
	//		return true;
	//	default:
	//		return false;
	//	}
	//}

	//bool MeshInfoElement::isTri
	//()
	//{
	//	switch (edges.size())
	//	{
	//	case 3:
	//		return true;
	//	default:
	//		return false;
	//	}
	//}

	////---------------------------- GEO INFO SUBGEOMETRY -----------
	//MeshInfoSubGeometry::MeshInfoSubGeometry
	//()
	//:written(false)
	//{};
	//
	//void MeshInfoSubGeometry::assignElement
	//(MeshInfoElement& e)
	//{
	//	elements.push_back(e);
	//};

	//bool MeshInfoSubGeometry::isWritten
	//()
	//{
	//	return written;
	//}

	////---------------------------- GEO INFO EDGE ----------------
	//MeshInfoEdge::MeshInfoEdge
	//(MeshInfoNode& node1,MeshInfoNode& node2, MeshInfoElement& e)
	//:node1(node1), node2(node2), interior(true)//, primary(true), written(false)
	//{
	//	assignElement(e);
	//};

	//void MeshInfoEdge::setAsBoundaryEdge
	//()
	//{
	//	interior = false;
	//	checked = true;
	//}

	//bool MeshInfoEdge::isSimilar
	//(MeshInfoEdge& other)
	//{
	//	/*--> only if order is not important, here we assume it should be done in speciefied order
	//	if(&node1==&other.node1 && &node2==&other.node2)
	//	{
	//		other.primary = false;
	//		checked = true;
	//		other.checked= true;
	//	}
	//	*/
	//	if(&node1==&other.node2 && &node2==&other.node1)
	//	{
	//		return true;
	//	}
	//	else
	//	{
	//		return false;
	//	}
	//}


	//void MeshInfoEdge::setAsChecked
	//()
	//{
	//	checked = true;
	//}

	//void MeshInfoEdge:: setAsBoundaryEdge
	//()
	//{
	//	interior = false;
	//}


	//bool MeshInfoEdge::isInterior
	//()
	//{
	//	return interior;
	//};
	//
	//void MeshInfoEdge::setAsWritten
	//(vector<int> globalWrittenInteriorNodes)
	//{
	//	this->written = true;
	//	this->globalWrittenInteriorNodes = globalWrittenInteriorNodes;
	//}
	///*
	//void MeshInfoEdge::setAsPrimary
	//()
	//{
	//	primary=true;
	//}
	//
	//void MeshInfoEdge::setAsSecondary
	//()
	//{
	//	primary= false;
	//}
	//
	//bool MeshInfoEdge::isPrimary
	//()
	//{
	//	return primary;
	//};
	//*/
	////---------------------------- GEO INFO NODE ----------------
	//MeshInfoNode::MeshInfoNode():globalWrittenNode(-1)
	//{};
	//
	//void MeshInfoNode::setAsWritten(int globalWrittenNode)
	//{
	//	this->written = true;
	//	this->globalWrittenNode = globalWrittenNode;
	//}

	//---------------------------- ADDITIONAL FUNCTION -----------
//	int *** generateGlobalAdresses
//	(vector<MeshInfoElement&> elements)
//	{
//	}

}
}
