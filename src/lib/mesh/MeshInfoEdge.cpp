#include"MeshInfoEdge.h"
#include"MeshInfoElement.h"
#include <algorithm>

namespace SEM
{
namespace mesh
{
	MeshInfoEdge::MeshInfoEdge
	(
		MeshInfoNode& node1, 
		MeshInfoNode& node2, 
		MeshInfoElement& e
	):node1(node1), node2(node2), checked(false), interior(true), initialMaster(e)//, primary(true), written(false)
	{
	};
	
	MeshInfoNode&  MeshInfoEdge::getFirstNode
	(
	)
	{
		return node1;
	}

	MeshInfoNode&  MeshInfoEdge::getSecondNode
	(
	)
	{
		return node2;
	}


	void MeshInfoEdge::setAsBoundaryEdge
	(
	)
	{
		interior = false;
		checked = true;
	}

	bool MeshInfoEdge::isSimilar
	(
		MeshInfoEdge& other
	)
	{
		/*--> only if order is not important, here we assume it should be done in speciefied order
		if(&node1==&other.node1 && &node2==&other.node2)
		{
			other.primary = false;
			checked = true;
			other.checked= true;
		}
		*/
		if(&node1==&other.node2 && &node2==&other.node1)
		{
			return true;
		}
		else
		{
			return false;
		}
	}


	void MeshInfoEdge::setAsChecked
	(
	)
	{
		checked = true;
	}

	bool MeshInfoEdge::isInterior
	(
	)
	{
		return interior;
	};
	
	bool MeshInfoEdge::isChecked
	(
	)
	{
		return checked;
	}


	void MeshInfoEdge::setAsWritten
	(
		const vector<int> &globalWrittenInteriorNodes
	)
	{
		this->written = true;
		this->globalWrittenInteriorNodes = globalWrittenInteriorNodes;
		std::reverse(this->globalWrittenInteriorNodes.begin(),this->globalWrittenInteriorNodes.end());
	}

	MeshInfoElement& MeshInfoEdge::getInitialMasterElement
	(
	)
	{
		return initialMaster;
	}


	vector<int> MeshInfoEdge::getGlobalWrittenNodes
	(
	)
	{
		return globalWrittenInteriorNodes;
	}

}
}