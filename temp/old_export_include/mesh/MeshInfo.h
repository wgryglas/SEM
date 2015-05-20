
#ifndef _MeshInfo_H_
#define _MeshInfo_H_

#include <vector>

#include"MeshInfoElement.h"

using namespace std;

namespace SEM
{
namespace MESH
{
	////---------------------------- GEO INFO ELEMENT -----------
	//class MeshInfoElement 
	//{
	//public :
	//	
	//	vector<MeshInfoEdge&> edges;
	//	vector<MeshInfoNode&> nodes;

	//	MeshInfoElement();

	//	void addNode(MeshInfoNode& n);
	//	void addEdge(MeshInfoEdge& e);
	//	void swapEdge(MeshInfoEdge& from, MeshInfoEdge& to);

	//	bool isQuad();
	//	bool isTri();

	//};

	////---------------------------- GEO INFO SUBGEOMETRY ---
	//class MeshInfoSubGeometry
	//{
	//protected:
	//	bool written;

	//public:
	//	vector<MeshInfoElement&> elements;

	//	MeshInfoSubGeometry();

	//	void assignElement(MeshInfoElement& e);
	//	bool isWritten();
	//};

	////---------------------------- GEO INFO EDGE -----------
	//class MeshInfoEdge: public MeshInfoSubGeometry
	//{
	//private: 
	//	MeshInfoNode& node1, node2;
	//	bool interior;
	//	// bool primary;
	//	bool checked;
	//	vector<int> globalWrittenInteriorNodes;

	//public:
	//	MeshInfoEdge(MeshInfoNode& node1,MeshInfoNode& node2, MeshInfoElement& e);

	//	bool isSimilar(MeshInfoEdge& other);
	//	void setAsBoundaryEdge();

	//	void setAsWritten(vector<int> globalWrittenInteriorNodes );

	//	void setAsChecked();
	//	//void setAsPrimary();
	//	//void setAsSecondary();

	//	bool isChecked();
	//	bool isInterior();
	//	//bool isPrimary();
	//};

	////---------------------------- GEO INFO NODE -----------
	//class MeshInfoNode: public MeshInfoSubGeometry
	//{
	//private:
	//	int globalWrittenNode;
	//public:
	//	MeshInfoNode();
	//	void setAsWritten(int globalWrittenNode);
	//};



	//--------------------------- ADDITIONAL FUNCTIONS ---------
	int *** generateGlobalAdresses(vector<MeshInfoElement&> elements);

}
}

#endif _MeshInfo_H_