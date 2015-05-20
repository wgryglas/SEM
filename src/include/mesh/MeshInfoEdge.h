
#ifndef _MeshInfoEdge_H_
#define _MeshInfoEdge_H_

#include<vector>

#include"MeshInfoSubGeometry.h"
#include"MeshInfoNode.h"


namespace SEM
{
namespace mesh
{
    class MeshInfoElement;

    class MeshInfoEdge: public MeshInfoSubGeometry
    {
    private: 
        MeshInfoNode& node1;
        MeshInfoNode& node2;
        
        MeshInfoElement& initialMaster;
        
        bool interior;
        bool checked;
        vector<int> globalWrittenInteriorNodes;

    public:
        MeshInfoEdge
        (
            MeshInfoNode& node1,
            MeshInfoNode& node2, 
            MeshInfoElement& element
        );

        bool isSimilar(MeshInfoEdge& other);

        void setAsBoundaryEdge();

        void setAsWritten(const vector<int> &globalWrittenInteriorNodes );

        void setAsChecked();
        
        bool isChecked();
        
        bool isInterior();
        
        MeshInfoElement& getInitialMasterElement();

        MeshInfoNode& getFirstNode();

        MeshInfoNode& getSecondNode();

        vector<int> getGlobalWrittenNodes();

    };
}
}

#endif //_MeshInfoEdge_H_
