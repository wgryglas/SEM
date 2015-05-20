#include"MeshInfoElement.h"

#include<iostream>

#include "iomanagment/InfoStream.h"

namespace SEM
{
namespace mesh
{
    MeshInfoElement::MeshInfoElement(int id)
        : m_id(id), nodes(vector<MeshInfoNode*>(0)), edges(vector<MeshInfoEdge*>(0))
    {
    }

    void MeshInfoElement::addNode(MeshInfoNode* n)
    {
        nodes.push_back(n);
    }

    void MeshInfoElement::addEdge(MeshInfoEdge* e)
    {
        edges.push_back(e);
    }

    void MeshInfoElement::swapEdge(MeshInfoEdge& from,MeshInfoEdge& to)
    {
        vector<MeshInfoEdge*>::iterator i;
        int index=-1;
        int j=0;
        for(i=edges.begin();i!=edges.end();++i,j++)
            if(*i==&from)
                index=j;

        if(index==-1)
        {
            ErrorInFunction<<"Edge couldn't be swaped, "
                            <<"because object 'from' haven't "
                            <<"been found in MeshInfoElemnt object"<<iomanagment::endProgram;
        }
        else
        {
            edges[index]=&to;
        }
    }

    bool MeshInfoElement::isQuad()
    {
        return edges.size()==4;
    }

    bool MeshInfoElement::isTri()
    {
        return edges.size()==3;
    }


}
}
