
#ifndef _SpectralElement_H_
#define _SpectralElement_H_

#include "SpectralElement.h"

namespace SEM{ namespace mesh{
    
    DEFINE_IMPLEMENTATION_FACTORY(SpectralElement,unsigned int)
	
    SpectralElement::SpectralElement(): NSpec(5),type("none")
	{
    }

    SpectralElement:: SpectralElement(unsigned int NSpec, const std::string &type)
        : NSpec(NSpec), type(type)
	{
	}

    int SpectralElement::getNumberOfInteriorNodes()
	{
		return numberOfInteriorNodes;
	}

    int SpectralElement::getNumberOfInteriorInEdgeNodes(int edgeNumber)
	{
		return numberOfInteriorInEdgesNodes[edgeNumber];
	}


    int SpectralElement::getSpectralSize() const
	{
		return NSpec;
	}

    std::vector<double> SpectralElement::getSpectralNodes()
	{
		return nodes;
    }

    std::vector<double> SpectralElement:: getSpectralWeights()
	{
		return weights;
    }
    
    numArray< Scalar > SpectralElement::computeInterpolationBarycentricWeights(const std::vector< Scalar >& nodes) 
    {
        numArray<Scalar> weights(nodes.size(),1.);
        
        for(size_t j=1; j<nodes.size(); ++j)
        {
            for(size_t k=0; k<j; ++k)
            {
                weights[k] *= nodes[k]-nodes[j];
                weights[j] *= nodes[j]-nodes[k];
            }    
        }
        
        return 1./weights;
    }
    

    std::ostream& operator <<( std::ostream &o, const SpectralElement &se)
    {
        o<<"Spectral Element"<<std::endl;
        o<<"nuber of spectral nodes="<<se.NSpec<<std::endl;
        o<<"type"<<se.type<<std::endl;
    }


} //mesh
} //SEM

#endif //_SpectralElement_H_
