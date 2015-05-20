
#include "TCAvgPressure.h"

namespace SEM { namespace field {
    
    
    std::string TCAvgPressure::TYPE_NAME = "timeChangingAvgPressure";

//     TCAvgPressure::TCAvgPressure(const mesh::Boundary& edges) :
//     AvgPressure(edges)
//     {
//     }
    TCAvgPressure::TCAvgPressure(const iomanagment::Register& reg, const mesh::Boundary& edges) 
    : AvgPressure(reg,edges)
    {
    }
        
    TCAvgPressure::~TCAvgPressure() 
    {
    }

    std::string TCAvgPressure::typeName() const 
    {
        return TCAvgPressure::TYPE_NAME;
    }
    
    void TCAvgPressure::read(const iomanagment::Dictionary& fieldDict) 
    {
        m_values.update(fieldDict.subDictionary("boundary").subDictionary(name()).entry("value"));
        m_avgValue=m_values(Case::time().time() );
    }
    
    void TCAvgPressure::write(iomanagment::Dictionary& fieldDict) 
    {
        iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
        
        iomanagment::DictEntry *valEntry =new iomanagment::DictEntry("value");
        m_values.writeValues(*valEntry);
        
        iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
        bDict->add(typeEntry);
        bDict->add(valEntry);
        
        if(!fieldDict.hasSubDictionary("boundary"))
            fieldDict.add(new iomanagment::Dictionary("boundary"));
        
        fieldDict.subDictionary("boundary").add(bDict);
    }
    
    void TCAvgPressure::timeChanged(Scalar time) 
    {
        AvgPressure::timeChanged(time);
        m_avgValue=m_values(time );
    }
    
    
}//field
}//SEM