#include "InletOutletVelocity.h"


namespace SEM {
namespace field {

    
std::string InletOutletVelocity::TYPE_NAME  = "inletOutletVelocity";

const std::string VELOCITY_DO_NOTHING= "velocityDoNothing";

InletOutletVelocity::InletOutletVelocity(const SEM::mesh::Boundary& edges) 
: PatchField< SEM::Vector >(edges)
{
}

InletOutletVelocity::~InletOutletVelocity()
{
}    
    
void InletOutletVelocity::read(const iomanagment::Dictionary& fieldDict)
{
    //Do nothing, it's just a kind of flag to inform that here comes inlet/outlet of velocity.
    //no value required
}

void InletOutletVelocity::write(iomanagment::Dictionary& fieldDict)
{
    //Write type only, values are uninportant
    iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
    
    iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
    bDict->add(typeEntry);
    
    if(!fieldDict.hasSubDictionary("boundary"))
        fieldDict.add(new iomanagment::Dictionary("boundary"));
    
    fieldDict.subDictionary("boundary").add(bDict);    
}

std::string InletOutletVelocity::typeName() const 
{
    return InletOutletVelocity::TYPE_NAME;
}
std::string InletOutletVelocity::typeFamily() const 
{
    return VELOCITY_DO_NOTHING;
}
 
    
}//field
}//SEM

