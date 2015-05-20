
#include "AvgPressure.h"
#include "iomanagment/case.h"
#include "elements/RealElement.h"
#include "fieldMath/TimeExtrapolation.h"
#include "GeometricField.h"

namespace SEM { namespace field {
    
    
const std::string PRESSURE_DO_NOTHING="pressureDoNothing";
    
std::string AvgPressure::TYPE_NAME = "avgPressure";
 
// AvgPressure::AvgPressure(const SEM::mesh::Boundary& edges) 
// : baseType(edges)
// {
// }

AvgPressure::AvgPressure(const iomanagment::Register& reg, const mesh::Boundary& edges) 
:baseType(edges), m_velocity(*reg.object<GeometricField<Vector> >("U") ),
m_viscosity(Case::material()["nu"])
{
}

AvgPressure::~AvgPressure()
{
}

std::string AvgPressure::typeName() const
{
    return AvgPressure::TYPE_NAME;
}

std::string AvgPressure::typeFamily() const 
{
    return DIRICHLET;
}

void AvgPressure::read(const iomanagment::Dictionary& fieldDict)
{
    fieldDict.subDictionary("boundary").subDictionary(name()).entry("value")>>m_avgValue;
}

void AvgPressure::write(iomanagment::Dictionary& fieldDict)
{
    iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
    
    iomanagment::DictEntry *valEntry =new iomanagment::DictEntry("value");
    *valEntry<<m_avgValue;
    
    iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
    bDict->add(typeEntry);
    bDict->add(valEntry);
    
    if(!fieldDict.hasSubDictionary("boundary"))
        fieldDict.add(new iomanagment::Dictionary("boundary"));
    
    fieldDict.subDictionary("boundary").add(bDict);
}

Scalar AvgPressure::avgPressure() const
{
    return m_avgValue;
}

void AvgPressure::updateValues() 
{
    for(size_t e=0; e<boundaryEdges().size(); ++e)
    {
        //get necessery objects and information
        const mesh::RealElement & element= boundaryEdges().element(e);
        
        Vector normal = boundaryEdges().normal(e);
        Scalar nx2 = normal.x()*normal.x();
        Scalar nxny = normal.x()*normal.y();
        Scalar ny2 = normal.y()*normal.y();
        
        const auto & mapL=boundaryEdges()[e].lNodesIdsInElementMatrix();
        
        //Calculate greadient of extrapolated velocity field
        auto localVelocity = m_velocity.slice(element.indexVectorMask());
        
        numArray<Vector> gradU;
        element.grad(xCmps(localVelocity),gradU);

        numArray<Vector> gradV;
        element.grad(yCmps(localVelocity),gradV);
         
        //calculate final values to be set for pressure
        // which is as folows: p=nu*(n*Grad(U)*n-div(U))+B(t)
        (*this)[e] = m_viscosity*
                        ( xCmps(gradU.slice(mapL))*(nx2-1.)  +  
                          yCmps(gradU.slice(mapL))*nxny       + 
                          xCmps(gradV.slice(mapL))*nxny       + 
                          yCmps(gradV.slice(mapL))*(ny2 -1.)
                        ) + m_avgValue;
    }
}
    
void AvgPressure::timeChanged(Scalar time) 
{
    updateValues();
}




    
}//field
}//SEM

