#include "fields/TCAvgPressure_Stab.h"

#include "fields/GeometricField.h"

#include "utilities/ArrayFunctions.h"

namespace SEM { namespace field {
    
    
    std::string TCAvgPressure_Stab::TYPE_NAME = "timeChangingAvgPressure_stab";
    
    TCAvgPressure_Stab::TCAvgPressure_Stab(const iomanagment::Register& reg, const mesh::Boundary& edges) 
    : TCAvgPressure(reg,edges), m_delta(1e-3), m_Uref(1.)
    {
    }
        
    TCAvgPressure_Stab::~TCAvgPressure_Stab() 
    {
    }

    std::string TCAvgPressure_Stab::typeName() const 
    {
        return TCAvgPressure_Stab::TYPE_NAME;
    }
    
    void TCAvgPressure_Stab::read(const iomanagment::Dictionary& fieldDict) 
    {
        TCAvgPressure::read(fieldDict);
        
        if( fieldDict.subDictionary("boundary").subDictionary(name()).hasEntry("delta") )
            fieldDict.subDictionary("boundary").subDictionary(name()).entry("delta") >> m_delta;
    
        if( fieldDict.subDictionary("boundary").subDictionary(name()).hasEntry("Uref") )
            fieldDict.subDictionary("boundary").subDictionary(name()).entry("Uref") >> m_Uref;
        
    }
    
    void TCAvgPressure_Stab::write(iomanagment::Dictionary& fieldDict) 
    {
        TCAvgPressure::write(fieldDict);
        
        iomanagment::DictEntry *deltaEntry = new iomanagment::DictEntry("delta");
        (*deltaEntry)<<m_delta;
        
        iomanagment::DictEntry *refUEntry = new iomanagment::DictEntry("Uref");
        (*refUEntry)<<m_Uref;
        
        fieldDict.subDictionary("boundary").subDictionary(name()).add(deltaEntry);
        fieldDict.subDictionary("boundary").subDictionary(name()).add(refUEntry);
    }
    
    void TCAvgPressure_Stab::updateValues()
    {
        for(size_t e=0; e<boundaryEdges().size(); ++e)
        {
            //get necessery objects and information
            const mesh::RealElement & element= boundaryEdges().element(e);
            
            Vector normal = boundaryEdges().normal(e);
            Scalar nx2 = normal.x()*normal.x();
            Scalar nxny = normal.x()*normal.y();
            Scalar ny2 = normal.y()*normal.y();
            
            const auto & mapL = boundaryEdges()[e].lNodesIdsInElementMatrix();
            
            //Calculate greadient of extrapolated velocity field
            auto localVelocity = m_velocity.slice(element.indexVectorMask());
            
            numArray<Vector> gradU;
            element.grad(xCmps(localVelocity),gradU);

            numArray<Vector> gradV;
            element.grad(yCmps(localVelocity),gradV);
            
            //calculate final values to be set for pressure
            // which is as folows: p=nu*n*Grad(U)*n+B(t)
            auto edgeGradU = gradU.slice(mapL);
            auto edgeGradV = gradV.slice(mapL);
            auto edgeVel   = localVelocity.slice(mapL);
            
            using array::magSqrt;
            using array::tanh;
            
            numArray<Scalar> stab = 0.25*magSqrt(edgeVel)*( 1.-tanh( dotProd(normal,edgeVel)/(delta()*Uref()) ) );
            
            std::cout<<"STABILIZATION AVG="<<array::avg(stab)<<std::endl;
            
            (*this)[e] = m_viscosity*
                            ( xCmps(edgeGradU)*nx2  
                             +( yCmps(edgeGradU) + xCmps(edgeGradV) )*nxny
                             +yCmps(edgeGradV)*ny2
                            ) 
                        - stab//0.25*magSqrt(edgeVel)*( 1.-tanh( dotProd(normal,edgeVel)/(delta()*Uref()) ) )
                        + m_avgValue;
        }
    }

    
    
    
}//field
}//SEM

