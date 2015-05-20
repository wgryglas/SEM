/*
 * 
 */

#ifndef INLETOUTLETVELOCITY_H
#define INLETOUTLETVELOCITY_H

#include <string>

#include "components/Basic.h"
#include "PatchField.h"

namespace SEM { namespace field { 


extern const std::string VELOCITY_DO_NOTHING;
    
class InletOutletVelocity : public PatchField<Vector>
{

    typedef PatchField<Vector> baseType;
public:
    static std::string TYPE_NAME;
    
    InletOutletVelocity(const mesh::Boundary & edges);
    
    //auto copy constructor 
    
    ~InletOutletVelocity();
    
    std::string typeName() const;
    
    std::string typeFamily() const;
    
    //IO
    void read(const iomanagment::Dictionary &fieldDict);
    
    void write(iomanagment::Dictionary &fieldDict);
    
    //assigment handling (in fact not necessery)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, Vector, InletOutletVelocity,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, Vector, InletOutletVelocity,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, Vector, InletOutletVelocity,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, Vector, InletOutletVelocity,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, Vector, InletOutletVelocity,baseType)
};

}//field
}//SEM

#endif // INLETOUTLETVELOCITY_H
