#ifndef WALLPRESSUREPATCH_H
#define WALLPRESSUREPATCH_H

#include "PatchField.h"
#include "components/Basic.h"
#include "iomanagment/Register.h"
#include "utilities/numArray.h"

namespace SEM{ namespace field {

template<typename T>    
class GeometricField;

extern const std::string WALL_PRESSURE;
    
class WallPressurePatch: public PatchField<Scalar>
{
    typedef GeometricField<Vector> VectorField;
    typedef PatchField<Scalar> baseType;

    WallPressurePatch& operator=(const WallPressurePatch& other);
    bool operator==(const WallPressurePatch& other);
    
//     const VectorField & m_velField;
//     const VectorField & m_residueField;
    
public:
    static std::string TYPE_NAME;
    
    using baseType::boundaryEdges;
    
//     WallPressurePatch(const VectorField & velocity,const VectorField & residue, const mesh::Boundary & edges);
    
//     WallPressurePatch(const iomanagment::Register & reg,const mesh::Boundary & edges);
    
    WallPressurePatch(const mesh::Boundary & edges);
    
    //auto copy constructor 
    
    ~WallPressurePatch();

    std::string typeName() const {return WallPressurePatch::TYPE_NAME;}
    std::string typeFamily() const {return WALL_PRESSURE;}
    
    //IO
    void read(const iomanagment::Dictionary &fieldDict);
    
    void write(iomanagment::Dictionary &fieldDict);
    
    //assigment handling (in fact not necessery)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE( =, Scalar, WallPressurePatch,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(+=, Scalar, WallPressurePatch,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(-=, Scalar, WallPressurePatch,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(*=, Scalar, WallPressurePatch,baseType)
    ALLOW_EXP_TEMP_UNARY_OP_IN_DERIVED_FROM_NUMARRAYBASE(/=, Scalar, WallPressurePatch,baseType)
};

}
}

#endif // WALLPRESSUREPATCH_H
