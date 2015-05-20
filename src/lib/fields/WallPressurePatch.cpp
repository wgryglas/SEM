#include "WallPressurePatch.h"

#include "iomanagment/case.h"
#include "components/Basic.h"
#include "utilities/ArrayFunctions.h"

namespace SEM  { namespace field {
    
const std::string WALL_PRESSURE = "wallPressure";    
std::string WallPressurePatch::TYPE_NAME = "wallPressure";
    
// WallPressurePatch::WallPressurePatch(const VectorField& velocity, const VectorField & residue, const mesh::Boundary& edges) 
// : baseType(edges), m_velField(velocity), m_residueField(residue)
// {
// }
// 
// WallPressurePatch::WallPressurePatch(const iomanagment::Register& reg, const mesh::Boundary& edges) 
// : baseType(edges), m_velField(*reg.object<VectorField>("U")), m_residueField(*reg.object<VectorField>("Residue"))
// {
// }

WallPressurePatch::WallPressurePatch(const mesh::Boundary& edges) 
: baseType(edges,0,0)
{
}


WallPressurePatch::~WallPressurePatch() 
{
}

void WallPressurePatch::read(const iomanagment::Dictionary& fieldDict) 
{
    //do nothing
}

void WallPressurePatch::write(iomanagment::Dictionary& fieldDict) 
{
    //Write type only, values are uninportant
    iomanagment::DictEntry *typeEntry =new iomanagment::DictEntry("type",typeName());
    
    iomanagment::Dictionary *bDict = new iomanagment::Dictionary(name());
    bDict->add(typeEntry);
    
    if(!fieldDict.hasSubDictionary("boundary"))
        fieldDict.add(new iomanagment::Dictionary("boundary"));
    
    fieldDict.subDictionary("boundary").add(bDict);
}

// numArray< Scalar > WallPressurePatch::weakForm(size_t edge) const 
// {
//     if(m_velField.patch(name()).dirichletCoeff()!=1.0)
//     {
//         ErrorInFunction<<"velocity patch "<<name()<<" is not consistient with "<<typeName()<<std::endl
//                        <<typeName()<<" requires fixedValue patch on velocity"<<iomanagment::endProgram;
//     }
//     
//     Scalar nu                               = Case::material()["nu"];
//     std::vector<Scalar> timeCoeffs          = selectDiscretizationScheme(Case::solutionControl().timeDiscretization());
//     Scalar dt                               = Case::time().timeStep();
//     const mesh::BoundaryEdge & bEdge        = boundaryEdges()[edge];
//     const numArray<Scalar> & integralCoeffs = bEdge.neumanBCPreValue();
//     Vector edgeVelocity                     = m_velField.patch(name())[edge];
//     Vector edgeNormal                       = boundaryEdges().normal(edge);
//     auto edgeResidue                        = m_residueField.slice(bEdge.gNodesIds());   
//     
//     
//     //Calculate extrapolated velocity in element
//     const std::vector<int> & elementMap = boundaryEdges().element(edge).indexVectorMask();
//     numArray<Vector> elementExtrapVelocity = 3.*(m_velField.cachedField(0).slice(elementMap) - m_velField.cachedField(1).slice(elementMap) ) + m_velField.cachedField(2).slice(elementMap);
//     
//     //Assign values for all nodes in edge element
//     numArray<Scalar> value = -nu*boundaryEdges().element(edge).boundInt_NxRotRotU(elementExtrapVelocity, bEdge.edgeId() );
//     
//     //Assign values for nodes coresponding to edge
//     value.slice(bEdge.lNodesIdsInElementMatrix()) -= dotProd(edgeNormal, timeCoeffs[0]/dt*edgeVelocity - edgeResidue)*integralCoeffs;
//     
//     return value;
// }



}//field
}//SEM