#ifndef PATCHSELECTION_H
#define PATCHSELECTION_H

#include <map>
#include <string>
#include <functional>

#include "PatchField.h"
#include "iomanagment/Register.h"
#include "mesh/Mesh.h"

#include "fields/FixedValuePatch.h"
#include "fields/FixedGradientPatch.h"
#include "fields/WallPressurePatch.h"
#include "fields/AvgPressure.h"
#include "fields/InletOutletVelocity.h"
#include "fields/TimeChangingValuePatch.h"
#include "fields/TimeChangingGradientPatch.h"
#include "fields/TCAvgPressure.h"
#include "fields/TCAvgPressure_Stab.h"
#include "fields/ExpressionValuePatch.h"
#include "fields/ExpressionGradientPatch.h"
#include "fields/InterpolatedFixedValuePatch.h"


namespace SEM { namespace field {

    
template<typename T, class C>
PatchField<T>* allArgs(const iomanagment::Register &reg, const mesh::Boundary &edges) 
{
    return new C(reg, edges);
}
template<typename T, class C>
PatchField<T>* edgesOnly(const iomanagment::Register &reg, const mesh::Boundary &edges) 
{
    return new C(edges);
}

    
template<typename T>
struct PatchFactory
{
    typedef std::string Key;
    typedef std::function<PatchField<T>*(const iomanagment::Register&,const mesh::Boundary&)> Creator;
    
    typedef std::map<Key,Creator> Selectors;
    
    PatchFactory();
    PatchField<T>* select(const std::string &patch, const iomanagment::Register &reg, const mesh::Boundary &edges) const;
private:
    Selectors m_selectors;
};


template<typename T>
struct TypeSpecificFactory
{
    static void configure(typename PatchFactory<T>::Selectors & selectors)
    {
        //filled by specialization
    }
};

template<>
struct TypeSpecificFactory<Scalar>
{
    static void configure(typename PatchFactory<Scalar>::Selectors & selectors)
    {
        selectors[WallPressurePatch::TYPE_NAME] = edgesOnly<Scalar,WallPressurePatch>;
        selectors[AvgPressure::TYPE_NAME] = allArgs<Scalar,AvgPressure>;
        selectors[TCAvgPressure::TYPE_NAME] = allArgs<Scalar,TCAvgPressure>;
        selectors[TCAvgPressure_Stab::TYPE_NAME] = allArgs<Scalar,TCAvgPressure_Stab>;
    }
};

template<>
struct TypeSpecificFactory<Vector>
{
    static void configure(typename PatchFactory<Vector>::Selectors & selectors)
    {
        selectors[InletOutletVelocity::TYPE_NAME] = edgesOnly<Vector,InletOutletVelocity>;
    }
};


template<typename T>
PatchFactory<T>::PatchFactory()
{
    m_selectors[FixedValuePatch<T>::TYPE_NAME]=edgesOnly<T,FixedValuePatch<T> >;
    m_selectors[FixedGradientPatch<T>::TYPE_NAME]=edgesOnly<T,FixedGradientPatch<T> >;
    m_selectors[TimeChangingValuePatch<T>::TYPE_NAME]=edgesOnly<T,TimeChangingValuePatch<T> >;
    m_selectors[TimeChangingGradientPatch<T>::TYPE_NAME]=edgesOnly<T,TimeChangingGradientPatch<T> >;
    m_selectors[ExpressionValuePatch<T>::TYPE_NAME]=edgesOnly<T,ExpressionValuePatch<T> >;
    m_selectors[ExpressionGradientPatch<T>::TYPE_NAME]=edgesOnly<T,ExpressionGradientPatch<T> >;
    m_selectors[InterpolatedFixedValuePatch<T>::TYPE_NAME]=edgesOnly<T,InterpolatedFixedValuePatch<T> >;
    
    TypeSpecificFactory<T>::configure(m_selectors);
}

template<typename T>
PatchField<T>* PatchFactory<T>::select(const std::string &patch, const iomanagment::Register &reg, const mesh::Boundary &edges) const
{
    auto creator = m_selectors.find(patch);
    
    if(creator != m_selectors.end())
    {
        return creator->second(reg,edges);
    }
    else
    {
        using namespace SEM::iomanagment;
        ErrorInFunction<<"Selected patch "<< patch << " is not known to program"<<endProgram;
        return nullptr;
    }
}


}
}

#endif // PATCHSELECTION_H
