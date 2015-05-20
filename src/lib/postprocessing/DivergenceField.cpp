
#include "fields/Fields.h"
#include "fieldMath/DivergenceOperator.h"
#include "PostObject.h"

namespace SEM {
    
    class DivergenceCalculator : public PostObject
    {
        field::VectorField m_field;
        field::ScalarField m_divergence;
    public:
        DivergenceCalculator(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
        : m_field(dict.entry("field").value(),mesh,iomanagment::ALWAYS,iomanagment::NO_WRITE), 
          m_divergence("div"+dict.entry("field").value(), mesh, iomanagment::NO_READ,iomanagment::AUTO)
        {
        }
        void calculate()
        {
            m_divergence = div(m_field);
            m_divergence.registryFile()->fireWriting();
        }
    };
    
    REGISTER_IMPLEMENATION(PostObject,DivergenceCalculator,"divergenceField")
    
}

