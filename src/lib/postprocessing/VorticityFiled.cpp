
#include "fields/Fields.h"
#include "fieldMath/RotOperator.h"
#include "PostObject.h"

namespace SEM {
    
    class vorticityCalculator : public PostObject
    {
        field::VectorField m_field;
        field::ScalarField m_vorticity;
    public:
        vorticityCalculator(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
        : m_field(dict.entry("field").value(),mesh,iomanagment::ALWAYS,iomanagment::NO_WRITE), 
        m_vorticity("curl"+dict.entry("field").value(), mesh, iomanagment::NO_READ,iomanagment::AUTO)
        {
        }
        void calculate()
        {
            std::cout << m_field.mesh().size() <<std::endl;
            m_vorticity = rot(m_field);
            m_vorticity.registryFile()->fireWriting();
        }
    };
    
    REGISTER_IMPLEMENATION(PostObject,vorticityCalculator,"vorticityField")
    
}

