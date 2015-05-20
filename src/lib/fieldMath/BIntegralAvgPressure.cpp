#include "BIntegralAvgPressure.h"

namespace SEM {
    
    BIntegralAvgPressure::BIntegralAvgPressure(const field::GeometricField< Scalar >& pressureField,  const field::GeometricField<Vector> & velocityField) 
    : m_pressure(pressureField), m_velocity(velocityField)
    {
    }
    
    
    boost::shared_ptr< SEM::las::DiscretOperator<Vector, BIntegralAvgPressure > > 
    bInt_avgPressure(const field::GeometricField< Scalar >& pressureField,  const field::GeometricField<Vector> & velocityField) 
    {
        return boost::shared_ptr<las::DiscretOperator<Vector,BIntegralAvgPressure> >
        (
            new BIntegralAvgPressure(pressureField, velocityField)
        );
    }
    
}//SEM
