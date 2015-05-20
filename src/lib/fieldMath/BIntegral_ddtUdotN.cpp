
#include "BIntegral_ddtUdotN.h"
namespace SEM {
    
    bool AcceptAllBoundary::operator()(const std::string& toCheck) const 
    {
        SEM_UNUSED(toCheck) return true;
    }

    CheckSingleBoundary::CheckSingleBoundary(const std::string& family) 
    : m_family(family)
    {
    }

    bool CheckSingleBoundary::operator()(const std::string& toCheck) const 
    {
        return m_family == toCheck;
    }
    
    boost::shared_ptr< SEM::las::DiscretOperator< SEM::Scalar, SEM::BIntegral_ddtUdotN< SEM::CheckSingleBoundary > > > bInt_ddtUdotN(const SEM::field::GeometricField< SEM::Vector >& field, const std::string& boundaryFamily) 
    {
        return boost::shared_ptr<las::DiscretOperator<Scalar,BIntegral_ddtUdotN<CheckSingleBoundary> > >
        (
            new BIntegral_ddtUdotN<CheckSingleBoundary>(field,CheckSingleBoundary(boundaryFamily))
        );
    }
    
    boost::shared_ptr< SEM::las::DiscretOperator< SEM::Scalar, SEM::BIntegral_ddtUdotN< SEM::AcceptAllBoundary > > > bInt_ddtUdotN(const SEM::field::GeometricField< SEM::Vector >& field) 
    {
        return boost::shared_ptr<las::DiscretOperator<Scalar,BIntegral_ddtUdotN<AcceptAllBoundary> > >
        (
            new BIntegral_ddtUdotN<AcceptAllBoundary>(field,AcceptAllBoundary())
        );
    }
    
    
}//SEM

