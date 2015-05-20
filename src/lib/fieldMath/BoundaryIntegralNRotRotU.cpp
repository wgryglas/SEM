#include "BoundaryIntegralNRotRotU.h"


namespace SEM {
    
namespace weak {
    
boost::shared_ptr< las::DiscretOperator< Scalar, BoundaryIntegralNRotRotU <AcceptAllBoundary>> > bInt_rotUdotNxGradPhi(const field::GeometricField< Vector >& velocity, Scalar coeff, bool extrap) 
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,BoundaryIntegralNRotRotU<AcceptAllBoundary> > >
    (
        new BoundaryIntegralNRotRotU<AcceptAllBoundary>(velocity,AcceptAllBoundary(), coeff, extrap)
    );
}
boost::shared_ptr< las::DiscretOperator< Scalar, BoundaryIntegralNRotRotU< CheckSingleBoundary > > > bInt_rotUdotNxGradPhi(const field::GeometricField< Vector >& velocity, const string& patchTypeFamily, Scalar coeff, bool extrap) 
{
    return boost::shared_ptr<las::DiscretOperator<Scalar,BoundaryIntegralNRotRotU< CheckSingleBoundary > > >
    (
        new BoundaryIntegralNRotRotU< CheckSingleBoundary >(velocity,CheckSingleBoundary(patchTypeFamily), coeff, extrap )
    );
}

}//weak


}//SEM