#ifndef _ddt_H_
#define _ddt_H_


#include <vector>

#include <boost/smart_ptr/make_shared.hpp>

#include "fields/GeometricField.h"
#include "iomanagment/InfoStream.h"
#include "iomanagment/case.h"
#include "solver/DiscretOperator.h"
#include "solver/DiscretEquation.h"
#include "TimeDerivativeDiscreteOperator.h"

namespace SEM {
    
/** *************************************************************************
 * \brief selectDiscretizationScheme
 * prepares apropriate coefficients for selected discretization scheme.
 ****************************************************************************/
inline std::vector<double> selectDiscretizationScheme(TimeDiscretization scheme)
{
    std::vector<double> coeffs;
    switch(scheme)
    {
        case Euler:
            coeffs.resize(2);
            coeffs[0]=1.0;
            coeffs[1]=-1.0;
            break;
        case BDF2:
            coeffs.resize(3);
            coeffs[0]=1.5;
            coeffs[1]=-2.0;
            coeffs[2]=0.5;
            break;
        default:
            break;
    }
    
    return coeffs;
}

template<typename T> 
inline std::vector<double> selectDiscretizationScheme(const field::GeometricField<T> & f,TimeDiscretization scheme = Case::solutionControl().timeDiscretization() )
{
    std::vector<double> coeffs;
    
    if(scheme > f.oldFieldsNumber())
    {
        scheme = Euler;
    }
    
    switch(scheme)
    {
        case Euler:
            coeffs.resize(2);
            coeffs[0]=1.0;
            coeffs[1]=-1.0;
            break;
        case BDF2:
            coeffs.resize(3);
            coeffs[0]=1.5;
            coeffs[1]=-2.0;
            coeffs[2]=0.5;
            break;
        default:
            break;
    }
    
    return coeffs;
}


/** *************************************************************************
 * \brief ddt 
 *  Function which prepares TimeFirstDerivativeBuilder::equationPair from field, and 
 *  discretization scheme, time step are read from file indirectly-by 
 *  SolutionControl class object stored in Case as singleton.
 *  TimeFirstDerivativeBuilder::equationPair--->builder connected to field
 * \return object which is able to apply this part of equation to metrix and
 *         rhs vector
 ****************************************************************************/
// template<typename T>
// typename TimeDerivativeDiscreteOperator<T>::ref ddt(field::GeometricField<T> &f)
// {
//     return typename TimeDerivativeDiscreteOperator<T>::ref
//     (
//         new TimeDerivativeDiscreteOperator<T>(f,Case::time().timeStep(), selectDiscretizationScheme(Case::solutionControl().timeDiscretization() ) )
//     );
// }


/** *************************************************************************
 * \brief ddt 
 *  Function which prepares TimeFirstDerivativeBuilder::equationPair from field, time step 
 *  and discretization scheme.
 * TimeFirstDerivativeBuilder::equationPair--->builder connected to field
 * \param f filed which shall be discretized
 * \param dt time step
 * \param scheme time discretization scheme (currently BDF2 | Euler)
 * \return object which is able to apply this part of equation to metrix and
 *         rhs vector
 ****************************************************************************/
// template<typename T>
// boost::shared_ptr<TimeDerivativeDiscreteOperator<T> > ddt(field::GeometricField<T> &f,double dt, TimeDiscretization scheme)
// {
//     return typename TimeDerivativeDiscreteOperator<T>::ref
//     (
//         new TimeDerivativeDiscreteOperator<T>(f, dt, selectDiscretizationScheme(scheme) )
//     );
// }

}//SEM



#endif //_ddt_H_



