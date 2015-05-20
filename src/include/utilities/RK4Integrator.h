#ifndef RK4INTEGRATOR_H
#define RK4INTEGRATOR_H

#include <functional>
#include "components/Basic.h"

/// \class RK4Integrator 
/// Integrator for ordinary differential equation
/// It's template class, where argument is type of 
/// entity that has to be integrated. IType can be
/// any type only if this type has
/// defined operator "+", multiplication with 
/// scalar value and copy constructor.
/// -------------------------------------------------------
/// This class would work properly with scalar type
/// or most SEM defined lists (eg. numArray, DiscontinousField etc.)

namespace SEM {

template<typename IType>
struct RK4Integrator
{
    typedef std::function<IType(Scalar, const IType &)> RHSFunction;    
public:
    RK4Integrator(std::function<IType(Scalar, const IType &)> rhs)
    : m_rhsFunction(rhs)
    {
    }
    
    /// \brief integrate 
    /// function integrates rhs function for independent variable from start to end with
    /// solution at start equal to y0. As result returns solution value 
    /// at end.
    /// \var start -  independent variable begin value
    /// \var end -  independent variable end value
    /// \var numberOfSteps - number of integration steps to be performed
    /// \var y0 - solution for independent variable start
    /// \return value integrated at time equal to end
    IType integrate(Scalar start, Scalar step, size_t numberOfSteps, const IType & y0) const
    {
        IType y = y0;
        for(size_t t=0; t<numberOfSteps; ++t)
        {
            y = integrateStep(start,step,y);
            start +=step;
        }
        return y;
    }
    
    /// \brief integrate 
    /// function integrates rhs function for independent variable from start to end with
    /// solution at start equal to Y[0]. As result resutrns
    /// resized list Y filled with values for each integration
    /// step.
    /// \var start -  independent variable begin value
    /// \var end -  independent variable end value
    /// \var numberOfSteps - number of solution points
    /// \var Y - solution storage. Y[0] is assumed to be initial condition.
    /// \return all integration values list.
    void integrate(Scalar start, Scalar end, size_t numberOfSteps, std::vector<IType>& Y) const
    {
        Scalar step = (end-start)/(numberOfSteps-1);
        Y.resize(numberOfSteps);
        
        for(size_t t=0; t<numberOfSteps; ++t)
        {
            Y[t+1] = integrateStep(start, step, Y[t]);
            start+=step;
        }
    }
private:
    RHSFunction m_rhsFunction;
    
    IType integrateStep(const Scalar &time, const Scalar &step, const IType & y0) const
    {
        IType k1 = m_rhsFunction(time,y0);
        IType k2 = m_rhsFunction(time+step/2,y0+step/2*k1);
        IType k3 = m_rhsFunction(time+step/2,y0+step/2*k2);
        IType k4 = m_rhsFunction(time+step,y0+step*k3);
        
        return y0 + step/6 * ( k1 + 2.*k2 + 2.*k3 + k4 );
    }
};

}//SEM

#endif // RK4INTEGRATOR_H
