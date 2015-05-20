#include "NonlinearOperatorSplitting.h"

namespace SEM {
    
field::DiscontinousField< Vector > convDerivativeEvaluation(Scalar time, const field::DiscontinousField< Vector >& f) 
{
    using array::matrixMul;
    using field::DiscontinousField;
    
    DiscontinousField<Vector> deriv(f.mesh());
    for(size_t e=0; e<f.elementsNumber(); ++e)
    {
        //compute convective operator matrix
        numArray2D<Scalar> derivMatrix = deriv.mesh()[e].convDerivMatrix(f[e]);

        //compute conv. derivative of velocity field
        xCmps(deriv[e]) = -matrixMul( derivMatrix, xCmps(f[e]));
        yCmps(deriv[e]) = -matrixMul( derivMatrix, yCmps(f[e]));
    }
    return deriv;
}

field::DiscontinousField<Vector> oifsConvDerivative(const field::GeometricField< Vector >& f, Scalar timeStep, unsigned int prevTime) 
{
    RK4Integrator<field::DiscontinousField<Vector> >  integrator(convDerivativeEvaluation);
    unsigned int division =5;
    size_t nStep = (prevTime+1)*division;
    return integrator.integrate(0.,timeStep/division, nStep, field::DiscontinousField<Vector>(f.oldField(prevTime),f.mesh()) );
}

}//SEM