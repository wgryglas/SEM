
#ifndef NONLINEAROPERATORSPLITTING_H
#define NONLINEAROPERATORSPLITTING_H

/// \file NonlinearOperatorSplitting
/// banch of tools wich would allow 
/// explicit evaluation of Navier-Stokes nonlinear
/// term for time dependent PDE 
/// according to Operator-Integrations-Factor-Splitting (OIFS).
///----------------------------------------------------
/// The aim of below tools is to make evaluation of corrected
/// old time results(for time steping method like BDF) which would
/// put explicitly nonlinear term into PDE equation. Such a equation would be now
/// linearized. Old times correction is made by integration special ODE problem
/// with RK4 method.
/// --------------------------------------------------
/// OIFS short description:
/// du/dt = N[u] + L[u]   (1) 
/// whre L-linear operator, N-nonlinear operator(in case of N-S eq. it's convective derivative of velocity)
///
/// Assume that we evaluate time derivative by some BDF sheme:
/// du/dt = 1/dt*( g1*u(t+dt) + g2*u(t) + g3*u(t-dt) + ... ) (2)
///
/// then equation (1) can be now rewritten to below form:
/// g1/dt*u(t+dt) = N[u] + L[u] - 1./dt(g2*u(t) + g3*u(t-dt) + ...) (3)
///
/// When we apply OIFS idea to (3) (we wont describe derivation of this idea, easy to find in internet)
/// then we get linearized equation (3) to below form:
/// 
/// g1/dt*u(t+dt) = L[u] - 1./dt(g2*v(t) + g3*v(t-dt) + ...) (4)
/// 
/// whre v is special value of u(t) evaluated by solving below 
/// ODE problem:
/// dv(h)/dh = N[v(h)]
/// v(h=0) = u(t*), h = {t*,t+dt}
/// t*- depends on which "v" term from (4) we calculate, eg t*=t, t*=t-dt, ....
///
/// Generaly speaking : if we want to replace nonlinear transient PDE with
/// linearized one, then we modify time derivative discret opeartor by substituting
/// solution old times with new values. Then solve new PDE equation, 
/// which no longer consists of nonlinear term.
/// New values of old times are obtained by integrating
/// ODE which tooks as RHS function nonlinear term of PDE.
/// ---------------------------------------------------------------


#include "fields/GeometricField.h"
#include "mesh/Mesh.h"
#include "utilities/ArrayFunctions.h"
#include "utilities/RK4Integrator.h"
#include "fields/DiscontinousField.h"

namespace SEM { 

 
field::DiscontinousField<Vector> convDerivativeEvaluation(Scalar time,const field::DiscontinousField<Vector>& f);

field::DiscontinousField<Vector> oifsConvDerivative(const field::GeometricField<Vector> & f, Scalar timeStep, unsigned int prevTime);



}//SEM




#endif // NONLINEAROPERATORSPLITTING_H
