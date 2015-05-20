#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "fieldMath/f.h"
#include "fieldMath/postCalculation.h"
#include "fieldMath/weakGradV.h"

#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "utilities/VectorUtils.h"
#include "fields/DiscontinousField.h"
#include "fieldMath/WeakFormDefinition.h"
#include "fieldMath/DivergenceOperator.h"
#include "fieldMath/GradientOperator.h"
#include "fieldMath/weakDivW.h"
#include "fieldMath/BoundaryIntegralNRotRotU.h"
#include "fieldMath/ConvectiveDerivative.h"
#include "fieldMath/NonlinearOperatorSplitting.h"
#include "fieldMath/BIntegral_ddtUdotN.h"
#include "fieldMath/BIntegralAvgPressure.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using SEM::weak::phi;
    using SEM::weak::ddt;
    using SEM::field::grad;
    using SEM::weak::ddn;
    using SEM::weak::grad;
    using SEM::weak::a;
    using SEM::weak::bInt_rotUdotNxGradPhi;
    using SEM::weak::bInt_UdotNPhi;
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
    //     //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
     
    //Create p field
    ScalarField p("p",mesh,NO_READ);
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    Scalar eps = 1e-2;
    Scalar alpha1=1./eps;
    Scalar alpha2=1./eps;

    DiscontinousField<Scalar> divU(mesh);
    divU = div(U);
    
    DiscontinousField<Scalar> ddtDivU(mesh);
    ddtDivU = 0;

    DiscontinousField<Scalar> disP(mesh);
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
	//Fixed euler
        //std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(U);
        
        for(int s=0; s<2; ++s)
        {
            auto rhs =  a(U.oldField(0),phi) + a( (dt*nu)*grad(U), grad(phi) ) -  a( dt*cDeriv(U,U), phi ) - a( dt*disP, grad(phi) ) - a( (dt*alpha1)*ddtDivU+(dt*alpha2)*divU, grad(phi) );
//            auto rhs =  a(eps/dt*U.oldField(0),phi) + a( (eps*nu)*grad(U), grad(phi) )  - a( eps*cDeriv(U,U), phi ) - a( eps*disP, grad(phi) );
	   
            las::solve( a( U, phi ) == rhs );
            
	    ddtDivU = -divU/dt;
            divU = div(U);
            ddtDivU += divU/dt;
            
            disP = disP  - (alpha1*divU + alpha2*ddtDivU);
        }
        
        p = disP;
	
        Case::time().fireWriting();
    }
    
    return 0;
}