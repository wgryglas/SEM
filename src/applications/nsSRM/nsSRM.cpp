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
    Scalar alpha1=1;
    Scalar alpha2=1;

    DiscontinousField<Scalar> divU(mesh);
    divU = div(U);
    
    DiscontinousField<Scalar> ddtDivU(mesh);
    ddtDivU = 0;
       
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(U);
        
        for(int i=0; i<timeCoeffs.size(); ++i)
        {
            DiscontinousField<Vector> Residue = a(eps*nu*grad(U),grad(phi));
            Residue -= a(cDeriv(eps*U,U),phi);
            Residue -= a(eps*p, grad(phi));
            Residue += a(divU+ddtDivU,grad(phi));
            
            for(unsigned int i=1; i<timeCoeffs.size(); ++i)
            {
                Residue-=eps*timeCoeffs[i]/dt*U.oldField(i-1);
            }
            
            U = U + dt/timeCoeffs[0]* Residue;
            
            divU = div(U);
            ddtDivU = ddt(divU);
            
            p = p  -  1./eps * (alpha1*divU + alpha2*ddtDivU);
        }
        
        Case::time().fireWriting();
    }
    
    return 0;
}