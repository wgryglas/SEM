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
#include "fieldMath/BIntegral_ddtUdotN.h"

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
    ScalarField p("p",mesh, READ_IF_PRESENT);
    ScalarField p_real("p_real",mesh,NO_READ);
    VectorField Residue("Residue",mesh,NO_READ,NO_WRITE);
    
    ScalarField divU("divU",mesh,NO_READ);
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    
    //initial conditions:
    p_real = p;
    p_real.registryFile()->fireWriting();
    
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        //pressure equation
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( U );
        Residue = 0.;
        for(unsigned int i=1; i<timeCoeffs.size(); ++i)
        {
            Residue-=timeCoeffs[i]/dt*U.oldField(i-1);
        }
        las::solve
        (
            a(grad(p),grad(phi)) == a(Residue,grad(phi))  + bInt_rotUdotNxGradPhi(U,nu,true) - bInt_UdotNPhi(U,timeCoeffs[0]/dt)
        );
        
        // Velocity equation
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi)) == a(p,grad(phi))
        );
        
        p_real = p - nu*div(U);
        divU = div(U);
        
        Case::time().fireWriting();
    }
    
    return 0;
}