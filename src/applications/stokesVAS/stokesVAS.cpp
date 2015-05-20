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
#include "fieldMath/BIntegralAvgPressure.h"
#include "fieldMath/BIntegral_ddtUdotN.h"
#include "fields/InletOutletVelocity.h"
#include "fieldMath/NonlinearOperatorSplitting.h"
#include "fieldMath/TimeExtrapolation.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using iomanagment::NO_READ;
    using iomanagment::NO_WRITE;
    using iomanagment::READ_IF_PRESENT;
    using iomanagment::AUTO;
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
    
    //--------------------- CREATE FIELDS -----------------------------------
    VectorField U("U", mesh);
    ScalarField p("p", mesh,READ_IF_PRESENT);
    ScalarField p_ext("p_ext",mesh, NO_READ,NO_WRITE);
    
    //--------------------- GET COEFFS -----------------------------------
    std::vector<Scalar> timeCoeffs = selectDiscretizationScheme( Case::solutionControl().timeDiscretization() );
    Scalar dt = Case::time().timeStep();
    
    Scalar Pi = 4.*atan(1.);
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        p_ext = extInTime(p);
        
        las::solve
        (
            a(ddt(U),phi) + a(nu*grad(U),grad(phi)) == a(p_ext,grad(phi))
        );
        
        las::solve
        (
            a(grad(p),grad(phi)) == bInt_rotUdotNxGradPhi(U,nu) - bInt_ddtUdotN(U)
        );
        
        Case::time().fireWriting();
    }
    
    return 0;
}