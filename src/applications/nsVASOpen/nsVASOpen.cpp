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
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using SEM::iomanagment::NO_READ;
    using SEM::iomanagment::NO_WRITE;
    using SEM::weak::phi;
    using SEM::weak::ddt;
    using SEM::field::grad;
    using SEM::weak::ddn;
    using SEM::weak::grad;
    using SEM::weak::a;
    using SEM::weak::bInt_rotUdotNxGradPhi;
    using SEM::weak::bInt_UdotNPhi;
    using SEM::bInt_avgPressure;
    using SEM::bInt_ddtUdotN;
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
    
//     //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
    
    ScalarField p("p",mesh, READ_IF_PRESENT);
    
    //ScalarField p_ext("p_ext",mesh,NO_READ,NO_WRITE);
    
    //--------------------- GET COEFFS -----------------------------------
    
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    int nCorrectors = 5;
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        std::vector<Scalar> timeCoeffs = selectDiscretizationScheme(U);
        
        DiscontinousField<Vector> oldTimeTerms(mesh);
        oldTimeTerms = 0.;
        for(unsigned int i=1; i<timeCoeffs.size();++i)
        {
            oldTimeTerms += timeCoeffs[i]/dt*oifsConvDerivative(U,dt,i-1);
        }
        
        //p=extInTime(p);
        
        for(int i=0; i<nCorrectors; ++i)
        {
            las::solve
            (
                a(timeCoeffs[0]/dt,U,phi) + a(nu*grad(U),grad(phi)) == a(p,grad(phi)) - a(oldTimeTerms,phi) - bInt_avgPressure(p,U)
            );
            
            p.updateBoundaryConditions();
            
            las::solve
            (
                a(grad(p),grad(phi)) == a(-cDeriv(U,U),grad(phi)) + bInt_rotUdotNxGradPhi(U,WALL_PRESSURE,nu) - bInt_ddtUdotN(U)
            );
        }
        
        if(nCorrectors != 1)
            --nCorrectors;
        
        Case::time().fireWriting();
        
    }
    
    return 0;
}