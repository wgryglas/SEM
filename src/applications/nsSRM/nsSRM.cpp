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



void computeNewTimestep
(
    SEM::field::GeometricField<SEM::Vector> & U, 
    const SEM::field::DiscontinousField<SEM::Scalar> & ps, 
    const SEM::field::DiscontinousField<SEM::Scalar> & divU_prev, 
    SEM::Scalar nu, 
    SEM::Scalar dt,
    SEM::Scalar alpha1,
    SEM::Scalar alpha2,
    SEM::Scalar eps
)
{
    using namespace SEM;
    using namespace SEM::field;
    using namespace SEM::mesh;
    
    const Mesh & elements = U.mesh();
    
//     VectorField tmp("tmp", elements, iomanagment::NO_READ, iomanagment::NO_WRITE);
//     
// //    numArray<Scalar> coeffs(4);
// //    coeffs[0] = 1; coeffs[1]=0.5; coeffs[2]=0.5; coeffs[3]=1;
//     numArray<Scalar> coeffs(1);
//     coeffs[0] = 1; 
//     
//     numArray<numArray<Vector> > k(coeffs.size());
//     for(int i=0; i<k.size(); ++i) k[i].resize(U.size(),0.);
//     
//     for(int i=0; i<coeffs.size(); ++i)
//     {    
//         if(i>0)
//             tmp = U + coeffs[i]*k[i-1];
//         else
//             tmp = U;
//         
//         DiscontinousField<Scalar> divU = div(tmp);
//         DiscontinousField<Scalar> forcing = ps - (alpha1/eps)*(divU-divU_prev)/dt - (alpha2/eps)*divU;
//         
//         for( size_t e=0; e < elements.size(); ++e )
//         {
//             numArray2D<Scalar> eStiff=elements[e].stiffMatrix(nu);
//             const numArray<Scalar>& localPs = forcing.element(e);
//             numArray<Vector> localPsGradW;
//             elements[e].weakGradV(localPs,localPsGradW);
//             numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
//             
//             for(int dim=0; dim<2; ++dim)
//             {
//                 auto cmpRhs = CmpTraits<Vector>::cmpArray(k[i],dim);
//                 auto localCmpRhs = cmpRhs.slice(elements[e].indexVectorMask());
//                 
//                 auto cmpU = CmpTraits<Vector>::cmpArray(tmp,dim);
//                 auto localCmpU = cmpU.slice(elements[e].indexVectorMask());
//     
//             
//                 auto cmpLocalPsGradW = CmpTraits<Vector>::cmpArray(localPsGradW,dim);
//                 
//                 //localCmpRhs  = (dt*SEM::array::matrixMul(eStiff,localCmpU) + dt*cmpLocalPsGradW) / massVector;
//                 localCmpRhs  = (-dt*SEM::array::matrixMul(eStiff,localCmpU) ) / massVector;
//             }
//         }
//     }
//     
//     Scalar fact=0;
//     for(int i=0; i<coeffs.size(); ++i) fact += 1./coeffs[i];
//     
//     for(int i=0; i<coeffs.size(); ++i)
//     {
//         U = U + k[i]/(coeffs[i]*fact);
//     }
//     
//     las::applyDirichletBCToSolution(U);
    
     numArray<Vector> rhs(U.size(), 0.);
     numArray<Vector> mass(U.size(), 0.);
     
     DiscontinousField<Scalar> divU = div(U);
     DiscontinousField<Scalar> forcing = ps - (alpha1/eps)*(divU-divU_prev)/dt - (alpha2/eps)*divU;
     //DiscontinousField<Vector> conv  = cDeriv(U,U);

     for( size_t e=0; e < elements.size(); ++e )
     {
        numArray2D<Scalar> eStiff=elements[e].stiffMatrix(nu);
        numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
        numArray<Vector> localPsGradW;
        //const numArray<Vector> & localConv = conv.element(e);
        elements[e].weakGradV(forcing.element(e),localPsGradW);
        
        numArray2D<Scalar> derivMat = elements[e].convDerivMatrix(U.element(e));
        
        for(int dim=0; dim<2; ++dim)
        {
            auto cmpMass = CmpTraits<Vector>::cmpArray(mass,dim);
            auto localCmpMass = cmpMass.slice(elements[e].indexVectorMask());
            
            auto cmpRhs = CmpTraits<Vector>::cmpArray(rhs,dim);
            auto localCmpRhs = cmpRhs.slice(elements[e].indexVectorMask());
            
            auto cmpU = CmpTraits<Vector>::cmpArray(U,dim);
            auto localCmpU = cmpU.slice(elements[e].indexVectorMask());
            
            auto cmpLocalPsGradW = CmpTraits<Vector>::cmpArray(localPsGradW,dim);
            
            //auto cmpLocalConv = CmpTraits<Vector>::cmpArray(localConv,dim);
 
            
            numArray<Scalar> convection(SEM::array::matrixMul(derivMat,localCmpU));
            
           //iomanagment::write(SEM::array::max(SEM::array::abs(convection) ), std::cout) <<std::endl;
            
            
            localCmpRhs  +=  cmpLocalPsGradW - SEM::array::matrixMul(eStiff,localCmpU) - convection * localCmpMass;
            
            localCmpMass += massVector;
        }
     }
     
     U = U.oldField(0) + dt* rhs / mass;
     
     las::applyDirichletBCToSolution(U);
     
     iomanagment::write( SEM::array::max(U) ,std::cout<<"max U:")<<std::endl;
}



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
   
    //--------------------- CREATE FIELDS -----------------------------------
    //Create U field
    VectorField U("U",mesh);
    
     
    //Create p field
    ScalarField p("p",mesh,NO_READ);
    
    VectorField ddtU("ddtU",mesh,NO_READ,NO_WRITE);
    
    ScalarField divergence("divU", mesh, NO_READ);
    
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    Scalar eps = 1e6;
    Scalar alpha1=1;
    Scalar alpha2=1;
    
    las::applyDirichletBCToSolution(U);
    (*const_cast<numArray<Vector>*>(&U.oldField(0))) = U;

    DiscontinousField<Scalar> divU(mesh);
    divU = div(U);
    DiscontinousField<Scalar> divU_old(divU);
    
    DiscontinousField<Scalar> ddtDivU(mesh);
    ddtDivU = 0;

    DiscontinousField<Scalar> disP(mesh);
    
      
    
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
        
        for(int s=0; s< 4; ++s) //selectDiscretizationScheme(U).size()
        {
            
            //DiscontinousField<Vector> convU = dt*cDeriv(U,U);
            //auto rhs =  a(U.oldField(0),phi) + a( (dt*nu)*grad(U), grad(phi) )  - a( dt*disP, grad(phi) ) + a( (dt*alpha1)*ddtDivU+(dt*alpha2)*divU, grad(phi) ); //-  a( convU, phi )
//            auto rhs =  a(eps/dt*U.oldField(0),phi) + a( (eps*nu)*grad(U), grad(phi) )  - a( eps*cDeriv(U,U), phi ) - a( eps*disP, grad(phi) );
	   
            //las::applyDirichletBCToSolution(U);
            
            
            computeNewTimestep(U, disP, divU, nu, dt, alpha1, alpha2, eps);

            
//             las::solve
//             (
//                 a(ddt(U), phi )  == a( -nu*grad(U), grad(phi)) //a( disP -(alpha1/eps*ddtDivU+alpha2/eps*divU), grad(phi) ) //
//             );
            
//             las::solve
//             (
//                     a(ddt(U), phi) == a(-nu*grad(U), grad(phi))
//             );
            
            //las::solve( a( U, phi ) == rhs );
            
            las::applyDirichletBCToSolution(U);
            
            divU = div(U);
            ddtU = ddt(U);
            ddtDivU += div(ddtU);
            
            disP = disP  - ( (alpha1/eps)*divU + (alpha2/eps)*ddtDivU);
        }
        
        divU_old = divU;
        
        p = disP;
	divergence = divU;
        
        Case::time().fireWriting();
        
    }
    
    
    
    return 0;
}