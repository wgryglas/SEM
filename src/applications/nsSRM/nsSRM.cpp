#include "sstream"

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
    //const SEM::field::DiscontinousField<SEM::Scalar> & divU_prev, 
    SEM::Scalar nu, 
    SEM::Scalar dt,
    SEM::Scalar alpha1,
    SEM::Scalar alpha2,
    SEM::Scalar eps,
    bool predictTime = false
)
{
    using namespace SEM;
    using namespace SEM::field;
    using namespace SEM::mesh;
    
    const Mesh & elements = U.mesh();

// ---------------------------- rk4 --------------------------------------------
/*
     VectorField tmp("tmp", elements, iomanagment::NO_READ, iomanagment::NO_WRITE);
     numArray<Vector> mass(U.size(), 0.);
    numArray<Scalar> coeffs(4);
    coeffs[0] = 1; coeffs[1]=0.5; coeffs[2]=0.5; coeffs[3]=1;
//     numArray<Scalar> coeffs(1);
//     coeffs[0] = 1; 
     
     numArray< numArray<Vector> > k(coeffs.size());
     for(int i=0; i<k.size(); ++i) k[i].resize(U.size(),0.);
     
     for(int i=0; i<coeffs.size(); ++i)
     {    
         if(i>0)
             tmp = U + coeffs[i] * dt * k[i-1] / mass;
         else
             tmp = U;
         
         DiscontinousField<Scalar> divU = div(tmp);
         DiscontinousField<Scalar> forcing = ps - (alpha1/eps)*(divU-divU_prev)/dt - (alpha2/eps)*divU;

         for( size_t e=0; e < elements.size(); ++e )
         {
	     numArray2D<Scalar> eStiff=elements[e].stiffMatrix(nu);
             numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
             numArray<Vector> localPsGradW;
             elements[e].weakGradV(forcing.element(e),localPsGradW);
             numArray2D<Scalar> derivMat = elements[e].convDerivMatrix(U.element(e));

	     for(int dim=0; dim<2; ++dim)
	     {
		  auto cmpMass = CmpTraits<Vector>::cmpArray(mass,dim);
		  auto localCmpMass = cmpMass.slice(elements[e].indexVectorMask());
		    
		  auto cmpRhs = CmpTraits<Vector>::cmpArray(k[i],dim);
		  auto localCmpRhs = cmpRhs.slice(elements[e].indexVectorMask());
		    
		  auto cmpU = CmpTraits<Vector>::cmpArray(tmp,dim);
		  auto localCmpU = cmpU.slice(elements[e].indexVectorMask());
		    
		  auto cmpLocalPsGradW = CmpTraits<Vector>::cmpArray(localPsGradW,dim);
		    
		  numArray<Scalar> convection(SEM::array::matrixMul(derivMat,localCmpU));
		    
		  localCmpRhs  +=  cmpLocalPsGradW - SEM::array::matrixMul(eStiff,localCmpU) - convection * massVector;
		  
                  if(i == 0) localCmpMass += massVector;
	     }
	}
     }
 
     Scalar fact=0;
     for(int i=0; i<coeffs.size(); ++i) fact += 1./coeffs[i];
     
     U = U.oldField(0);
     for(int i=0; i<coeffs.size(); ++i)
     {
         U += ( dt/(coeffs[i]*fact) ) * k[i] / mass;
     }
     
*/
// ------------------------------------ euler -----------------------------------

     numArray<Vector> rhs(U.size(), 0.);
     numArray<Vector> mass(U.size(), 0.);
     
     //DiscontinousField<Scalar> divU = div(U);
     static VectorField ddtU("ddtU", U.mesh(), SEM::iomanagment::NO_READ, SEM::iomanagment::NO_WRITE);
     
     if(predictTime)
        ddtU = (U.oldField(0)-U.oldField(1) )/dt;
     else
        ddtU = (U - U.oldField(0)) / dt;
     
     DiscontinousField<Scalar> forcing = ps - (alpha1/eps)*div(ddtU) - (alpha2/eps)*div(U);
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
            
            numArray<Scalar> convection(SEM::array::matrixMul(derivMat,localCmpU));
            
            localCmpRhs  +=  cmpLocalPsGradW - SEM::array::matrixMul(eStiff,localCmpU) - convection * massVector;
            
            localCmpMass += massVector;
        }
     }
     
     U = U.oldField(0) + dt* rhs / mass;
     
//------------------ adams-bashford mehtod --------------------------------------
/*    
    numArray< numArray<Vector> > rhs(2);
    for(int i=0; i<2; ++i) rhs[i].resize(U.size(),0.);

    numArray<Vector> mass(U.size(), 0.);
    
    DiscontinousField<Scalar> divU = div(U);
    DiscontinousField<Scalar> forcing = ps;
    for(int i=0; i<2; ++i)
    {        
        if(i == 0)
           forcing -= (alpha1/eps)*(divU-divU_prev)/dt + (alpha2/eps)*divU;
        else
           forcing = ps - (alpha1/eps)*(div(U) - divU)/dt + (alpha2/eps)*div(U);

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
                
                auto cmpRhs = CmpTraits<Vector>::cmpArray(rhs[i],dim);
                auto localCmpRhs = cmpRhs.slice(elements[e].indexVectorMask());
                
                auto cmpU = CmpTraits<Vector>::cmpArray(U,dim);
                auto localCmpU = cmpU.slice(elements[e].indexVectorMask());
                
                auto cmpLocalPsGradW = CmpTraits<Vector>::cmpArray(localPsGradW,dim);
                
                numArray<Scalar> convection(SEM::array::matrixMul(derivMat,localCmpU));
                
                localCmpRhs  +=  cmpLocalPsGradW - SEM::array::matrixMul(eStiff,localCmpU) - convection * massVector;
                
                if(i==0) localCmpMass += massVector;
            }
        }

        if(i==0) U = U.oldField(0) + dt* rhs[0] / mass;VectorField ddtU("ddtU", U.mesh(), SEM::iomanagment::NO_READ, SEM::iomanagment::NO_WRITE);
        else     U = U.oldField(0) + dt* 0.5*(rhs[0] + rhs[1]) / mass;
    }
*/
    las::applyDirichletBCToSolution(U);

    iomanagment::write( SEM::array::max(U) ,std::cout<<"max U:")<<std::endl;
}

/*
void projectVelocity(SEM::field::VectorField &U)
{
    using namespace SEM;
    using namespace SEM::field;
    using namespace SEM::mesh;


    const Mesh & elements = U.mesh();
    
    double alpha = 1;
    
    DiscontinousField<Scalar> q(elements);
    q = 0.;
    numArray<Vector> rhs(U.size(), 0.);
    numArray<Vector> mass(U.size(), 0.);
    
    for(int i=0; i<5; ++i)
    {
        DiscontinousField<Scalar> forcing = alpha*div(U) - q; 
        
        for( size_t e=0; e < elements.size(); ++e )
        {
            numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
            numArray<Vector> localPsGradW;
            elements[e].weakGradV(forcing.element(e),localPsGradW);
            
            for(int dim=0; dim<2; ++dim)
            {
                auto cmpMass = CmpTraits<Vector>::cmpArray(mass,dim);
                auto localCmpMass = cmpMass.slice(elements[e].indexVectorMask());
                
                auto cmpLocalPsGradW = CmpTraits<Vector>::cmpArray(localPsGradW,dim);
                
                auto cmpRhs = CmpTraits<Vector>::cmpArray(rhs,dim);
                auto localCmpRhs = cmpRhs.slice(elements[e].indexVectorMask());
                
                localCmpRhs  +=  cmpLocalPsGradW ;
                
                if(i == 0)
                    localCmpMass += massVector;
            }
        }
        
        U = U  + rhs / mass;
        q = q - alpha * div(U);
    }
}
*/

void advanceInTime(SEM::field::VectorField &U, SEM::field::DiscontinousField<SEM::Scalar> &p_prev, SEM::field::DiscontinousField<SEM::Scalar> &p_next, SEM::Scalar dt, SEM::Scalar nu, SEM::Scalar alpha1, SEM::Scalar alpha2)
{
     using namespace SEM;
     using namespace SEM::field;
     using namespace SEM::mesh;

     const Mesh & elements = U.mesh();

     static VectorField ddtU("ddtU", U.mesh(), SEM::iomanagment::NO_READ, SEM::iomanagment::NO_WRITE);

     numArray<Vector> mass(U.size(), 0.);
     numArray<Vector> rhs(U.size(), 0.);
     
     ddtU = (U.oldField(0)-U.oldField(1) )/dt;
        
     DiscontinousField<Scalar> forcing = p_prev - alpha1*div(ddtU) - alpha2*div(U);
 
     for( size_t e=0; e < elements.size(); ++e )
     {
        numArray2D<Scalar> eStiff=elements[e].stiffMatrix(nu);
        numArray2D<Scalar>::const_mappedArray massVector = elements[e].massMatrix().sliceArray( elements[e].localNodesInMatrix() );
        numArray<Vector> localPsGradW;
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
            
            numArray<Scalar> convection(SEM::array::matrixMul(derivMat,localCmpU));
            
            localCmpRhs  +=  cmpLocalPsGradW - SEM::array::matrixMul(eStiff,localCmpU) - convection * massVector;
            
            localCmpMass += massVector;
        }
     }
     
     ddtU = rhs / mass;

     U = U.oldField(0) + dt*ddtU;
     p_next = p_prev - ( alpha1*div(ddtU) + alpha2*div(U)); 
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
    using boost::shared_ptr;

    int srmOrder = 4;
    
    //--------------------- SETUP CASE --------------------------------------
    Case::setup(argc, args);
    
    //--------------------- CREATE MESH -------------------------------------
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    mesh.writeSpectralNodes(Case::nodesPath());
   
    //--------------------- CREATE FIELDS -----------------------------------
    numArray<shared_ptr<VectorField> > Us(4);
    numArray<shared_ptr<DiscontinousField<Scalar> > > ps_old(4);
    numArray<shared_ptr<DiscontinousField<Scalar> > > ps_new(4);
    ScalarField p("p",mesh,NO_READ);
    ScalarField divergence("divU",mesh,NO_READ);

    for(int i=0; i<srmOrder; ++i)
    {
        if(i!=srmOrder-1)
        {
            std::stringstream uName("U");
            uName << "_"<<i;
            Us[i]= boost::shared_ptr<VectorField>( new VectorField(uName.str(),mesh,NO_READ, NO_WRITE) );
        }
        else
        {
            Us[i] = boost::shared_ptr<VectorField>( new VectorField("U",mesh) );
        }

        ps_old[i] = boost::shared_ptr<DiscontinousField<Scalar> >( new DiscontinousField<Scalar>(mesh) );
        ps_new[i] = boost::shared_ptr<DiscontinousField<Scalar> >( new DiscontinousField<Scalar>(mesh) );
    }

    
    las::applyDirichletBCToSolution(*Us[srmOrder-1]);
    (*const_cast<numArray<Vector>*>(&Us[srmOrder-1]->oldField(0))) = *Us[srmOrder-1];

    for(int i=0; i<srmOrder-1; ++i)
    {
        (*Us[i]) == (*Us[srmOrder-1]);    
        (*const_cast<numArray<Vector>*>(&Us[i]->oldField(0))) = *Us[srmOrder-1];
    }


        
    //--------------------- GET COEFFS -----------------------------------
    Scalar dt = Case::time().timeStep();
    
    Scalar nu = Case::material()["nu"];
    
    //--------------------- SOLVE EQUATION -----------------------------------
    Scalar eps = 1e6;
    Scalar alpha1=1;
    Scalar alpha2=1;

    DiscontinousField<Scalar> zeroP(mesh);
    zeroP = 0;
       
    //solve time loop
    while(!Case::time().end())
    {
        ++Case::time();
       
        for(int s=0; s< srmOrder; ++s) 
        {
            if(s==0)
            {
                advanceInTime(*Us[s], zeroP, *ps_new[s], dt, nu, alpha1, alpha2);
            }
            else
            {
                advanceInTime(*Us[s], *ps_old[s-1], *ps_new[s], dt, nu, alpha1, alpha2);
            }
        }
        
        p = *ps_new[srmOrder-1];
	    divergence = div(*Us[srmOrder-1]);
        
        for(int s=0; s<srmOrder; ++s)
            (*ps_old[s]) = (*ps_new[s]);

        Case::time().fireWriting();
    }
    
    return 0;
}
