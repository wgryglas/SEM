#include "iomanagment/case.h"
#include "mesh/Mesh.h"
#include "fields/Fields.h"
#include "fieldMath/ddt.h"
#include "fieldMath/laplacian.h"
#include "solver/Solver.h"
#include "utilities/ArrayFunctions.h"
#include "fieldMath/postCalculation.h"
#include "solver/EquationAssignment.h"
#include "fieldMath/WeakFormDefinition.h"
#include "fieldMath/DivergenceOperator.h"
#include "fieldMath/RotOperator.h"
#include "fieldMath/ConvectiveDerivative.h"

int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::field;
    using namespace SEM::mesh;
    using namespace SEM::weak;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    
    //Create T field
    ScalarField T
    (
        Case::time(),
        new RegistryFile
        (
            Case::time(),
            "T",
            Case::time().localPath(),
            iomanagment::NO_READ,
            iomanagment::AUTO
        ),
        mesh
    );
    
    //Create post grad T field
    VectorField U
    (
        Case::time(),
        new RegistryFile
        (
            Case::time(),
            "U",
            Case::time().localPath(),
            iomanagment::NO_READ,
            iomanagment::AUTO
        ),
        mesh
    );
    
    VectorField gradT(mesh);
    ScalarField divU
    (
        Case::time(),
        new RegistryFile
        (
            Case::time(),
            "divU",
            Case::time().localPath(),
            iomanagment::NO_READ,
            iomanagment::AUTO
        ),
        mesh
    );
    
    
    
    
    numArray<Scalar> X = xCmps(mesh.spectralNodes());
    numArray<Scalar> Y = yCmps(mesh.spectralNodes());
    
    using SEM::array::sin;
    using SEM::array::cos;
    using SEM::array::min;    
    using SEM::array::max;
    using SEM::array::abs;
    using SEM::array::mag;
    
    T = sin(X)*cos(Y);
    xCmps(U) = sin(X)*cos(Y);
    yCmps(U) =-cos(X)*sin(Y);
    
    gradT = grad(T);
    numArray<Vector> gradTsol(X.size());
    xCmps(gradTsol) = cos(X)*cos(Y);
    yCmps(gradTsol) =-1.*sin(X)*sin(Y);
    Scalar gradTError = max(mag(gradTsol-gradT));

    divU = div(U);
    Scalar divUError = max(abs(divU));
    
    std::cout << "gradient Linf error = "<<gradTError <<std::endl;
    std::cout << "divergence Linf error = "<<divUError <<std::endl;
    std::cout << "divergence avg value = "<<SEM::array::avg(divU) <<std::endl;
    
//     DiscontinousField<Scalar> dU = div(U);
//     std::cout << std::endl;
//     for(size_t e=0; e< mesh.size(); ++e)
//     {
//         iomanagment::write(dU[e],std::cout)<<std::endl;
//     }
//     
//     T = X;
//     DiscontinousField<Vector> gT = grad(T);
//     std::cout<<"grd T"<<std::endl;
//     for(size_t e=0; e< mesh.size(); ++e)
//     {
//         auto xL = X.slice(mesh[e].indexVectorMask());
//         auto yL = Y.slice(mesh[e].indexVectorMask());
//         numArray<Scalar> gradXErr = cos(xL)*cos(yL) - xCmps(gT[e]);
//         iomanagment::write(gradXErr,std::cout)<<std::endl;
//     }
    
    
    
    
    divU.registryFile()->fireWriting();
    
    // Check integrated weak gradient of 2 diffrent forms: direct diff on field T, and  
    // moving grad onto test function
    std::vector<numArray<Scalar> > rhsVector1;
    rhsVector1.resize(2);
    rhsVector1[0].resize(T.size(),0.);
    rhsVector1[1].resize(T.size(),0.);
    
    std::vector<numArray<Scalar> > rhsVector2;
    rhsVector2.resize(2);
    rhsVector2[0].resize(T.size(),0.);
    rhsVector2[1].resize(T.size(),0.);
    
    a(grad(T),phi)->buildExplicit( mesh, rhsVector1, las::AddEqAssigment() );
//     a(gradTsol,phi)->buildExplicit( mesh, rhsVector1, las::AddEqAssigment() );
    a(T,grad(phi))->buildExplicit( mesh, rhsVector2, las::AddEqAssigment() );
    
    std::cout<<"max diff between a(grad(T),phi) and a(T,grad(Phi)): ("<<max(abs(rhsVector1[0]-rhsVector2[0])) << ", "<<max(abs(rhsVector1[1]-rhsVector2[1]))<<std::endl;
    std::cout<<"/////////////////////////"<<std::endl;
    
    numArray<bool> boundaryMask(T.size(),false);
    for(auto itr=mesh.boundaryMesh().begin();itr!=mesh.boundaryMesh().end();++itr)
    {
        auto edgeItr=itr->begin();
        
        for(auto edgeItr=itr->begin(); edgeItr!=itr->end(); ++edgeItr)
        {
            boundaryMask.slice(edgeItr->gNodesIds())=true;
        }
    }
    
//     std::cout <<"diff.x=int(T*grad(phi)).x + int(grad(T),phi).x > 1e-5"<<std::endl;     
//     for(size_t i=0; i<T.size(); ++i)
//     {
//         Scalar diff = rhsVector2[0][i] + rhsVector1[0][i];
//         if(std::abs(diff) > 1e-5 )
//         {
//            std::cout <<"diff.x("<<i<<")="<<diff<<", boundary? "<<boundaryMask[i]<<",int(T*grad(phi)).x="<<rhsVector2[0][i]<<", int(grad(T)*phi).x="<<rhsVector1[0][i]<<std::endl;
//         }
//     }
    
        std::cout <<"diff.x=int(T*grad(phi)).x + int(grad(T),phi).x > 1e-5 outside boundary "<<std::endl; 
        for(size_t e=0; e<mesh.size(); ++e)
        {
            const auto &mask = mesh[e].indexVectorMask();
            std::cout <<"element "<<e<<", with "<<mask.size()<<" nodes"<<std::endl;
            for(size_t n=0; n<mask.size(); ++n)
            {
                Scalar diff = rhsVector2[0][mask[n]] + rhsVector1[0][mask[n]];
                if(std::abs(diff) > 1e-5 && !boundaryMask[mask[n]])
                {
                    std::cout <<"diff("<<n+1<<").x="<<diff<<",int(T*grad(phi)).x="<<rhsVector2[0][mask[n]]<<", int(grad(T)*phi).x="<<rhsVector1[0][mask[n]]<<std::endl;
                }
                
                diff = rhsVector2[1][mask[n]] + rhsVector1[1][mask[n]];
                if(std::abs(diff) > 1e-5 && !boundaryMask[mask[n]])
                { 
                    std::cout <<"diff("<<n+1<<").y="<<diff<<",int(T*grad(phi)).y="<<rhsVector2[1][mask[n]]<<", int(grad(T)*phi).y="<<rhsVector1[1][mask[n]]<<std::endl;
                }
            }
        }
    
    
    // rotation test
    
    DiscontinousField<Scalar> rotU = rot(GeometricFieldElementProxy<Vector>(U));
    
    numArray<Scalar> rotUSol= 2.*sin(X)*sin(Y);
    for(size_t e=0; e<mesh.size(); ++e)
    {
        numArray<Scalar> err = rotUSol.slice(mesh[e].indexVectorMask()) - rotU[e];
        if(SEM::array::max(SEM::array::abs(err))>1e-5)
            iomanagment::write(err,std::cout << "element "<<e<<"\n")<<std::endl;
    }
    
    
    //convective deriv test:
    std::cout<<"Convective derivative test:"<<std::endl;
    DiscontinousField<Vector> cDer=cDeriv(U,U);
    numArray<Scalar> cDerSol_x = sin(X)*cos(X);
    numArray<Scalar> cDerSol_y = sin(Y)*cos(Y);
    for(size_t e=0; e< mesh.size(); ++e)
    {
        auto map =mesh[e].indexVectorMask();
        numArray<Scalar> errX = cDerSol_x.slice(map) - xCmps(cDer[e]);
        numArray<Scalar> errY = cDerSol_y.slice(map) - yCmps(cDer[e]);
        
        if(SEM::array::max(SEM::array::abs(errX)) > 1e-5)
        {
            std::cout<<"error x in element "<<e<<std::endl;
            for(size_t n=0; n<errX.size(); ++n)
            {
                std::cout<<"err="<<errX[n]<<", cDerSolx="<<cDerSol_x[map[n]]<<", cDerNum="<<cDer[e][n].x()<<std::endl;
            }
        }
        
        if(SEM::array::max(SEM::array::abs(errY)) > 1e-5)
        {
            std::cout<<"error y in element "<<e<<std::endl;
            for(size_t n=0; n<errX.size(); ++n)
            {
                std::cout<<"err="<<errY[n]<<", cDerSolx="<<cDerSol_y[map[n]]<<", cDerNum="<<cDer[e][n].y()<<std::endl;
            }
        }
        
    }
    
    
    
    return 0;
    
}




