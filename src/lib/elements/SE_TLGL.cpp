
#include <cmath>

#include "SE_TLGL.h"

#include "utilities/ArrayFunctions.h"
#include "utilities/VectorUtils.h"

namespace SEM{ namespace mesh{

    REGISTER_IMPLEMENATION(SpectralElement,SE_TLGL,"SE_TLGL")
    
    SE_TLGL::SE_TLGL(): SpectralElement(5,"Tri")
    {
        this->generateSpectralNodesAndWeights();
        this->numberOfInteriorNodes = generateNumberOfInteriorNodes();
        this->numberOfInteriorInEdgesNodes = generateNumberOfInteriorInEdgesNodes();
        this->generateDerivativeMatrix();
        generateLocalIndexesInMatrix(m_localIndexesInMatrix);
        m_barycentricInterpWeights = computeInterpolationBarycentricWeights(nodes);
    }

    SE_TLGL::SE_TLGL(unsigned int NSpec): SpectralElement(NSpec,"Tri")
    {
        generateSpectralNodesAndWeights();
        numberOfInteriorNodes = generateNumberOfInteriorNodes();
        numberOfInteriorInEdgesNodes = generateNumberOfInteriorInEdgesNodes();
        this->generateDerivativeMatrix();
        generateLocalIndexesInMatrix(m_localIndexesInMatrix);
        m_barycentricInterpWeights = computeInterpolationBarycentricWeights(nodes);
    }

    SE_TLGL::~SE_TLGL()
    {
    }

    unsigned int SE_TLGL::numberOfRealNodes() const 
    {
        return 3;
    }
    
    void SE_TLGL::generateSpectralNodesRealCoordinates
    (
    Matrix<Vector>::type & rSCoords,
    const std::vector<Vector> & rCoords
    )
    {
        using std::vector;
        using boost::array;

        rSCoords.resize(NSpec);

        for (int i=0;i<NSpec;i++)
        {
            //First column is bigger-> colapsed edge to vertex
            i!=0 ? rSCoords[i].resize(NSpec-1) : rSCoords[i].resize(NSpec);
                
            for(int j=0;j<rSCoords[i].size();++j)
            {
                rSCoords[i][j][0] = (1+nodes[i])*(1-nodes[j])*(rCoords[1][0]-rCoords[0][0])/4 + (1+nodes[j])*(rCoords[2][0]-rCoords[0][0])/2 + rCoords[0][0];
                rSCoords[i][j][1] = (1+nodes[i])*(1-nodes[j])*(rCoords[1][1]-rCoords[0][1])/4 + (1+nodes[j])*(rCoords[2][1]-rCoords[0][1])/2 + rCoords[0][1];
            }
        }
//         rSCoords[0][NSpec-1][0] = rCoords[2][0];
//         rSCoords[0][NSpec-1][1] = rCoords[2][1];
    }

    void SE_TLGL::generateSpectralNodesAndWeights()
    {
        nodes.resize(NSpec);
        weights.resize(NSpec);

        const double pi = std::atan(1.)*4;
        int i,j;

        int Nmax=10; //maximal number of Newton solver's iterations to finde roots
        double eps=0.00001;//tolerance for finding poitns

        double q;
        double dq;
        double LN;

        if(NSpec==2)
        {
            nodes[0]=-1.;
            nodes[1]=1.;
            weights[0]=1.;
            weights[1]=weights[0];
        }
        else
        {
            nodes[0]=-1.;
            weights[0]=2./((NSpec-1.)*NSpec);
            nodes[NSpec-1]=1.;
            weights[NSpec-1]=weights[0];


            for (j=1;j<=NSpec/2-1.;j++)
            {
                nodes[j] = -std::cos((j+1./4)* pi/(NSpec-1) - 3./(8.*(NSpec-1.)*pi*(j+1./4)));
                //Newton iteration:
                for(i=0;i<Nmax;i++)
                {
                    qdqLN(nodes[j],NSpec,&q,&dq,&LN);
                    nodes[j] = nodes[j] - q/dq;
                    if (std::abs(q/dq)<eps) break;
                }
                qdqLN(nodes[j],NSpec,&q,&dq,&LN);
                nodes[NSpec-1-j] =-nodes[j];
                weights[j] = 2./((NSpec-1.)*NSpec*LN*LN);
                weights[NSpec-1-j]=weights[j];
            }
        }

        if ((NSpec-1)%2==0)
        {
            qdqLN(0.0,NSpec,&q,&dq,&LN);
            nodes[(NSpec-1)/2]=0;
            weights[(NSpec-1)/2]=2./((NSpec-1.)*NSpec*LN*LN);
        }
    }

    void SE_TLGL::generateDerivativeMatrix()
    {
        D_Matrix.resize(NSpec);
        int i,j;
        for (i=0;i<NSpec;i++)
        {
            D_Matrix[i] =  std::vector<double>(NSpec);
            
            for(j=0;j<NSpec;j++)
            {
                if (i!=j)
                {
                    D_Matrix[i][j] = WJ(j)/WJ(i)*(1./(nodes[i]-nodes[j]));
                    //Diagonal using "Negative Sum Trick"
                    D_Matrix[i][i]-=D_Matrix[i][j];
                }
            }
        }
    }
    
    Vector SE_TLGL::mapToLocal(const Vector& point, const std::vector< Vector >& nodes) const 
    {
        Vector Ksi;
        Scalar numerator   = ( point.x() - nodes[0].x() )*( nodes[2].y() - nodes[0].y() ) - ( point.y() - nodes[0].y() )*( nodes[2].x() - nodes[0].x() );
        Scalar denominator = ( point.x() - nodes[2].x() )*( nodes[1].y() - nodes[0].y() ) - ( point.y() - nodes[2].y() )*( nodes[1].x() - nodes[0].x() );
        Ksi.x() = 2. * numerator/denominator - 1.;
        
        numerator   = ( point.x() - nodes[0].x() )*( nodes[1].y() - nodes[0].y() ) - ( point.y() - nodes[0].y() )*( nodes[1].x() - nodes[0].x() );
        denominator = ( nodes[2].x() - nodes[0].x() )*( nodes[1].y() - nodes[0].y() ) - ( nodes[2].y() - nodes[0].y() )*( nodes[1].x() - nodes[0].x() );
        Ksi.y() = 2. * numerator/denominator - 1.;
        
        return Ksi;
    }
    
    numArray< Scalar > SE_TLGL::computeInterpolationCoefficients(Scalar ksi, Scalar etha) const 
    {
        numArray<Scalar> coeffs((NSpec-1)*NSpec+1,0.);
        
        size_t ksiNode=-1;
        size_t ethaNode =-1;
        
        for(size_t i=0; i<NSpec; ++i)
        {
            if(std::abs(nodes[i]-ksi) < 1e-6)
            {
                ksiNode=i;
            }
            if(std::abs(nodes[i]-etha) < 1e-6)
            {
                ethaNode=i;
            }
        }
        
        
        Scalar ksiSum = 0;
        Scalar ethaSum = 0;
        for(size_t k=0; k<NSpec; ++k)
        {
            ksiSum += m_barycentricInterpWeights[k] / (ksi-nodes[k]);
            ethaSum += m_barycentricInterpWeights[k] / (etha-nodes[k]);
        }
        
        //Calculate coeffs for all nodes despite colapsed one
        for(size_t i=0; i<NSpec; ++i)
        {
            Scalar ksiCoeff=0.;
            if(ksiNode!=-1)
            {
                if(ksiNode == i) ksiCoeff =1.;
            }
            else
            {
                ksiCoeff = m_barycentricInterpWeights[i] /( (ksi-nodes[i])*ksiSum );
            }
            
            for(size_t j=0; j<NSpec-1; ++j)
            {
                Scalar ethaCoeff=0.;
                if(ethaNode!=-1)
                {
                    if(ethaNode == j) ethaCoeff =1.;
                }
                else
                {
                    ethaCoeff = m_barycentricInterpWeights[j] /( (etha-nodes[j])*ethaSum );
                }   
                
                coeffs[i + j*NSpec] = ksiCoeff*ethaCoeff;
            }
        }
        
        //colapsed node 
        if(ethaNode == NSpec-1)
        {
            coeffs[(NSpec-1)*NSpec] = 1.;
        }
        else
        {
            coeffs[(NSpec-1)*NSpec] = m_barycentricInterpWeights[NSpec-1] /( (etha-nodes[NSpec-1])*ethaSum );
        }
        
        return coeffs;
    }
    

    std::vector<int> SE_TLGL::getEdgeLocalNodesIdInMatrix(int edgeId) const
    {
        std::vector<int> edge(NSpec);
        switch (edgeId)
        {
            case 0:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = i;
                }
                break;
            case 1:
                for(int i=0;i<NSpec-1;i++)
                {
                    edge[i] = (i+1)*NSpec-1;
                }
                edge[NSpec-1] = (NSpec-1)*NSpec;
                break;
            case 2:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = NSpec*(NSpec-1-i);
                }
                break;
        }

        return edge;
    }

    std::vector<int> SE_TLGL::getEdgeGloblaNodesId(int edgeId, const Matrix<int>::type & mask) const
    {
        std::vector<int> edge(NSpec);

        switch (edgeId)
        {
            case 0:
                for(int i=0;i<NSpec;++i)
                {
                    edge[i] = mask[i][0];
                }
                break;
            case 1:
                for(int i=0;i<NSpec-1;++i)
                {
                    edge[i] = mask[NSpec-1][i];
                }
                edge[NSpec-1] = mask[0][NSpec-1];
                break;
            case 2:
                for(int i=0;i<NSpec;++i)
                {
                    edge[i] = mask[0][NSpec-1-i];
                }
                break;
        }

        return edge;
    }

    std::vector<boost::array<int,2> > SE_TLGL::getEdgeLocalNodesId(int edgeId) const
    {
        std::vector<boost::array<int,2> > edge(NSpec);
        switch (edgeId)
        {
            case 0:
                for(int i=0;i<NSpec;++i)
                {
                    edge[i][0] = i;
                    edge[i][1] = 0;
                }
                break;
            case 1:
                for(int i=0;i<NSpec-1;++i)
                {
                    edge[i][0] = NSpec-1;
                    edge[i][1] = i;
                }
                edge[NSpec-1][0]=0;
                edge[NSpec-1][0]=NSpec-1;

                break;
            case 2:
                for(int i=0;i<NSpec;++i)
                {
                    edge[i][0] = 0;
                    edge[i][1] = NSpec-1-i;
                }

                break;
            default: 
                break;
        }

        return edge;
    }

    void SE_TLGL::defineMaskArraySize(Matrix<int>::type & mask)
    {
        mask.resize(NSpec);
        
        for( int i=0;i<NSpec;i++ )
        {
            i!=0 ? mask[i].resize(NSpec-1) : mask[i].resize(NSpec);
        }

    }

    void SE_TLGL::allocateMaskArray_VertexNode
    (
        int vertexNr,
        int nodeGlobalAdress,
        Matrix<int>::type & addressArray
    )
    {
        switch(vertexNr)
        {
        case 0:
            addressArray[0][0]=nodeGlobalAdress;
            break;
        case 1:
            addressArray[NSpec-1][0]=nodeGlobalAdress;
            break;
        case 2:
            addressArray[0][NSpec-1]=nodeGlobalAdress;
            break;
        }
    }

    void SE_TLGL::allocateMaskArray_NodesInEdge//Only interior in edge!!!
    (
        int edgeNr,
        std::vector<int> edgeGlobalAdress,
        Matrix<int>::type & addressArray
    )
    {
        int i;

        switch(edgeNr) //counterclockwise numeration on each edge ( in this way it is given also )
        {
        case 0:
            for(i=1;i<NSpec-1;++i)// from left to right, on bottom edge
            {
                addressArray[i][0]=edgeGlobalAdress[i-1];
            }
            break;
        case 1:
            for(i=1;i<NSpec-1;++i)// from bottom to top, on right edge
            {
                addressArray[NSpec-1][i]=edgeGlobalAdress[i-1];
            }
            break;
        case 2:
            for(i=1;i<NSpec-1;++i)// from top to bottm, on left edge
            {
                addressArray[0][(NSpec-1)-i]=edgeGlobalAdress[i-1];
            }
            break;
        }
    }

    void SE_TLGL::allocateMaskArray_InteriorNodes(int startAddress, Matrix<int>::type & addressArray)
    {
        for(int i=1;i<NSpec-1;i++)
        {
            for(int j=1;j<NSpec-1;j++)
            {
                addressArray[j][i]=startAddress;
                startAddress++;
            }
        }
    }

    void SE_TLGL::transformMatrixMaskToVectorMask(const Matrix<int>::type &mask, std::vector<int> &vecMask) const
    {
        vecMask.clear();
        vecMask.reserve(NSpec*(NSpec-1)+1);

        for(int i=0;i<NSpec-1;++i) //column in element nodes
            for(int j=0; j<NSpec;++j) //row in element nodes
                vecMask.push_back(mask[j][i]);

        vecMask.push_back(mask[0][NSpec-1]);
    }
    
    Vector SE_TLGL::normal(unsigned int edgeId, const std::vector<Vector> &rC) const
    {
        Vector dS;
        switch(edgeId)
        {
            case 0:
                dS = rC[1]-rC[0];
                break;
            case 1:
                dS = rC[2]-rC[1];
                break;
            case 2:
                dS = rC[0]-rC[2];
                break;
            default:
                ErrorInFunction << "Can't obtain normal for edge="<<edgeId<<", while triangle element has only 3 edges"<<iomanagment::endProgram;
        }
        dS.normalize();
        return Vector(dS.y(),-dS.x());
    }
    
    //Fill up the "H-Matrix"---> this is the evaluation of this parts of equation in laplacian discret week form:
    // H_Matrix[i][j][0] = (ksix*ksix)*J*W[i]*W[j];          H_Matrix[i][j][3] = (ksiy*ksiy)*J*weights[i]*weights[j];
    // H_Matrix[i][j][1] = (ksix*etax)*J*W[i]*W[j];          H_Matrix[i][j][4] = (ksiy*etay)*J*weights[i]*weights[j];
    // H_Matrix[i][j][2] = (etax*etax)*J*W[i]*W[j];          H_Matrix[i][j][5] = (etay*etay)*J*weights[i]*weights[j];
    // i,j - element spectral node index
    // W - weights
    // J - jacobian at i,j
    void SE_TLGL::generateH_Matrix(HMatrix & H_Matrix, const std::vector<Vector>& rC)
    {
        using std::vector;
        using boost::array;

        double G1,G2,G3,G4,G5,G6;

        double jakCoeff = ( (rC[1].x()-rC[0].x())*(rC[2].y()-rC[0].y()) - (rC[2].x()-rC[0].x())*(rC[1].y()-rC[0].y()) )/8;

        H_Matrix.resize(NSpec);
        for(int i=0;i<NSpec;++i)
        {
            H_Matrix[i].resize(NSpec);

            for(int j=0;j<NSpec;++j)//note: j!=NSpec-1 --> nodes[j] is never here equal to 0
            {
                if(j==NSpec-1)
                {
                    G1 = 0;//NEVER USED VALUE (need to be due to matrix dimmensions required by other coeficients which are used)
                    G2 = 0;//NEVER USED VALUE (need to be due to matrix dimmensions required by other coeficients which are used)
                    G2 = 0;//NEVER USED VALUE (need to be due to matrix dimmensions required by other coeficients which are used)
                    G3 = 0;//NEVER USED VALUE (need to be due to matrix dimmensions required by other coeficients which are used)
                    G4 = 0;//NEVER USED VALUE (need to be due to matrix dimmensions required by other coeficients which are used)
                    G5 = 0;//It is used! At j=Nspec-1 it is equal to 0. For clarity of equations we left this value in stiffnes matrix evaluation
                    G6 = 0;//It is used! At j=Nspec-1 it is equal to 0. For clarity of equations we left this value in stiffnes matrix evaluation
                }
                else
                { 
                    double xksi = (rC[1].x()-rC[0].x())*(1.-nodes[j])/4;
                    double yksi = (rC[1].y()-rC[0].y())*(1.-nodes[j])/4;
                    
                    double xetha=-(rC[1].x()-rC[0].x())*(1.+nodes[i])/4 + (rC[2].x()-rC[0].x())/2;
                    double yetha=-(rC[1].y()-rC[0].y())*(1.+nodes[i])/4 + (rC[2].y()-rC[0].y())/2;
                    
                    double J = std::abs( (1.-nodes[j])*jakCoeff );
                    
                    
                    G1 = (yetha*yetha)/J;
                    G2 =-(yetha*yksi) /J;
                    G3 = (yksi*yksi)  /J;
                    G4 = (xetha*xetha)/J;
                    G5 =-(xetha*xksi) /J;
                    G6 = (xksi*xksi)  /J;
                }

                //Obliczenie elementow macierzy H:
                double Wij = weights[i]*weights[j];
                
                H_Matrix[i][j][0] = G1 * Wij;
                H_Matrix[i][j][1] = G2 * Wij;
                H_Matrix[i][j][2] = G3 * Wij;
                
                H_Matrix[i][j][3] = G4 * Wij;
                H_Matrix[i][j][4] = G5 * Wij;
                H_Matrix[i][j][5] = G6 * Wij;

            }
        }
    }

    void SE_TLGL::generateM_Matrix(numArray2D<double> &M_Matrix, const std::vector<Vector> & rC ) const
    {
        M_Matrix.clear();
        M_Matrix.resize(NSpec);

        for(int p=0;p<NSpec;p++)
        {
            p!=0 ? M_Matrix[p].resize(NSpec-1) : M_Matrix[p].resize(NSpec);
        }

        for(int p=0;p<M_Matrix.size();p++)
        {
            for(int q=0;q<M_Matrix[p].size();q++)
            {
                // note :for p=0,q=N -->J==0 --> M[0][N] = 0;
                M_Matrix[p][q] = std::abs( jacobian(p,q,rC) )*weights[p]*weights[q]; 
            }
        }
    }

    void SE_TLGL::generateStiff_Matrix
    (
        Matrix4<double>::type & stiffMastrix,
        double scalarValue,
        const HMatrix & H_Matrix
    ) const
    {
        using std::vector;
        using boost::array;

        //Init matrix
        stiffMastrix.clear();
        stiffMastrix.resize(NSpec);
        
        //Resize matrix:
        for(int p=0; p<stiffMastrix.size(); ++p)
        {
            p!=0 ? stiffMastrix[p].resize(NSpec-1) : stiffMastrix[p].resize(NSpec);
            for(int q=0; q< stiffMastrix[p].size(); ++q)
            {
                stiffMastrix[p][q].resize(NSpec);
                for(int i=0; i< NSpec; ++i )
                {
                    i!=0 ? stiffMastrix[p][q][i].resize(NSpec-1) : stiffMastrix[p][q][i].resize(NSpec,0.);
                }
            }
        }

        //Fill matrix
        for(int p=0; p< stiffMastrix.size(); ++p)
        {
            for(int q=0; q< stiffMastrix[p].size(); ++q)
            {
                for(int i=0; i<NSpec; ++i)
                {
                    for(int j=0; j<NSpec; ++j)
                    {
                        if(q!=NSpec-1)
                        {
                            stiffMastrix[p][q][j][q] += scalarValue*(H_Matrix[i][q][0]+H_Matrix[i][q][3])*D_Matrix[i][p]*D_Matrix[i][j];
                            if(j!=NSpec-1)
                            {
                                stiffMastrix[p][q][p][j] += scalarValue*(H_Matrix[p][i][2]+H_Matrix[p][i][5])*D_Matrix[i][q]*D_Matrix[i][j];
                                stiffMastrix[p][q][i][j] += scalarValue*(H_Matrix[p][j][1]+H_Matrix[p][j][4])*D_Matrix[j][q]*D_Matrix[p][i];
                                stiffMastrix[p][q][i][j] += scalarValue*(H_Matrix[i][q][1]+H_Matrix[i][q][4])*D_Matrix[i][p]*D_Matrix[q][j];
                            }
                            else
                            {
                                stiffMastrix[p][q][0][NSpec-1] += scalarValue*(H_Matrix[p][i][2]+H_Matrix[p][i][5])*D_Matrix[i][q]*D_Matrix[i][NSpec-1];
                                stiffMastrix[p][q][0][NSpec-1] += scalarValue*(H_Matrix[i][q][1]+H_Matrix[i][q][4])*D_Matrix[i][p]*D_Matrix[q][NSpec-1];
                            }
                        }
                        else //q == NSpec-1
                        {
                            for(int k=0; k<NSpec; ++k)
                            {
                                if(k!=NSpec-1)
                                {
                                    stiffMastrix[p][q][i][k] += scalarValue*(H_Matrix[i][j][2]+H_Matrix[i][j][5])*D_Matrix[j][NSpec-1]*D_Matrix[j][k];
                                    stiffMastrix[p][q][j][k] += scalarValue*(H_Matrix[i][k][1]+H_Matrix[i][k][4])*D_Matrix[k][NSpec-1]*D_Matrix[i][j];
                                }
                                else
                                {
                                    stiffMastrix[p][q][0][NSpec-1] += scalarValue*(H_Matrix[i][j][2]+H_Matrix[i][j][5])*D_Matrix[j][NSpec-1]*D_Matrix[j][NSpec-1];
                                }
                            }//end k
                            
                        }
                        
                    }//end j
                }//end i
            }//end q
        }//end p

    }

    numArray2D<Scalar> SE_TLGL::generateStiff_Matrix
    (
            const double &scalarValue,
            const HMatrix &H_Matrix
    ) const
    {
        Matrix4<double>::type stiff_4d; 
        generateStiff_Matrix(stiff_4d,scalarValue, H_Matrix);

        //Resize matrix
        numArray2D<Scalar> matrix(NSpec*(NSpec-1)+1,NSpec*(NSpec-1)+1);

        for(int p=0; p<NSpec; ++p)
        {
            int size1;

            p!=0 ? size1=NSpec-1 : size1 =NSpec;

            for(int q=0; q<size1; ++q)
            {
                int row = q*NSpec + p;

                for(int i=0; i<NSpec; ++i)
                {
                    int size2;
                    i!=0 ? size2=NSpec-1 : size2 = NSpec;
                    
                    for(int j=0; j<size2; ++j)
                    {
                        int col = j*NSpec + i;
                        matrix[row][col] = stiff_4d[p][q][i][j];
                    }
                }
            }
        }
        return matrix;
    }
    
    Scalar SE_TLGL::x_ksi(int i, int j, const std::vector< Vector >& rC) const
    {
        return (rC[1].x()-rC[0].x())*(1.-nodes[j])/4.;
    }
    
    Scalar SE_TLGL::y_ksi(int i, int j, const std::vector< Vector >& rC) const
    {
        return (rC[1].y()-rC[0].y())*(1.-nodes[j])/4.;
    }
    
    Scalar SE_TLGL::x_etha(int i, int j, const std::vector< Vector >& rC) const
    {
        return -(rC[1].x()-rC[0].x())*(1.+nodes[i])/4. + (rC[2].x()-rC[0].x())/2.;
    }
    
    Scalar SE_TLGL::y_etha(int i, int j, const std::vector< Vector >& rC) const
    {
        return -(rC[1].y()-rC[0].y())*(1.+nodes[i])/4. + (rC[2].y()-rC[0].y())/2.;
    }
    
    Scalar SE_TLGL::jacobianCoeff(const std::vector< Vector >& rC) const 
    {
        return ( (rC[1].x()-rC[0].x())*(rC[2].y()-rC[0].y()) - (rC[2].x()-rC[0].x())*(rC[1].y()-rC[0].y()) )/8.;
    }
    
    Scalar SE_TLGL::jacobian(int i, int j, const std::vector< Vector >& rC) const
    {
        return (1.-nodes[j])*jacobianCoeff(rC);
    }
 
    Scalar SE_TLGL::ksi_x(int i, int j, const std::vector< Vector >& rC) const
    {
        if(j == NSpec-1)
            return 0.;
        else
            return y_etha(i,j,rC)/jacobian(i,j,rC);
    }
 
 
    Scalar SE_TLGL::ksi_y(int i, int j, const std::vector< Vector >& rC) const
    {
        if(j == NSpec-1)
            return 0.;
        else
            return -x_etha(i,j,rC)/jacobian(i,j,rC);
    }
    
    Scalar SE_TLGL::etha_x(int i, int j, const std::vector< Vector >& rC) const 
    {
        return - (rC[1].y()-rC[0].y())/(4.*jacobianCoeff(rC));
    }
    
    Scalar SE_TLGL::etha_y(int i, int j, const std::vector< Vector >& rC) const
    {
        return (rC[1].x()-rC[0].x())/(4.*jacobianCoeff(rC));
    }
    
    Scalar SE_TLGL::rotU(const numArray< Vector >& u, std::size_t ksi, std::size_t etha, const std::vector< Vector >& rC) const 
    {
        Scalar dvdksi=0;
        Scalar dudksi=0;
        Scalar dvdetha=0;
        Scalar dudetha=0;
        
        if(etha!=NSpec-1)
        {
            for(size_t k=0; k<NSpec; ++k)
            {
                dvdksi += u[k+etha*NSpec].y()*D_Matrix[ksi][k];
                dudksi += u[k+etha*NSpec].x()*D_Matrix[ksi][k];
            }
        }
        
        for(size_t k=0; k<NSpec; ++k)
        {
            dvdetha += ( k==NSpec-1 ? u[k*NSpec].y() : u[ksi+k*NSpec].y() ) * D_Matrix[etha][k];
            dudetha += ( k==NSpec-1 ? u[k*NSpec].x() : u[ksi+k*NSpec].x() ) * D_Matrix[etha][k];
        }
        
        Scalar l11 = ksi_x(ksi,etha,rC);
        Scalar l21 = etha_x(ksi,etha,rC);
        Scalar l12 = ksi_y(ksi,etha,rC);
        Scalar l22 = etha_y(ksi,etha,rC);
        
        return (dvdksi*l11 + dvdetha*l21) - (dudksi*l12 + dudetha*l22);
    }
    
    numArray< Scalar > SE_TLGL::rot(const numArray< Vector >& U, const std::vector< Vector >& rCoords) const 
    {
        numArray<Scalar> result(NSpec*(NSpec-1)+1,0.);
        for(size_t ksi=0; ksi<NSpec; ++ksi)
        {
            size_t s2 = ksi==0 ? NSpec : NSpec-1;
            for(size_t etha=0; etha<s2; ++etha)
            {
                result[ksi+etha*NSpec] = rotU(U,ksi,etha,rCoords);
            }
        }
        return result;
    }
    
    Scalar SE_TLGL::gamma(std::size_t node, unsigned int edge, const std::vector< Vector >& rC) const 
    {
        switch(edge)
        {
            case 0://bottom
            {
                Scalar xKsi=x_ksi(node,0,rC);
                Scalar yKsi=y_ksi(node,0,rC);
                return std::sqrt(xKsi*xKsi + yKsi*yKsi);
            }
            case 1://right
            {
                Scalar xEtha=x_etha(NSpec-1,node,rC);
                Scalar yEtha=y_etha(NSpec-1,node,rC);
                return std::sqrt( xEtha*xEtha + yEtha*yEtha );
            }
            case 2://left = reversed nodes order in edge acording to edge definition
            {
                Scalar xEtha=x_etha(0,NSpec-1-node,rC);
                Scalar yEtha=y_etha(0,NSpec-1-node,rC);
                return std::sqrt( xEtha*xEtha + yEtha*yEtha );
            }
            default:
                ErrorInFunction << "Wrong edge id="<<edge
                <<" for calulation edge integral jacobian - triangle element has only 3 edges" 
                <<iomanagment::endProgram;
        }
    }
    
    void SE_TLGL::generateWeekVddx
    (
        const numArray< Scalar >& u,
        const std::vector< Vector >& rCoords, 
        numArray< Scalar >& result
    ) const 
    {
        result.resize(NSpec *(NSpec-1) + 1, 0.);
        
        for(size_t p=0; p<NSpec; ++p)
        {
            size_t qLen = p==0? NSpec : NSpec-1;
            
            for(size_t q=0; q<qLen; ++q)
            {
                Scalar ksiSum=0;
                Scalar ethaSum=0;
                
                if(q != NSpec-1)
                {
                    for(size_t k=0; k<NSpec; ++k )
                    {
                        Scalar J = std::abs(jacobian(k,q,rCoords));
                        Scalar H11 =  ksi_x(k,q,rCoords)*J*weights[k]*weights[q];
                        ksiSum += u[k+q*NSpec]*D_Matrix[k][p]*H11;
                        
                        J= std::abs(jacobian(p,k,rCoords));
                        Scalar H21 = etha_x(p,k,rCoords) * J *weights[p] * weights[k];
                        
                        //when j=N jacobian is equal to 0, so integral is well defined even if value in top node is undefined
                        if(k!=NSpec-1)
                            ethaSum += u[p+k*NSpec] * D_Matrix[k][q] * H21;
                        else
                            ethaSum += u[0+k*NSpec] * D_Matrix[k][q] * H21; 
                    }
                }
                else //p=0, q=NSpec-1
                {
                    for(size_t i=0; i<NSpec; ++i)
                    {
                        for(size_t j=0; j<NSpec; ++j)
                        {
                            //when j=N jacobian is equal to 0, so integral is well defined even if value in top node is undefined
                            Scalar J = std::abs(jacobian(i,j,rCoords));
                            Scalar H21 = etha_x(i,j,rCoords) * J * weights[i] * weights[j];
                            
                            ethaSum += ( j!=NSpec-1 ? u[i+j*NSpec] : u[j*NSpec] ) * D_Matrix[j][q]*H21;
                        }
                    }
                }
                
                result[p+q*NSpec] = ksiSum + ethaSum;
            }
        }
    }
    
    void SE_TLGL::generateWeekVddy
    (
        const numArray< Scalar >& u,
        const std::vector< Vector >& rCoords, 
        numArray< Scalar >& result
    ) const 
    {
        result.resize(NSpec *(NSpec-1) + 1, 0.);
        
        for(size_t p=0; p<NSpec; ++p)
        {
            size_t qLen = p==0? NSpec : NSpec-1;
            
            for(size_t q=0; q<qLen; ++q)
            {
                Scalar ksiSum=0;
                Scalar ethaSum=0;
                
                if(q != NSpec-1)
                {
                    for(size_t k=0; k<NSpec; ++k )
                    {
                        Scalar J = std::abs(jacobian(k,q,rCoords));
                        Scalar H12 =  ksi_y(k,q,rCoords)*J*weights[k]*weights[q];
                        ksiSum += u[k+q*NSpec]*D_Matrix[k][p]*H12;
                        
                        //when j=N jacobian is equal to 0, so integral is well defined even if value in top node is undefined
                        J= std::abs(jacobian(p,k,rCoords));
                        Scalar H22 = etha_y(p,k,rCoords) * J *weights[p] * weights[k];
                        
                        if(k!=NSpec-1)
                            ethaSum += u[p+k*NSpec] * D_Matrix[k][q] * H22;
                        else
                            ethaSum += u[0+k*NSpec] * D_Matrix[k][q] * H22; 
                    }
                }
                else //p=0, q=NSpec-1
                {
                    for(size_t i=0; i<NSpec; ++i)
                    {
                        for(size_t j=0; j<NSpec; ++j)
                        {
                            //when j=N jacobian is equal to 0, so integral is well defined even if value in top node is undefined
                            Scalar J = std::abs(jacobian(i,j,rCoords));
                            Scalar H22 = etha_y(i,j,rCoords) * J * weights[i] * weights[j];
                            
                            ethaSum += ( j!=NSpec-1 ? u[i+j*NSpec] : u[j*NSpec] ) * D_Matrix[j][q]*H22;
                        }
                    }
                }
                
                result[p+q*NSpec] = ksiSum + ethaSum;
            }
        }
    }
    
    void SE_TLGL::generate_ddx
    (
        const numArray< Scalar >& u,
        const std::vector< Vector >& rCoords, 
        numArray< Scalar >& result
    ) const 
    {
        result.resize(NSpec*(NSpec-1)+1,0.);
        
        for(size_t i=0; i<NSpec; ++i)
        {
            size_t size2 = i !=0 ? NSpec-1 : NSpec;
            for(size_t j=0; j<size2; ++j)
            {
                Scalar sumKsi=0;
                Scalar sumEtha=0;
                
                if(j!=NSpec-1)
                {
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        sumKsi += u[j*NSpec+k]*D_Matrix[i][k];
                        sumEtha += ( k!=NSpec-1 ? u[k*NSpec+i] : u[k*NSpec] )*D_Matrix[j][k];
                    }
                    result[j*NSpec+i] = sumKsi * ksi_x(i,j,rCoords)  +  sumEtha * etha_x(i,j,rCoords); 
                }
                else//j==NSpec-1 --> base function depends only on etha, 
                    // and value singular because we don't know direction from which du/detha is evaluated
                    // Here we took value from left edge. This value in implementation
                    // shall not be directly used in calculations. This value shall be cleared by zero jacobian in 
                    // integral evaluation
                {
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        sumEtha += ( k!=NSpec-1 ? u[k*NSpec+i] : u[k*NSpec] )*D_Matrix[j][k];
                    }
                    result[j*NSpec] = sumEtha * etha_x(i,j,rCoords);
                }
            }
        }
    }
    
    void SE_TLGL::generate_ddy
    (
        const numArray< Scalar >& u,
        const std::vector< Vector >& rCoords, 
        numArray< Scalar >& result
    ) const 
    {
        result.resize(NSpec*(NSpec-1)+1,0.);
        
        for(size_t i=0; i<NSpec; ++i)
        {
            size_t size2 = i !=0 ? NSpec-1 : NSpec;
            
            for(size_t j=0; j<size2; ++j)
            {
                Scalar sumKsi=0;
                Scalar sumEtha=0;
                
                if(j!=NSpec-1)
                {
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        sumKsi += u[j*NSpec+k]*D_Matrix[i][k];
                        sumEtha += ( k!=NSpec-1 ? u[k*NSpec+i] : u[k*NSpec] )*D_Matrix[j][k];
                    }
                    result[j*NSpec+i] = sumKsi * ksi_y(i,j,rCoords)  +  sumEtha * etha_y(i,j,rCoords); 
                }
                else //j==NSpec-1 --> base function depends only on etha, 
                    // and value singular because we don't know direction from which du/detha is evaluated
                    // Here we took value from left edge. This value in implementation
                    // shall not be directly used in calculations. This value shall be cleared by zero jacobian in 
                    // integral evaluation
                {
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        sumEtha += ( k!=NSpec-1 ? u[k*NSpec+i] : u[k*NSpec] )*D_Matrix[j][k];
                    }
                    result[j*NSpec] += sumEtha * etha_y(i,j,rCoords) ;
                }
            }
        }
    }

    numArray< Scalar > SE_TLGL::boundInt_NxRotRotU(const numArray< Vector >& U, const std::vector< Vector >& rCoords, const size_t& edge) const 
    {
        numArray<Scalar> integral( (NSpec-1)*NSpec + 1 ,0.);
        
        Vector edgeNormal = normal(edge, rCoords);
        switch(edge)
        {
            case 0:
            {
                for(size_t p=0; p<NSpec; ++p)
                {
                    size_t size2 = p==0 ? NSpec : NSpec-1;
                    for(size_t q=0; q<size2; ++q)
                    {
                        size_t index = p + q*NSpec;
                        
                        if(q!=NSpec-1)
                        {
                            if(q==0)
                            {
                                for(size_t i=0; i<NSpec; ++i)
                                {
                                    Scalar l11=ksi_x(i,0,rCoords);
                                    Scalar l12=ksi_y(i,0,rCoords);
                                    Scalar nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                                    integral[index] +=rotU(U,i,0,rCoords)*D_Matrix[i][p]*nKsi*gamma(i,edge,rCoords)*weights[i];
                                }
                            }
                            Scalar l21=etha_x(p,0,rCoords);
                            Scalar l22=etha_y(p,0,rCoords);
                            Scalar nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                            integral[index] += rotU(U,p,0,rCoords)*D_Matrix[0][q]*nEtha*gamma(p,edge,rCoords)*weights[p];
                        }
                        else
                        {
                            Scalar sum=0;
                            for(size_t i=0; i<NSpec; ++i)
                            {
                                Scalar l21=etha_x(i,0,rCoords);
                                Scalar l22=etha_y(i,0,rCoords);
                                Scalar nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                                sum += rotU(U,i,0,rCoords)*nEtha*gamma(i,edge,rCoords)*weights[i];
                            }
                            integral[index] += sum*D_Matrix[0][NSpec-1];
                        }
                    }
                }
            }
                break;
            case 1:
            case 2:
            {
                size_t ksiId = edge == 1? NSpec-1 : 0;
                bool invertNodes = edge == 1 ? false : true; //inform if edge definitions assumes inverted nodes order
                
                for(size_t p=0; p<NSpec; ++p)
                {
                    size_t size2 = (p==0) ? NSpec : NSpec-1;
                    
                    for(size_t q=0; q<size2; ++q)
                    {   
                        size_t index = p + q*NSpec;
                        
                        if(q!=NSpec-1)
                        {
                            Scalar l11 = ksi_x(ksiId,q,rCoords);
                            Scalar l12 = ksi_y(ksiId,q,rCoords);
                            Scalar nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                            
                            Scalar G = invertNodes ? gamma(NSpec-1-q,edge,rCoords) : gamma(q,edge,rCoords);
                            
                            integral[index] += rotU(U,ksiId,q,rCoords) * D_Matrix[ksiId][p] * nKsi * G * weights[q];
                            
                            if(p == ksiId)
                            {
                                for(size_t j=0; j<NSpec; ++j)
                                {
                                    Scalar l21 = etha_x(ksiId,j,rCoords);
                                    Scalar l22 = etha_y(ksiId,j,rCoords);
                                    Scalar nEtha=edgeNormal.x()*l22 - edgeNormal.y()*l21;
                                    
                                    G = invertNodes ? gamma(NSpec-1-j,edge,rCoords) : gamma(j,edge,rCoords);
                                    
                                    //Scalar rotVel = j==NSpec-1 ? rotU(U,0,j,rCoords) : rotU(U,ksiId,j,rCoords);
                                    Scalar rotVel = rotU(U,ksiId,j,rCoords);
                                    integral[index] += rotVel * D_Matrix[j][q] *nEtha * G * weights[j];
                                }
                            }
                        }
                        else //p=0, q=NSpec-1
                        {
                            for(size_t j=0; j<NSpec; ++j)
                            {
                                Scalar l21 = etha_x(ksiId,j,rCoords);
                                Scalar l22 = etha_y(ksiId,j,rCoords);
                                Scalar nEtha=edgeNormal.x()*l22 - edgeNormal.y()*l21;
                                
                                Scalar G = invertNodes ? gamma(NSpec-1-j,edge,rCoords) : gamma(j,edge,rCoords);
                                
//                                 Scalar rotVel = j==NSpec-1 ? rotU(U,0,j,rCoords) : rotU(U,ksiId,j,rCoords);
                                Scalar rotVel = rotU(U,ksiId,j,rCoords);
                                integral[index] += rotVel * D_Matrix[j][NSpec-1] * nEtha * G * weights[j];
                            }
                        }
                    }
                }
            }
                break;
            default:
                ErrorInFunction << "Wrong edge id="<<edge<<" used with triangle element"<<iomanagment::endProgram;
        }
       
        return integral;
    }
    
    numArray2D< Scalar > SE_TLGL::convMatrix(const numArray< Vector >& U, const std::vector< Vector >& rC) const 
    {
        size_t size = NSpec*(NSpec-1)+1;
        numArray2D<Scalar> deriv(NSpec*(NSpec-1)+1,NSpec*(NSpec-1)+1,0.);
        
        for(size_t i=0; i<NSpec; ++i)
        {
            size_t size2 = i !=0 ? NSpec-1 : NSpec;
            
            for(size_t j=0; j<size2; ++j)
            {
                size_t row = i+j*NSpec;
                
                if(j!=NSpec-1)
                {
                    Scalar ksiCoeff = ksi_x(i,j,rC)*U[row].x()+ksi_y(i,j,rC)*U[row].y();
                    Scalar ethaCoeff = etha_x(i,j,rC)*U[row].x()+etha_y(i,j,rC)*U[row].y();
                    
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        deriv[row][k+j*NSpec] += ksiCoeff*D_Matrix[i][k];
                        
                        if(k!=NSpec-1)
                            deriv[row][i+k*NSpec] += ethaCoeff*D_Matrix[j][k];
                        else
                            deriv[row][k*NSpec] += ethaCoeff*D_Matrix[j][k];
                    }
                }
                else //j==NSpec-1 --> base function depends only on etha, 
                    // and value singular because we don't know direction from which du/detha is evaluated
                    // Here we took value from left edge. This value in implementation
                    // shall not be directly used at integration step. This value shall be cleared by zero jacobian in 
                    // integral evaluation
                {
                    Scalar ethaCoeff = etha_x(i,j,rC)*U[row].x()+etha_y(i,j,rC)*U[row].y();
                    
                    for(size_t k=0; k<NSpec; ++k)
                    {
                        if(k!=NSpec-1)
                            deriv[row][i+k*NSpec] += ethaCoeff*D_Matrix[j][k];
                        else
                            deriv[row][k*NSpec]  += ethaCoeff*D_Matrix[j][k];
                    }
                }
            }
        }
                
        return deriv;
        
    }
    
    void SE_TLGL::generateEdgeNeumanBCIntegralCoefficient
    (
        int edgeNumber,
        const std::vector<Vector> & rC,
        std::vector<double>& resultNodalVal
    )
    {
        resultNodalVal.resize(NSpec);
        switch(edgeNumber)
        {
        case 0:
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = gamma(i,edgeNumber,rC)*weights[i];
            }
            break;
        case 1:
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = gamma(i,edgeNumber,rC)*weights[i];
            }
            break;

        case 2:
            for(int i=0;i<NSpec;i++)
            {
                size_t node = NSpec-1-i;
                resultNodalVal[i] = gamma(node,edgeNumber,rC)*weights[node];
            }
            break;
        }
    }


    int SE_TLGL::generateNumberOfInteriorNodes()
    {
        return (NSpec-2)*(NSpec-2);
    }

    void SE_TLGL::generateLocalIndexesInMatrix(VectorToMatrixMap &map) const
    {
        map.resize(NSpec*(NSpec-1)+1);

        for(int i=0; i<NSpec; ++i)
        {
            int size;
            i!=0 ? size = NSpec -1 : size = NSpec;
            
            for(int j=0; j<size; ++j)
            {
                map[j*NSpec+i][0]=i;
                map[j*NSpec+i][1]=j;
            }
        }
    }

    std::vector<int> SE_TLGL::generateNumberOfInteriorInEdgesNodes()
    {
        return std::vector<int> (3,NSpec-2);
    }

    void SE_TLGL::qdqLN(double x, int N, double*q, double*dq, double*LN)
    {
        double LN2,LN1,LN11;
        double dLN,dLN2,dLN1,dLN11;
        int k;

        LN2=1;
        LN1=x;
        dLN2=0;
        dLN1=1;

        for (k=2;k<N;k++)
        {
            *LN=(2.*k-1.)/k*x*LN1-(k-1.)/k*LN2;
            dLN = dLN2+(2.*k-1.)*LN1;
            LN2=LN1;
            LN1=*LN;
            dLN2=dLN1;
            dLN1=dLN;
        }

        k=N+1;
        LN11=(2.*k-1.)/k*x*LN1-(k-1.)/k*LN2;
        dLN11=dLN2+(2.*k-1.)*LN1;
        *q=LN11-LN2;
        *dq=dLN11-dLN2;
    }

    double SE_TLGL::WJ(int j)
    {
        double w=1;
        for (int i=0;i<nodes.size();i++)
        {
            if (i!=j)
            {
                w*=1./(nodes[j]-nodes[i]);
            }
        }
        return w;
    }

    std::ostream& operator<<(std::ostream& o, SE_TLGL& e)
    {
        o<<"========== Spectral Element::TLGL=========="<<std::endl<<"number of nodes in element:\t"<<e.getSpectralSize()<<std::endl;
        o<<"Nodes:"<<std::endl<<"[";
        for(int i=0;i<e.getSpectralSize();i++)
        {
            o<<e.getSpectralNodes()[i];
            if(i!=e.getSpectralSize())
                o<<";";
        }
        o<<"]"<<std::endl<<"Weights"<<std::endl<<"[";
        for(int i=0;i<e.getSpectralSize();i++)
        {
            o<<e.getSpectralWeights()[i];
            if(i!=e.getSpectralSize())
                o<<";";
        }
        o<<"]"<<std::endl<<"=========================================="<<std::endl;

        return o;
    }
    
    

} //mesh
} //SEM
