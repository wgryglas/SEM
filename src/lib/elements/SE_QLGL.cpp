#include "SE_QLGL.h"

#include <cmath>

#include "utilities/ArrayFunctions.h"

namespace SEM{ namespace mesh{
    
    
    REGISTER_IMPLEMENATION(SpectralElement,SE_QLGL,"SE_QLGL")
    
    
    SE_QLGL::SE_QLGL():SpectralElement(5,"Quad")
    {
        this->generateSpectralNodesAndWeights();
        this->numberOfInteriorNodes = generateNumberOfInteriorNodes();
        this->numberOfInteriorInEdgesNodes = generateNumberOfInteriorInEdgesNodes();
        this->generateDerivativeMatrix();
        generateLocalIndexesInMatrix(m_localIndexesInMatrix);
        m_barycentricInterpWeights = computeInterpolationBarycentricWeights(nodes);
    }

    SE_QLGL::SE_QLGL(unsigned int NSpec):SpectralElement(NSpec,"Quad")
    {
        generateSpectralNodesAndWeights();
        numberOfInteriorNodes = generateNumberOfInteriorNodes();
        numberOfInteriorInEdgesNodes = generateNumberOfInteriorInEdgesNodes();
        this->generateDerivativeMatrix();
        generateLocalIndexesInMatrix(m_localIndexesInMatrix);
        m_barycentricInterpWeights = computeInterpolationBarycentricWeights(nodes);
    }

    SE_QLGL::~SE_QLGL()
    {
    }

    unsigned int SE_QLGL::numberOfRealNodes() const 
    {
        return 4;
    }
    
    void SE_QLGL::generateSpectralNodesRealCoordinates
    (
        Matrix<Vector>::type& rSCoords,
        const std::vector<Vector>& rCoords
    )
    {
        using std::vector;
        using boost::array;

        rSCoords.resize(NSpec);
        for (int i=0;i<NSpec;i++)
        {
            rSCoords[i] = vector<Vector>(NSpec);
            for(int j=0;j<NSpec;j++)
            {
                rSCoords[i][j][0] = 0.25*(
                                            rCoords[0][0]*(1.-nodes[i])*(1.-nodes[j]) + 
                                            rCoords[1][0]*(1.+nodes[i])*(1.-nodes[j]) + 
                                            rCoords[2][0]*(1.+nodes[i])*(1.+nodes[j]) + 
                                            rCoords[3][0]*(1.-nodes[i])*(1.+nodes[j])
                                            );
                
                rSCoords[i][j][1] = 0.25*(
                                            rCoords[0][1]*(1.-nodes[i])*(1.-nodes[j]) +
                                            rCoords[1][1]*(1.+nodes[i])*(1.-nodes[j]) +
                                            rCoords[2][1]*(1.+nodes[i])*(1.+nodes[j]) +
                                            rCoords[3][1]*(1.-nodes[i])*(1.+nodes[j])
                                            );
            }
        }
    }

    void SE_QLGL::generateSpectralNodesAndWeights()
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

    void SE_QLGL::generateDerivativeMatrix()
    {
        using std::vector;

        D_Matrix.resize(NSpec);
        int i,j;
        for (i=0;i<NSpec;i++)
        {
            D_Matrix[i] =  vector<double>(NSpec);
            
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


    Vector SE_QLGL::mapToLocal(const Vector& point, const std::vector< Vector >& nodes) const 
    {
        Vector E = 0.25*(-nodes[0] + nodes[1] + nodes[2] - nodes[3]);
        Vector F = 0.25*(-nodes[0] - nodes[1] + nodes[2] + nodes[3]);
        Vector G = 0.25*( nodes[0] - nodes[1] + nodes[2] - nodes[3]);
        Vector H =  point  -  0.25*(nodes[0] + nodes[1] + nodes[2] + nodes[3]);
        
        Scalar k0 = H.x()*E.y() - H.y()*E.x();
        Scalar k1 = E.x()*F.y() - E.y()*F.x() + H.x()*G.y() - H.y()*G.y();
        Scalar k2 = G.x()*F.y() - G.y()*F.x();
        
        Vector Ksi;
        
        if(std::abs(k2) < 1e-8)
        {
            Ksi.y() = - k0/k1;
            
            if(std::abs(E.x()) < 1e-8)
            {
                Ksi.x() = ( H.y() - F.y()*Ksi.y() ) / E.y();
            }
            else
                Ksi.x() = ( H.x() - F.x()*Ksi.y() ) / E.x();
        }
        else
        {
            Ksi.y() = ( -k1 - std::sqrt(k1*k1 - 4.*k0*k2) ) / (2. * k2);
            if(Ksi.y() < -1 || Ksi.y() > 1)
            {
                Ksi.y() = ( -k1 + std::sqrt(k1*k1 - 4.*k0*k2) ) / (2. * k2);
            }
            
            Ksi.x() = ( H.x() - F.x()*Ksi.y() ) / ( E.x() + G.x()*Ksi.y() );
        }
        
        return Ksi;
    }
    
    numArray< Scalar > SE_QLGL::computeInterpolationCoefficients(Scalar ksi, Scalar etha) const 
    {
        numArray<Scalar> coeffs(NSpec*NSpec,0.);
        
        
        size_t ksiNode=-1;
        size_t ethaNode =-1;
        
        for(size_t i=0; i<NSpec; ++i)
        {
            if(std::abs(nodes[i]-ksi) < 1e-8)
            {
                ksiNode=i;
            }
            if(std::abs(nodes[i]-etha) < 1e-8)
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
            
            for(size_t j=0; j<NSpec; ++j)
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
        
        return coeffs;
    }
    
    double SE_QLGL::x_ksi(int ksi, int etha, const std::vector<Vector>& rC) const
    {
        // return 1./4*(-rC[0][0]*(1.-nodes[etha])+rC[1][0]*(1.-nodes[etha])+rC[2][0]*(1.+nodes[etha])-rC[3][0]*(1.+nodes[etha]));
        
        return 0.25*(              (-rC[0].x()+rC[1].x()+rC[2].x()-rC[3].x()) 
                        + nodes[etha]*( rC[0].x()-rC[1].x()+rC[2].x()-rC[3].x()) );
    }

    double SE_QLGL::x_etha(int ksi, int etha, const std::vector<Vector>& rC) const
    {
        // return  1./4*(-realCoords[0][0]*(1.-nodes[ksi])-rC[1][0]*(1.+nodes[ksi])+rC[2][0]*(1.+nodes[ksi])+rC[3][0]*(1.-nodes[ksi]));
        
        return  0.25*(             (-rC[0].x()-rC[1].x()+rC[2].x()+rC[3].x()) 
                        + nodes[ksi]*( rC[0].x()-rC[1].x()+rC[2].x()-rC[3].x()) );
    }

    double SE_QLGL::y_ksi(int ksi, int etha, const std::vector<Vector>& rC) const
    {
        //return 1./4*(-rC[0][1]*(1.-nodes[etha])+rC[1][1]*(1.-nodes[etha])+rC[2][1]*(1.+nodes[etha])-rC[3][1]*(1.+nodes[etha]));

        return 0.25*(              (-rC[0].y()+rC[1].y()+rC[2].y()-rC[3].y()) 
                        + nodes[etha]*( rC[0].y()-rC[1].y()+rC[2].y()-rC[3].y()) );
    }

    double SE_QLGL::y_etha(int ksi, int etha, const std::vector<Vector>& rC) const
    {
        //return 1./4*(-rC[0][1]*(1.-nodes[ksi])-rC[1][1]*(1.+nodes[ksi])+rC[2][1]*(1.+nodes[ksi])+rC[3][1]*(1.-nodes[ksi]));
        
        return 0.25*(             (-rC[0].y()-rC[1].y()+rC[2].y()+rC[3].y()) 
                        + nodes[ksi]*( rC[0].y()-rC[1].y()+rC[2].y()-rC[3].y()) );
    }

    double SE_QLGL::jacobian(double xksi, double xetha, double yksi, double yetha) const
    {
        return  xksi*yetha-xetha*yksi;
    }

    double SE_QLGL::jacobian(int ksi, int etha, const std::vector<Vector>& rC) const
    {
        return jacobian(
                            x_ksi (ksi,etha, rC),
                            x_etha(ksi,etha, rC),
                            y_ksi (ksi,etha, rC),
                            y_etha(ksi,etha, rC)
                        );
    }


    Vector SE_QLGL::normal(unsigned int edge, const std::vector< Vector >& rC) const 
    {
        Vector n;
        
        switch(edge)
        {
            case 0:
                n.x() = rC[1].y()-rC[0].y();
                n.y() = rC[0].x()-rC[1].x();
                break;
            case 1:
                n.x() = rC[2].y()-rC[1].y();
                n.y() = rC[1].x()-rC[2].x();
                break;
            case 2:
                n.x() = rC[3].y()-rC[2].y();
                n.y() = rC[2].x()-rC[3].x();
                break;
            case 3:
                n.x() = rC[0].y()-rC[3].y();
                n.y() = rC[3].x()-rC[0].x();
                break;
            default:
            {
                using namespace SEM::iomanagment;
                ErrorInFunction<<"wrong edge number"<<endProgram;
                return Vector();
            }
        }
        
        n.normalize();
        return n;
    }
    
    
    std::vector<int> SE_QLGL::getEdgeLocalNodesIdInMatrix(int edgeId) const
    {
        std::vector<int> edge(NSpec);
        int i=0;
        switch (edgeId)
        {
            case 0:
                for(i=0;i<NSpec;++i)
                {
                    edge[i] = i;
                }
                break;
            case 1:
                for(i=0;i<NSpec;++i)
                {
                    edge[i] = (i+1)*NSpec-1;
                }
                break;
            case 2:
                {
                    int firstInLastRow= NSpec*(NSpec-1);
                    for(i=0;i<NSpec;++i)
                    {
                        edge[i] = firstInLastRow + (NSpec-1-i);
                    }
                }
                break;
            case 3:
                for(i=0;i<NSpec;++i)
                {
                    edge[i] = NSpec*(NSpec-1-i);
                }
                break;
            default:
                break;
        }

        return edge;
    }

    std::vector<int> SE_QLGL::getEdgeGloblaNodesId(int edgeId, const Matrix<int>::type& mask) const
    {
        std::vector<int> edge(NSpec);

        switch (edgeId)
        {
            case 0:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = mask[i][0];
                }
                break;
            case 1:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = mask[NSpec-1][i];
                }
                break;
            case 2:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = mask[NSpec-1-i][NSpec-1];
                }
                break;
            case 3:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i] = mask[0][NSpec-i-1];
                }

                break;
            default: 
                break;
        }

        return edge;
    }

    std::vector<boost::array<int,2> > SE_QLGL::getEdgeLocalNodesId(int edgeId) const
    {
        std::vector<boost::array<int,2> > edge(NSpec);
        switch (edgeId)
        {
            case 0:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i][0] = i;
                    edge[i][1] = 0;
                }
                break;
            case 1:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i][0] = NSpec-1;
                    edge[i][1] = i;
                }
                break;
            case 2:
                for(int i=0;i<NSpec;i++)
                {
                    edge[i][0] = NSpec-1-i;
                    edge[i][1] = NSpec-1;
                }
                break;
            case 3:
                for(int i=0;i<NSpec;i++)
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

    void SE_QLGL::defineMaskArraySize(Matrix<int>::type& mask)
    {
        mask.resize(NSpec);
        
        for( int i=0;i<NSpec;i++ )
        {
            mask[i] = std::vector<int>(NSpec);
        }
    }

    void SE_QLGL::allocateMaskArray_VertexNode(int vertexNr, int nodeGlobalAdress, Matrix<int>::type & addressArray)
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
            addressArray[NSpec-1][NSpec-1]=nodeGlobalAdress;
            break;
        default:
            addressArray[0][NSpec-1] = nodeGlobalAdress;
            break;
        }
    }

    void SE_QLGL::allocateMaskArray_NodesInEdge
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
            for(i=1;i<NSpec-1;i++)// from left to right, on bottom edge
            {
                addressArray[i][0]=edgeGlobalAdress[i-1];
            }
            break;
        case 1:
            for(i=1;i<NSpec-1;i++)// from bottom to top, on right edge
            {
                addressArray[NSpec-1][i]=edgeGlobalAdress[i-1];
            }
            break;
        case 2:
            for(i=1;i<NSpec-1;i++)// from right to left, on top edge
            {
                addressArray[(NSpec-1)-i][NSpec-1]=edgeGlobalAdress[i-1];
            }
            break;
        default:
            for(i=1;i<NSpec-1;i++)// from top to bottm, on left edge
            {
                addressArray[0][(NSpec-1)-i]=edgeGlobalAdress[i-1];
            }
            break;
        }
    }

    void SE_QLGL::transformMatrixMaskToVectorMask(const Matrix<int>::type &mask, std::vector<int> &vecMask) const
    {
        vecMask.clear();
        vecMask.reserve(NSpec*NSpec);

        for(int i=0;i<NSpec;++i)
            for(int j=0; j<NSpec;++j)
                vecMask.push_back(mask[j][i]);
    }

    void SE_QLGL::allocateMaskArray_InteriorNodes(int startAddress, Matrix<int>::type & addressArray)
    {
        for(int i=1;i<NSpec-1;i++)
        {
            for(int j=1;j<NSpec-1;j++)
            {
                addressArray[j][i]=startAddress;
                ++startAddress;
            }
        }
    }

    void SE_QLGL::generateH_Matrix(HMatrix& H_Matrix, const std::vector<Vector>& realCoords)
    {
        using std::vector;
        using boost::array;

        double xksi,yksi,xeta,yeta;
        double ksix,ksiy,etax,etay;
        double J;

        H_Matrix.resize(NSpec);
        for(int i=0;i<NSpec;i++)
        {
            H_Matrix[i].resize(NSpec);
            for(int j=0;j<NSpec;j++)
            {
                // pochodne rzeczywistych wsp. po elementowych:
                xksi = x_ksi (i, j, realCoords);
                xeta = x_etha(i, j, realCoords);
                yksi = y_ksi (i, j, realCoords);
                yeta = y_etha(i, j, realCoords);
                
                // obliczenie jakobianu przeksztalecenia:
                J = jacobian(xksi,xeta,yksi,yeta);
                

                // przeliczenie pochodnych na odwrotne:
                ksix = yeta/J;
                ksiy =-xeta/J;
                etax =-yksi/J;
                etay = xksi/J;
                
                                
                Scalar JW = std::abs(J)*weights[i]*weights[j];
                
                // Obliczenie elementow macierzy H-(w jawnej formie po zsumowaniu symetryczna),
                // tutaj rozlozona na 6 elementow, ze wzgledu na ewentualne pole tensorowe :
                H_Matrix[i][j][0] = (ksix*ksix)*JW;
                H_Matrix[i][j][1] = (ksix*etax)*JW;
                H_Matrix[i][j][2] = (etax*etax)*JW;

                H_Matrix[i][j][3] = (ksiy*ksiy)*JW;
                H_Matrix[i][j][4] = (ksiy*etay)*JW;
                H_Matrix[i][j][5] = (etay*etay)*JW;

            }
        }

    }

    void SE_QLGL::generateM_Matrix(numArray2D<Scalar> &M_Matrix, const std::vector<Vector> & realCoords) const
    {
        using std::vector;

        double J;

        M_Matrix.clear();
        M_Matrix.resize(NSpec);

        for(int i=0;i<NSpec;++i)
        {
            M_Matrix[i].resize(NSpec);
            for(int j=0;j<NSpec;++j)
            {
                //Modul jakobianu:
                J=std::abs( jacobian(i, j, realCoords) );

                M_Matrix[i][j] = J*weights[i]*weights[j];
            }
        }
    }

    void SE_QLGL::generateStiff_Matrix
    (
        Matrix4<double>::type & Stiff_Matrix,
        double scalarValue,
        const HMatrix & H_Matrix
    ) const
    {
        using std::vector;
        using boost::array;

        //Init matrix
        Stiff_Matrix.clear();
        Stiff_Matrix.resize(NSpec);
        
        for(int p=0;p<NSpec;p++)
        {
            Stiff_Matrix[p].resize(NSpec);
            for(int q=0;q<NSpec;q++)
            {
                Stiff_Matrix[p][q].resize(NSpec);
                for(int i=0;i<NSpec;i++)
                {
                    Stiff_Matrix[p][q][i].resize(NSpec,0.);
                }
            }
        }

        //Fill matrix
        for(int p=0;p<NSpec;p++)
        {
            for(int q=0;q<NSpec;q++)
            {
                for(int i=0;i<NSpec;i++)
                {
                    for(int j=0;j<NSpec;j++)
                    {
                        Stiff_Matrix[p][q][j][q] += scalarValue*((H_Matrix[i][q][0]+H_Matrix[i][q][3])*D_Matrix[i][j]*D_Matrix[i][p]);
                        Stiff_Matrix[p][q][i][j] += scalarValue*((H_Matrix[p][j][1]+H_Matrix[p][j][4])*D_Matrix[p][i]*D_Matrix[j][q] + (H_Matrix[i][q][1]+H_Matrix[i][q][4])*D_Matrix[q][j]*D_Matrix[i][p]);
                        Stiff_Matrix[p][q][p][i] += scalarValue*((H_Matrix[p][j][2]+H_Matrix[p][j][5])*D_Matrix[j][i]*D_Matrix[j][q]);
                    }
                }
            }
        }

    }

    numArray2D<Scalar> SE_QLGL::generateStiff_Matrix
    (
        const double &scalarValue,
        const HMatrix &H_Matrix
    ) const
    {
        //Resize matrix
        numArray2D<Scalar> matrix(NSpec*NSpec,NSpec*NSpec,0.);
        
        for(int p=0; p<NSpec;++p) //column index in element nodes
        {
            for(int q=0; q<NSpec;++q) //row index in element nodes
            {
                int row = q*NSpec+p;//row in matrix
                for(int i=0; i<NSpec; ++i) //column index in element nodes
                {
                    for(int j=0;j<NSpec; j++) //row index in element nodes
                    {
                        matrix[row][q*NSpec+j] += scalarValue*((H_Matrix[i][q][0]+H_Matrix[i][q][3])*D_Matrix[i][j]*D_Matrix[i][p]);
                        matrix[row][j*NSpec+i] += scalarValue*((H_Matrix[p][j][1]+H_Matrix[p][j][4])*D_Matrix[p][i]*D_Matrix[j][q] + (H_Matrix[i][q][1]+H_Matrix[i][q][4])*D_Matrix[q][j]*D_Matrix[i][p]);
                        matrix[row][i*NSpec+p] += scalarValue*((H_Matrix[p][j][2]+H_Matrix[p][j][5])*D_Matrix[j][i]*D_Matrix[j][q]);
                    }
                }
            }
        }
        return matrix;
    }
    
    void SE_QLGL::generateWeekVddx(const numArray< Scalar >& u, const std::vector< Vector >& rCoords, numArray< Scalar >& result) const 
    {
        result.resize(NSpec*NSpec);
        result =0.;
        
        for(size_t p=0; p<NSpec; ++p)
        {
            for(size_t q=0; q<NSpec; ++q)
            {
                size_t row = q*NSpec + p;
                
                for(size_t k=0; k<NSpec; ++k)
                {
                    Scalar J = jacobian(k,q,rCoords);
                    Scalar H11 = y_etha(k,q,rCoords)*std::abs(J)/J*weights[k]*weights[q];
                    
                    J = jacobian(p,k,rCoords);
                    Scalar H21 = -y_ksi(p,k,rCoords)*std::abs(J)/J*weights[p]*weights[k];
                    
                    result[row]+=(u[q*NSpec+k]*D_Matrix[k][p]*H11 + u[k*NSpec+p]*D_Matrix[k][q]*H21);
                }
            }
        }
    }
    
    void SE_QLGL::generateWeekVddy(const numArray< Scalar >& u, const std::vector< Vector >& rCoords, numArray< Scalar >& result) const 
    {
        result.resize(NSpec*NSpec);
        result =0.;
        
        for(size_t p=0; p<NSpec; ++p)
        {
            for(size_t q=0; q<NSpec; ++q)
            {
                size_t row = q*NSpec + p;
                
                for(size_t k=0; k<NSpec; ++k)
                {
                    Scalar J = jacobian(k,q,rCoords);
                    Scalar H12 = -x_etha(k,q,rCoords)*std::abs(J)/J*weights[k]*weights[q];

                    J = jacobian(p,k,rCoords);
                    Scalar H22 = x_ksi(p,k,rCoords)*std::abs(J)/J*weights[p]*weights[k];
                    
                    result[row] += ( u[q*NSpec+k]*D_Matrix[k][p]*H12 + u[k*NSpec+p]*D_Matrix[k][q]*H22 );
                }
            }
        }
    }
    
    void SE_QLGL::generate_ddx(const numArray< Scalar >& u, const std::vector< Vector >& rCoords, numArray< Scalar >& result) const 
    {
        result.resize(NSpec*NSpec,0.);
        
        for(size_t i=0; i<NSpec; ++i)
        {
            for(size_t j=0; j<NSpec; ++j)
            {
                Scalar sumKsi=0;
                Scalar sumEtha=0;
                
                for(size_t k=0; k<NSpec; ++k)
                {
                    sumKsi += u[j*NSpec+k]*D_Matrix[i][k];
                    sumEtha+= u[k*NSpec+i]*D_Matrix[j][k];
                } //k
                
                Scalar J = jacobian(i,j,rCoords);
                Scalar ksi_x  = y_etha(i,j,rCoords) / J;
                Scalar etha_x =-y_ksi(i,j,rCoords)  / J;
                
                result[j*NSpec+i] = sumKsi * ksi_x  +  sumEtha * etha_x;
                
            } //j 
        }// i
    }
    void SE_QLGL::generate_ddy(const numArray< Scalar >& u, const std::vector< Vector >& rCoords, numArray< Scalar >& result) const 
    {
        result.resize(NSpec*NSpec,0.);
        
        for(size_t i=0; i<NSpec; ++i)
        {
            for(size_t j=0; j<NSpec; ++j)
            {
                Scalar sumKsi=0;
                Scalar sumEtha=0;
                
                for(size_t k=0; k<NSpec; ++k)
                {
                    sumKsi += u[j*NSpec+k]*D_Matrix[i][k];
                    sumEtha+= u[k*NSpec+i]*D_Matrix[j][k];
                } //k
                
                Scalar J = jacobian(i,j,rCoords);
                Scalar ksi_y  =-x_etha(i,j,rCoords) / J;
                Scalar etha_y = x_ksi(i,j,rCoords)  / J;
                
                result[j*NSpec+i] = sumKsi * ksi_y  +  sumEtha * etha_y;
                
            } //j 
        }// i
    }
    
    
    

    void SE_QLGL::generateEdgeNeumanBCIntegralCoefficient
    (
        int edgeNumber,
        const std::vector<Vector> & eC,
        std::vector<double>& resultNodalVal
    )
    {
        resultNodalVal.resize(NSpec);
        double J1, J2,JJ;

        switch(edgeNumber)
        {
        case 0:
            J1 = ( (-eC[0][0]+eC[1][0]+eC[2][0]-eC[3][0]) + nodes[0]*(eC[0][0]-eC[1][0]+eC[2][0]-eC[3][0]) )/4;
            J2 = ( (-eC[0][1]+eC[1][1]+eC[2][1]-eC[3][1]) + nodes[0]*(eC[0][1]-eC[1][1]+eC[2][1]-eC[3][1]) )/4;
            JJ = std::sqrt(J1*J1+J2*J2);
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = JJ*weights[i];//*weights[0]
            }
            break;
        case 1:
            J1 = ( (-eC[0][0]-eC[1][0]+eC[2][0]+eC[3][0]) + nodes[NSpec-1]*(eC[0][0]-eC[1][0]+eC[2][0]-eC[3][0]) )/4;
            J2 = ( (-eC[0][1]-eC[1][1]+eC[2][1]+eC[3][1]) + nodes[NSpec-1]*(eC[0][1]-eC[1][1]+eC[2][1]-eC[3][1]) )/4;
            JJ = std::sqrt(J1*J1+J2*J2);
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = JJ*weights[i];//*weights[NSpec-1]
            }
            break;
        case 2:
            J1 = ( (-eC[0][0]+eC[1][0]+eC[2][0]-eC[3][0]) + nodes[NSpec-1]*(eC[0][0]-eC[1][0]+eC[2][0]-eC[3][0]) )/4;
            J2 = ( (-eC[0][1]+eC[1][1]+eC[2][1]-eC[3][1]) + nodes[NSpec-1]*(eC[0][1]-eC[1][1]+eC[2][1]-eC[3][1]) )/4;
            JJ = std::sqrt(J1*J1+J2*J2);
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = JJ*weights[NSpec-1-i];//*weights[NSpec-1]
            }

            break;
        case 3:
            J1 = ( (-eC[0][0]-eC[1][0]+eC[2][0]+eC[3][0]) + nodes[0]*(eC[0][0]-eC[1][0]+eC[2][0]-eC[3][0]) )/4;
            J2 = ( (-eC[0][1]-eC[1][1]+eC[2][1]+eC[3][1]) + nodes[0]*(eC[0][1]-eC[1][1]+eC[2][1]-eC[3][1]) )/4;
            JJ = std::sqrt(J1*J1+J2*J2);
            for(int i=0;i<NSpec;i++)
            {
                resultNodalVal[i] = JJ*weights[NSpec-1-i];//*weights[0]
            }
            break;
        }
    }
    
    
    Scalar SE_QLGL::rotU(const numArray< Vector >& U,size_t ksi, size_t etha, const std::vector<Vector> & rCoords) const 
    {
        Scalar sumVx=0;
        Scalar sumVy=0;
        Scalar sumUx=0;
        Scalar sumUy=0;
        
        for(size_t k=0; k<NSpec; ++k)
        {
            sumVx += U[k+etha*NSpec].y() * D_Matrix[ksi][k];
            sumVy += U[ksi+k*NSpec].y()  * D_Matrix[etha][k];
            sumUx += U[k+etha*NSpec].x() * D_Matrix[ksi][k];
            sumUy += U[ksi+k*NSpec].x()  * D_Matrix[etha][k];
        }
        
        Scalar J = jacobian(ksi,etha,rCoords);
        Scalar l11=ksi_x(ksi,etha,rCoords,J);
        Scalar l12=ksi_y(ksi,etha,rCoords,J);
        Scalar l21=etha_x(ksi,etha,rCoords,J);
        Scalar l22=etha_y(ksi,etha,rCoords,J);
        
        return (sumVx*l11 + sumVy*l21 - sumUx*l12 - sumUy*l22);
    }
    
    numArray2D< Scalar > SE_QLGL::convMatrix(const numArray< Vector >& U, const std::vector< Vector >& rC) const 
    {
        numArray2D<Scalar> matrix(NSpec*NSpec,NSpec*NSpec,0.);
        Scalar ksiCoeff,ethaCoeff,J;
        
        for(size_t i=0; i<NSpec;++i)
        {
            for(size_t j=0; j<NSpec; ++j)
            {
                size_t row =i+j*NSpec;
                J = jacobian(i,j,rC);
                ksiCoeff = ksi_x(i,j,rC,J)*U[row].x()+ksi_y(i,j,rC,J)*U[row].y();
                ethaCoeff = etha_x(i,j,rC,J)*U[row].x()+etha_y(i,j,rC,J)*U[row].y();
                
                for(size_t k=0; k<NSpec; ++k)
                {
                    matrix[row][k+j*NSpec] += D_Matrix[i][k]*ksiCoeff;
                    matrix[row][i+k*NSpec] += D_Matrix[j][k]*ethaCoeff;
                }
            }
        }
        return matrix;
    }
    
    
    Scalar SE_QLGL::gamma(std::size_t node, unsigned int edge, const std::vector< Vector >& rCoords) const 
    {
        Scalar dx,dy;
        switch(edge)
        {
            case 0:
                dx = x_ksi(node,0,rCoords);
                dy = y_ksi(node,0,rCoords);
                break;
            case 1:
                dx = x_etha(NSpec-1,node,rCoords);
                dy = y_etha(NSpec-1,node,rCoords);
                break;
            case 2:
                node = NSpec-1-node;
                dx = x_ksi(node,NSpec-1,rCoords);
                dy = y_ksi(node,NSpec-1,rCoords);
                break;
            case 3:
                node = NSpec-1-node;
                dx = x_etha(0,node,rCoords);
                dy = y_etha(0,node,rCoords);
                break;
            default:
                ErrorInFunction << "wrong edge number " <<iomanagment::endProgram;
        }
        
        return std::sqrt(dx*dx + dy*dy);
    }
    
    numArray< Scalar > SE_QLGL::boundInt_NxRotRotU(const numArray< Vector >& U, const std::vector< Vector >& rCoords, const size_t& edge) const
    {
        numArray<Scalar> integral(NSpec*NSpec,0.);
        
        Scalar nKsi, nEtha, l11, l12, l21,l22, J;
        Vector edgeNormal = normal(edge, rCoords);
        
        switch( edge )
        {
            case 0:
            {
                for(size_t p=0; p<NSpec; ++p )
                {
                    for(size_t q=0; q<NSpec; ++q)
                    {
                        if(q==0)
                        {
                            for(size_t i=0; i<NSpec; ++i)
                            {
                                J = jacobian(i,0,rCoords);
                                l11=ksi_x(i,0,rCoords,J);
                                l12=ksi_y(i,0,rCoords,J);
                                nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                                
                                integral[p+q*NSpec] += rotU(U,i,0,rCoords)*D_Matrix[i][p]*nKsi*gamma(i,0,rCoords)*weights[i];
                            }
                        }
                        
                        J = jacobian(p,0,rCoords);
                        l21=etha_x(p,0,rCoords,J);
                        l22=etha_y(p,0,rCoords,J);
                        nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                        
                        integral[p+q*NSpec] += rotU(U,p,0,rCoords)*D_Matrix[0][q]*nEtha*gamma(p,0,rCoords)*weights[p];
                    }
                }
                break;
            }
            case 1:
                for(size_t p=0; p<NSpec; ++p )
                {
                    for(size_t q=0; q<NSpec; ++q)
                    {
                        J = jacobian(NSpec-1,q,rCoords);
                        l11=ksi_x(NSpec-1,q,rCoords,J);
                        l12=ksi_y(NSpec-1,q,rCoords,J);
                        nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                        
                        integral[p+q*NSpec] += rotU(U,NSpec-1,q,rCoords)*D_Matrix[NSpec-1][p]*nKsi*gamma(q,1,rCoords)*weights[q];
                        
                        if(p==NSpec-1)
                        {
                            for(size_t j=0; j<NSpec; ++j)
                            {
                                J = jacobian(NSpec-1,j,rCoords);
                                l21=etha_x(NSpec-1,j,rCoords,J);
                                l22=etha_y(NSpec-1,j,rCoords,J);
                                nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                                
                                integral[p+q*NSpec] += rotU(U,NSpec-1,j,rCoords)*D_Matrix[j][q]*nEtha*gamma(j,1,rCoords)*weights[j];
                            }
                        }
                    }
                }
                break;
            case 2:
                for(size_t p=0; p<NSpec; ++p )
                {
                    for(size_t q=0; q<NSpec; ++q)
                    {
                        if(q==NSpec-1)
                        {
                            for(size_t i=0; i<NSpec; ++i)
                            {
                                J = jacobian(i,NSpec-1,rCoords);
                                l11=ksi_x(i,NSpec-1,rCoords,J);
                                l12=ksi_y(i,NSpec-1,rCoords,J);
                                nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                                
                                integral[p+q*NSpec] += rotU(U,i,NSpec-1,rCoords)*D_Matrix[i][p]*nKsi*gamma(NSpec-1-i,2,rCoords)*weights[i];
                            }
                        }
                        
                        J = jacobian(p,NSpec-1,rCoords);
                        l21=etha_x(p,NSpec-1,rCoords,J);
                        l22=etha_y(p,NSpec-1,rCoords,J);
                        nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                        
                        integral[p+q*NSpec] += rotU(U,p,NSpec-1,rCoords)*D_Matrix[NSpec-1][q]*nEtha*gamma(NSpec-1-p,2,rCoords)*weights[p];
                    }
                }
                break;
            case 3:
                for(size_t p=0; p<NSpec; ++p )
                {
                    for(size_t q=0; q<NSpec; ++q)
                    {
                        J = jacobian(0,q,rCoords);
                        l11=ksi_x(0,q,rCoords,J);
                        l12=ksi_y(0,q,rCoords,J);
                        nKsi = edgeNormal.x()*l12 - edgeNormal.y()*l11;
                        
                        integral[p+q*NSpec] += rotU(U,0,q,rCoords)*D_Matrix[0][p]*nKsi*gamma(NSpec-1-q,3,rCoords)*weights[q];
                        
                        if(p==0)
                        {
                            for(size_t j=0; j<NSpec; ++j)
                            {
                                J = jacobian(0,j,rCoords);
                                l21=etha_x(0,j,rCoords,J);
                                l22=etha_y(0,j,rCoords,J);
                                nEtha = edgeNormal.x()*l22 - edgeNormal.y()*l21;
                                
                                integral[p+q*NSpec] += rotU(U,0,j,rCoords)*D_Matrix[j][q]*nEtha*gamma(NSpec-1-j,3,rCoords)*weights[j];
                            }
                        }
                    }
                }
                break;
            default:
                ErrorInFunction << "wrong edge id used with quad element"<<iomanagment::endProgram;
                break;
        }
        
        return integral;
    }
    
    
    int SE_QLGL::generateNumberOfInteriorNodes()
    {
        return (NSpec-2)*(NSpec-2);
    }

    void SE_QLGL::generateLocalIndexesInMatrix(VectorToMatrixMap &map) const
    {
        map.resize(NSpec*NSpec);
        for(int i=0;i<NSpec;++i)
        {
            for(int j=0;j<NSpec;++j)
            {
                map[j*NSpec+i][0]=i;
                map[j*NSpec+i][1]=j;
            }
        }
    }

    std::vector<int> SE_QLGL::generateNumberOfInteriorInEdgesNodes()
    {
        return std::vector<int> (4,NSpec-2);
    }

    void SE_QLGL::qdqLN(double x, int N, double*q, double*dq, double*LN)
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

    double SE_QLGL::WJ(int j)
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

    std::ostream& operator<<(std::ostream& o, SE_QLGL& e)
    {
        o<<"========== Spectral Element::LGL=========="<<std::endl<<"number of nodes in element:\t"<<e.getSpectralSize()<<std::endl;
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
    
    

}//mesh
}//SEM

