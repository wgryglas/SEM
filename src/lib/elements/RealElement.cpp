#include"RealElement.h"
#include "utilities/VectorUtils.h"

namespace SEM{ namespace mesh{

    RealElement::RealElement
	(
        const std::vector<Vector>& nodes,
		SpectralElement& sEnt,
        const Matrix<int>::type &indexMask
    ): m_nodes(nodes),  m_sEnt(sEnt), m_indexMask(indexMask)
	{
        sEnt.generateH_Matrix(m_hMatrix,nodes);
        sEnt.generateM_Matrix(m_massMatrix,nodes);
        sEnt.transformMatrixMaskToVectorMask(indexMask,m_indexVectorMask);
    }

    RealElement::RealElement(RealElement & other)
        : m_nodes(other.m_nodes),
          m_sEnt(other.m_sEnt),
          m_indexMask(other.m_indexMask),
          m_hMatrix(other.m_hMatrix),
          m_massMatrix(other.m_massMatrix),
          m_indexVectorMask(other.m_indexVectorMask)
    {
    }

    RealElement & RealElement::operator=(const RealElement& e)
    {
        m_nodes = e.m_nodes;
        m_sEnt  = e.m_sEnt;

        return *this;
    }
    
    numArray2D<Scalar> RealElement::stiffMatrix(const Scalar coeff) const 
    {
        return m_sEnt.generateStiff_Matrix(coeff,m_hMatrix);
    }
    
    bool RealElement::containsPoint(const Vector& point) const 
    {
        
        if( (point-m_nodes[m_nodes.size()-1]).lengthSqrt() <1e-8 ) return true;
        
        //check if lies on some edge
        for(size_t i=0; i<m_nodes.size(); ++i)
        {
            size_t i2 = i==(m_nodes.size()-1) ? 0 : i+1;
            Vector e1 = m_nodes[i2] - m_nodes[i];
            Vector e2 = point - m_nodes[i];
                
            if(e2.lengthSqrt() < 1e-8) return true;
            
            bool moreCond = point.x() > m_nodes[i].x() && point.y() > m_nodes[i].y();
            bool lessCond = point.x() < m_nodes[i+1].x() && point.y() < m_nodes[i+1].y();
            
            if( std::abs( crossProd(e1,e2) ) < 1e-8  &&  moreCond && lessCond )
                return true;
        }
        
        //check if is inside
        bool odd=false;
        for(size_t i=0, j=m_nodes.size()-1; i<m_nodes.size(); j = i++ )
        {
            bool condition1 = ( m_nodes[i].y() > point.y() ) != ( m_nodes[j].y() > point.y() ) ;
            
            Scalar dx1 = m_nodes[j].x() - m_nodes[i].x();
            Scalar dy1 = point.y() - m_nodes[i].y();
            Scalar dy2 = m_nodes[j].y() - m_nodes[i].y();
            bool condition2 =  point.x() < ( dx1 * dy1 / dy2 + m_nodes[i].x() );
            
            if( condition1 && condition2 ) odd = !odd;
        }
        return odd;
    }
    
    void RealElement::weakGradV(const numArray< Scalar >& u, numArray< Vector >& result) const 
    {
        numArray<Scalar> ux;
        numArray<Scalar> uy;
        
        m_sEnt.generateWeekVddx(u,m_nodes,ux);
        m_sEnt.generateWeekVddy(u,m_nodes,uy);
        
        result.resize(ux.size());
        
        xCmps(result) = ux;
        yCmps(result) = uy;
    }
    
    void RealElement::weakDivW(const numArray< Vector >& u, numArray< Scalar >& result) const 
    {
        numArray<Scalar> dudx;
        numArray<Scalar> dudy;
        
        numArray<Scalar> ux(xCmps(u));
        numArray<Scalar> uy(yCmps(u));
        
        m_sEnt.generateWeekVddx(ux,m_nodes,dudx);
        m_sEnt.generateWeekVddy(uy,m_nodes,dudy);
        
        result.resize(dudx.size());
        result = dudx + dudy;  
    }
    
    void RealElement::grad(const numArray< Scalar >& u, numArray< Vector >& result) const 
    {
        numArray<Scalar> ux;
        numArray<Scalar> uy;
        
        m_sEnt.generate_ddx(u,m_nodes,ux);
        m_sEnt.generate_ddy(u,m_nodes,uy);
        
        result.resize(ux.size());
        
        xCmps(result) = ux;
        yCmps(result) = uy;
    }
    
    void RealElement::div(const numArray< Vector >& u, numArray< Scalar >& result) const 
    {
        numArray<Scalar> dudx;
        numArray<Scalar> dudy;
        
        numArray<Scalar> ux(xCmps(u));
        numArray<Scalar> uy(yCmps(u));
        
        m_sEnt.generate_ddx(ux,m_nodes,dudx);
        m_sEnt.generate_ddy(uy,m_nodes,dudy);
        
        result.resize(dudx.size());
        result = dudx + dudy;   
    }
    numArray< Scalar > RealElement::boundInt_NxRotRotU(const numArray< Vector >& U, int edge) const 
    {
        return m_sEnt.boundInt_NxRotRotU(U,m_nodes,edge);
    }
    
    numArray2D< Scalar > RealElement::convDerivMatrix(const numArray< Vector >& U) const 
    {
        return m_sEnt.convMatrix(U,m_nodes);
    }
    
    numArray< Scalar > RealElement::rot(const numArray< Vector >& U) const 
    {
        return m_sEnt.rot(U,m_nodes);
    }


}//mesh
}//SEM

