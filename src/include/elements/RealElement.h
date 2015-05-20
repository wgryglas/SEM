
#ifndef _RealElement_H_
#define _RealElement_H_

#include<iostream>

#include "boost/array.hpp"

#include "utilities/Reference.h"
#include "components/Basic.h"
#include "SpectralElement.h"

namespace SEM{ namespace mesh{

    class RealElement
    {
        REFERENCE_TYPE(RealElement)
    public:
        typedef std::vector<std::vector<boost::array<double,6> > > HMatrix;

    protected:
        int m_id;
        
        /// \var m_nodes geometri element nodes
        std::vector<Vector> m_nodes;
        
        /// \var m_sEnt Spcectral element atached to this phisical element
        SpectralElement& m_sEnt;
        
        /// \var local index map - 2D array of global indexes, representing phisical 
        /// locations of node in spectral element
        Matrix<int>::type m_indexMask;
        
        /// \var vectoral index map - global indexes located in the same order as
        /// they would in eg. local stiffnes matrix row/column
        std::vector<int> m_indexVectorMask;
        
        /// \var m_hMatrix - cached values for further computations
        HMatrix m_hMatrix;
        
        /// \var m_massMatrix -mass matrix for this phisical element
        numArray2D<double> m_massMatrix;


    public:
        
        RealElement(RealElement& other);
        
        RealElement& operator=(const RealElement& other);

        RealElement( const std::vector<Vector> &m_nodes, SpectralElement& m_sEnt, const Matrix<int>::type &m_indexMask);

        /// Mass matrix getter
        inline const numArray2D<double> & massMatrix() const {return m_massMatrix; }

        /// H matrix getter
        inline const HMatrix & hMatrix() const { return m_hMatrix; }

        /// id getter
        inline int id() const { return m_id; }

        /// nodes getter
        inline const std::vector<Vector>& nodes() const { return m_nodes; }

        bool containsPoint(const Vector & point) const;
        
        numArray<Scalar> computeInterpCoeffs(const Vector & point) const
        {
            Vector Ksi = m_sEnt.mapToLocal(point, m_nodes);
            return m_sEnt.computeInterpolationCoefficients(Ksi.x(), Ksi.y());
        }
        
        /// indexMask getter
        /// \return matrix of indexes - ordered in the same way as nodes in
        /// 2D element ([row][column], stored row by row)
        inline const Matrix<int>::type& indexMask() const { return m_indexMask; }

        inline const VectorToMatrixMap & localNodesInMatrix() const {return m_sEnt.localNodesInMatrix();}

        /// indexVectorMask - getter
        /// \return std::vector of global node indexes stored in the same
        /// order as in local matrix row/column corresponding to node
        inline const std::vector<int>& indexVectorMask() const {return m_indexVectorMask;}

        /// spectral element getter
        inline SpectralElement& spectralElement() { return m_sEnt;}

        /// spectral element const getter
        inline const SpectralElement& spectralElement() const { return m_sEnt;}

        
        numArray2D<Scalar> stiffMatrix(const Scalar coeff) const;
        
        /// \brief weekGradV - function evaluates rhs vector of equation (in week form) of 
        /// gradient, by evaluation below expression:
        /// V[s(p,q)] = [integrate(u*d(Phi_pq)/dx),
        ///              integrate(u*d(Phi_pq)/dy)]
        /// where u-scalar field, Phi_pq - test function coresponding to p,q node,
        /// and s(p,q)->node location in local vector
        void weakGradV(const numArray<Scalar>& u, numArray<Vector>& result) const;
   
        /// \brief grad - function evaluates nodal values of gradient from scalar field u
        void grad(const numArray<Scalar>& u, numArray<Vector>& result) const;
        
        
        void weakDivW(const numArray<Vector>& u, numArray<Scalar>& result) const;
        
        /// \brief div - function evaluates nodal values of divergence from vector field u
        void div(const numArray<Vector>& u, numArray<Scalar>& result) const;

        numArray<Scalar> boundInt_NxRotRotU(const numArray<Vector> & U,int edge) const;
        
        numArray2D<Scalar> convDerivMatrix(const numArray<Vector>& U) const;
        
        numArray<Scalar> rot(const numArray<Vector> & U) const;
      
      
        Vector normal(unsigned int edgeId) const
        {
            return spectralElement().normal(edgeId,m_nodes);
        }
      
    };

}//mesh
}//SEM


#endif //_RealElement_H_
