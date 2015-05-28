#ifndef BASIC_H
#define BASIC_H

#include <iostream>
#include <vector>
#include <cmath>

#include "Eigen/Dense"

#include <boost/container/string.hpp>

#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/InfoStream.h"

namespace SEM
{
    //////////////////////////////////////////////////
    /// Define global parameter for space definition
    //////////////////////////////////////////////////
    const int SPACE_DIM = 2;

    //////////////////////////////////////////////////
    /// \typedef Scalar
    /// ---------------------------------------------
    /// Define scalar compoentent for calculation
    //////////////////////////////////////////////////
    typedef double Scalar;

    template<typename T, size_t Dim>
    struct TensorTX
    {
        static const unsigned int N = Dim*Dim;
    private:
        T data[N];
    public:
        
        TensorTX(const T& initVal=0)
        {
            for(int i=0; i<N; ++i)
                data[i]=initVal;
        }
        
        TensorTX(const TensorTX<T,Dim> &other)
        {
            for(int i=0; i<N; ++i)
                data[i]=other[i];
        }
        
        T& operator[](size_t index) {return data[index];}
        const T& operator[](size_t index) const {return data[index];}
        
        T& operator()(size_t i, size_t j) {return data[i*Dim+j]; }
        const T& operator()(size_t i, size_t j) const { return data[i*Dim+j]; }
        
    };
    
    /////////////////////////////////////////////////
    /// macro for creation operation between tensor
    /// and scalar value
    /////////////////////////////////////////////////
    #define TENSOR_SCALAR_OPERATOR(SIGN)                                             \
    template<typename T, size_t D>                                                   \
    TensorTX<T,D> operator SIGN (const T &lhs, const TensorTX<T,D> &rhs)             \
    {                                                                                \
        TensorTX<T,D> res;                                                           \
        for(int i=0; i<TensorTX<T,D>::N; i++)                                        \
            res[i] = lhs SIGN rhs[i];                                                \
        return res;                                                                  \
    }                                                                                \
    template<typename T, size_t D>                                                   \
    TensorTX<T,D> operator SIGN (const TensorTX<T,D> &lhs, const T &rhs)             \
    {                                                                                \
        TensorTX<T,D> res;                                                           \
        for(int i=0; i<TensorTX<T,D>::N; i++)                                        \
            res[i] = lhs[i] SIGN rhs;                                                \
            return res;                                                              \
    }                                                                                    
    
    TENSOR_SCALAR_OPERATOR(+)
    TENSOR_SCALAR_OPERATOR(-)
    TENSOR_SCALAR_OPERATOR(/)
    TENSOR_SCALAR_OPERATOR(*)
    
    typedef TensorTX<Scalar,SPACE_DIM> Tensor;
    
    
    
    
    
    
    //////////////////////////////////////////////////
    ///  \class Vector
    /// ---------------------------------------------
    /// Define vector component value
    /// note: inheritance is due to preparing
    /// interface. Any change of linear alg.
    /// library will require redefining Vector to allow
    /// compilation
    //////////////////////////////////////////////////
#define MAKE_VEC_UNARY_OPERATOR(OP)                           \
    VectorTX<T,N> & operator OP (const VectorTX<T,N> & other) \
    {                                                         \
        for(int i=0; i<size();i++)                            \
        {                                                     \
            data[i] OP other[i];                              \
        }                                                     \
        return *this;                                         \
    }                                                         \

    //////////////////////////////////////////////
    template<typename T,size_t N>
    class VectorTX
    {
        T data[N];

    public:
        VectorTX(const T& initVal=0)
        {
            for(size_t i=0; i<N; ++i)
                data[i]=initVal;
        }
        
        VectorTX(const T& x, const T& y)
        {
            data[0] = x;
            data[1] = y;
            
            for(size_t i=2; i<N; ++i)
                data[i] = 0;
        }

        VectorTX(const T& x, const T& y, const T& z)
        {
            data[0] = x;
            data[1] = y;
            data[2] = z;
            
            for(size_t i=3; i<N; ++i)
                data[i] = 0;
        }
        
        VectorTX(const VectorTX<T,N>& other)
        {
            for(int i=0;i<N;i++)
                data[i]=other[i];
        }

        VectorTX operator =(const VectorTX<T,N> &other)
        {
            for(int i=0; i<N;i++)
                data[i]=other[i];
        }
        
        const T& x() const { return data[0];}
        const T& y() const { return data[1];}
        const T& z() const { return data[2];}

        T& x() { return data[0];}
        T& y() { return data[1];}
        T& z() { return data[2];}
        
        T length() const
        {
            return std::sqrt(lengthSqrt());
        }
        
        T mag() const
        {
            return std::sqrt(lengthSqrt());
        }
        
        T lengthSqrt() const
        {
            T res = 0;
            for(size_t i=0; i<N; ++i)
            {
                res += data[i]*data[i];
            }
            
            return res;
        }
        
        void normalize()
        {
            T len = length();
            
            for(size_t i=0; i<N; ++i)
            {
                 data[i]/= len;
            }
        }
        
        size_t size() const
        {
            return N;
        }

        const T & operator [](size_t loc) const
        {
            return data[loc];
        }

        T & operator [](size_t loc)
        {
            return data[loc];
        }

        bool operator ==(const VectorTX<T,N>& other) const
        {
            for(int i=0;i<N;++i)
            {
                if(data[i]!=other.data[i])
                    return false;
            }
            return true;
        }

        bool operator !=(const VectorTX<T,N>& other) const
        {
            return !( operator ==(other) );
        }
        
        VectorTX<T,N> operator -() const
        {
            VectorTX<T,N> result;
            for(size_t i=0; i<N; ++i) result[i]=-(*this)[i];
            return result;
        }
        VectorTX<T,N> operator +() const
        {
            return *this;
        }
        
        bool operator<(const VectorTX<T,N>& other) const
        {
              return this->lengthSqrt() < other.lengthSqrt();
        }
        
        bool operator>(const VectorTX<T,N>& other) const
        {
              return this->lengthSqrt() > other.lengthSqrt();
        }
        

        MAKE_VEC_UNARY_OPERATOR(+=)
        MAKE_VEC_UNARY_OPERATOR(-=)
        MAKE_VEC_UNARY_OPERATOR(*=)
        MAKE_VEC_UNARY_OPERATOR(/=)

        static VectorTX<T,N> ZERO()
        {
            return VectorTX<T,N>(0);
        }

        static VectorTX<T,N> UNIT()
        {
            return VectorTX<T,N>(1);
        }
    };

    
    
    ////////////////////////////////////////////////
    /// macro for creation simple vector operators
    /// note: default operators are not mathematical
    /// exp. but simple elemnts operations
    ////////////////////////////////////////////////
    #define VEC_BINARY_OPERATOR(SIGN)                                                    \
        template<typename T, size_t N>                                                   \
        VectorTX<T,N> operator SIGN (const VectorTX<T,N> &lhs, const VectorTX<T,N> &rhs) \
        {                                                                                \
            VectorTX<T,N> res;                                                           \
            for(int i=0; i<N; i++)                                                       \
                res[i] = lhs[i] SIGN rhs[i];                                             \
            return res;                                                                  \
        }                                                                                \

    ////////////////////////////////////////////////
    /// Defined binary operators
    ///////////////////////////////////////////////
    VEC_BINARY_OPERATOR(+)
    VEC_BINARY_OPERATOR(-)
    VEC_BINARY_OPERATOR(*)
    VEC_BINARY_OPERATOR(/)

    
    /////////////////////////////////////////////////
    /// macro for creation operation between vector
    /// and scalar value
    /////////////////////////////////////////////////
    #define VEC_SCALAR_OPERATOR(SIGN)                                                \
    template<typename T, typename E, size_t N>                                       \
    VectorTX<T,N> operator SIGN (const E &lhs, const VectorTX<T,N> &rhs)             \
    {                                                                                \
        VectorTX<T,N> res;                                                           \
        for(int i=0; i<N; i++)                                                       \
            res[i] = lhs SIGN rhs[i];                                                \
        return res;                                                                  \
    }                                                                                \
    template<typename T, typename E, size_t N>                                       \
    VectorTX<T,N> operator SIGN (const VectorTX<T,N> &lhs, const E &rhs)             \
    {                                                                                \
        VectorTX<T,N> res;                                                           \
        for(int i=0; i<N; i++)                                                       \
            res[i] = lhs[i] SIGN rhs;                                                \
        return res;                                                                  \
    }                                                                                \

    ////////////////////////////////////
    /// define vector-scalar operators
    ////////////////////////////////////
    VEC_SCALAR_OPERATOR(+)
    VEC_SCALAR_OPERATOR(-)
    VEC_SCALAR_OPERATOR(/)
    VEC_SCALAR_OPERATOR(*)
    
    
    //////////////////////////////////////////////////
    /// Define function to operate on Vector 
    //////////////////////////////////////////////////
    template<typename T, size_t Dim>
    T dotProd(const VectorTX<T,Dim> & u, const VectorTX<T,Dim> & v)
    {
        T result =0;
        for(int i=0; i<Dim; ++i)
            result += u[i] * v[i];
        
        return result;
    }
    
    template<typename T>
    T crossProd(const VectorTX<T,2> & u, const VectorTX<T,2> & v)
    {
        return u[0]*v[1] - u[1]*v[0];
    }
    
    template<typename T>
    VectorTX<T,3> crossProd(const VectorTX<T,3> & u, const VectorTX<T,3> & v)
    {
        return VectorTX<T,3>(u[1]*v[2]-u[2]*v[1], u[2]*v[0]-u[0]*v[2], u[0]*v[1] - u[1]*v[0]);
    }
    
    
    //////////////////////////////////////////////////
    /// \typedef Vector
    /// ----------------------------------------------
    /// Define 'Vector' component for calculation
    //////////////////////////////////////////////////
    typedef VectorTX<Scalar,SPACE_DIM> Vector;


    //////////////////////////////////////////////////////
    /// Acctually it shall be removed, because
    /// own vector impl. is provided and set as 'Vector'
    /// type
    //////////////////////////////////////////////////////
    typedef Eigen::Matrix<Scalar,SPACE_DIM,1> VectorBase;
    class EigenVector : public VectorBase
    {
    public:
        EigenVector(Scalar x=0, Scalar y=0);
        EigenVector(const EigenVector& v);
        EigenVector & operator =(const EigenVector & v);
    };
    std::ostream & operator <<(std::ostream & stream, const EigenVector & v);
    std::istream & operator >> (std::istream & stream, EigenVector & v);

} //SEM


#endif // BASIC_H
