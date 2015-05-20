
#ifndef DYNAMICEXPRESSION_H
#define DYNAMICEXPRESSION_H

#include <string>

#include "components/Basic.h"
#include "exprtk/exprtk.hpp"
#include "boost/array.hpp"

namespace SEM {
    
    struct DynamicExpression : public exprtk::expression<Scalar>
    {
        DynamicExpression();
        
        DynamicExpression(const std::string &exprStr);
        
        void setX(const Scalar & x) { m_x = x;}
        void setY(const Scalar & y) { m_y = y;}
        void setT(const Scalar & t) { m_t = t;}
        
        void setExpression(const std::string &expr);
        
        bool isTransient() const;
        
        bool isSpatial() const;
        
    private:
        Scalar m_x;
        Scalar m_y;
        Scalar m_t;
        
        bool m_spatial;
        bool m_transient;
        
        void configure();
        void updateType(const std::string & exprStr);
        
    };
    
    template<size_t D>
    inline boost::array<DynamicExpression,D> build2DExpressions(const boost::array<std::string,D> & exprStrings)
    {
        boost::array<DynamicExpression,D> expressions;
        
        for(size_t i=0; i<D; ++i)
            expressions[i].setExpression(exprStrings[i]);
        
        return expressions;
    }
    
} //SEM







#endif // DYNAMICEXPRESSION_H
