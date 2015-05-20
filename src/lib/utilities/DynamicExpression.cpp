#include "DynamicExpression.h"
namespace SEM {
    
    DynamicExpression::DynamicExpression() 
    : exprtk::expression<Scalar>(),m_x(0), m_y(0), m_t(0), m_transient(false), m_spatial(false)
    {
        configure();
    }
    
    DynamicExpression::DynamicExpression(const std::string& exprStr) 
    : exprtk::expression<Scalar>(), m_x(0), m_y(0), m_t(0), m_transient(false), m_spatial(false)
    {
        configure();
        setExpression(exprStr);
    }
    
    void DynamicExpression::configure() 
    {
        exprtk::symbol_table<Scalar> symbol_table;
        symbol_table.add_variable("x",m_x);
        symbol_table.add_variable("y",m_y);
        symbol_table.add_variable("t",m_t);
        symbol_table.add_constants();
        register_symbol_table(symbol_table);
    }
    void DynamicExpression::setExpression(const std::string& expr) 
    {
        exprtk::parser<Scalar> parser;
        
        if(! parser.compile(expr,*this) )
        {
            ErrorInFunction << "Wrong expression: \n "<<parser.error()<<iomanagment::endProgram;
        }
        
        updateType(expr);
    }
    
    bool DynamicExpression::isTransient() const 
    {
        return m_transient;
    }
    
    bool DynamicExpression::isSpatial() const 
    {
        return m_spatial;
    }
    
    void DynamicExpression::updateType(const std::string& exprStr) 
    {
        if( exprStr.find('x')!=std::string::npos || exprStr.find('y')!=std::string::npos )
            m_spatial = true;
        
        if(exprStr.find('t') != std::string::npos)
            m_transient = true;                
    }

}//SEM