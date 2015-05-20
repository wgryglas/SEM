
#include <memory>
#include <map>
#include <fstream>

#include "boost/filesystem/path.hpp"
#include "boost/filesystem/fstream.hpp"
#include "iomanagment/case.h"
#include "postprocessing/PostObject.h"
#include "utilities/numArray.h"
#include "fields/GeometricField.h"
#include "utilities/ArrayFunctions.h"
#include "time/Time.h"


namespace SEM {
    
    template<typename T>
    class MinMaxCalculatorBase : public PostObject
    {
    protected:
      field::GeometricField<T> m_field;
      boost::filesystem::path m_filePath;
      Time & m_time;
    public:
        MinMaxCalculatorBase(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
        : m_field(dict.entry("field").value(),mesh, iomanagment::ALWAYS,iomanagment::AUTO),
          m_filePath( Case::path()/dict.entry("file").value() ),
          m_time(time)
        {
        }
    };
    
    class ScalarMinMaxCalculator : public MinMaxCalculatorBase<Scalar>
    {
    public:
        ScalarMinMaxCalculator(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
        : MinMaxCalculatorBase<Scalar>(dict,mesh,time)
        {
        }
        
        void calculate()
        {
            boost::filesystem::ofstream out(m_filePath,ios::app);
            out << m_time.time() << '\t' << SEM::array::min(m_field)<<'\t' << SEM::array::max(m_field)<<std::endl;
            out.close();
        }
    };
    
    class VectorMinMaxCalculator : public MinMaxCalculatorBase<Vector>
    {
    public:
        VectorMinMaxCalculator(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
        : MinMaxCalculatorBase<Vector>(dict,mesh,time)
        {
        }
        
        void calculate()
        {
            boost::filesystem::ofstream out(m_filePath,ios::app);
            out << m_time.time() << '\t' << SEM::array::min(SEM::array::mag(m_field))<<'\t' << SEM::array::max(SEM::array::mag(m_field))<<std::endl;
            out.close();
        }
    };
    
    REGISTER_IMPLEMENATION(PostObject,ScalarMinMaxCalculator,"scalarMinMax")
    REGISTER_IMPLEMENATION(PostObject,VectorMinMaxCalculator,"vectorMinMax")
    
}