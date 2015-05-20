
#ifndef VALUESONLINE_H
#define VALUESONLINE_H

#include <memory>
#include <map>
#include <fstream>
#include <math.h>

#include "boost/filesystem/path.hpp"

#include "postprocessing/PostObject.h"
#include "utilities/numArray.h"
#include "fields/GeometricField.h"
#include "utilities/ArrayFunctions.h"
#include "time/Time.h"


namespace SEM {
    
template<typename T>
struct ValueWriter
{
    virtual void  write(const numArray<Vector> & points, const numArray<T> & values, const boost::filesystem::path & file) const=0;
};

template<typename T>
struct DefaultWrite : public ValueWriter<T>
{
    void  write(const numArray<Vector> & points, const numArray<T> & values, const boost::filesystem::path & file) const
    {
        iomanagment::DictEntry * pointsEntry = new iomanagment::DictEntry("points");
        (*pointsEntry)<<points;
        
        iomanagment::DictEntry *valuesEntry = new iomanagment::DictEntry("values");
        
        (*valuesEntry) << values;
        
        iomanagment::Dictionary dict(  file.filename().string() );
        dict.add(pointsEntry);
        dict.add(valuesEntry);
        
        file << dict;
    }
};

template<typename T>
struct ColumnWrite : public ValueWriter<T>
{
    void  write(const numArray<Vector> & points, const numArray<T> & values, const boost::filesystem::path & file) const
    {
        std::ofstream out(file.c_str());
        
        size_t dimSize = CmpTraits<T>::dim();
        
        for(size_t i=0; i<points.size(); ++i)
        {
            out << points[i].x() << "\t" << points[i].y()<<"\t";
            
            for(size_t d=0; d<dimSize; ++d)
            {
                out << CmpTraits<T>::component(values[i],d);
                if(d!=dimSize-1)
                    out<<"\t";
                else
                    out<<std::endl;
            }
        }
        
        out.close();
    }
};

    
template<typename T>    
class ValuesOnLine :  public PostObject
{
    Time & m_time;
    const mesh::Mesh & m_mesh;
    
    std::auto_ptr<ValueWriter<T> > m_writer;
    std::string m_fileName;
    numArray<Vector> m_points;
    
    std::vector<std::pair<size_t,numArray<Scalar> > > m_interpolants;
    
    field::GeometricField<T> m_field;
    
public:
    ValuesOnLine(const iomanagment::Dictionary& dict, const mesh::Mesh& mesh, Time& time) 
    : m_time(time), m_mesh(mesh), m_field(dict.entry("field").value(), mesh, iomanagment::ALWAYS, iomanagment::NO_WRITE)
    {
        if(dict.hasEntry("writeType"))
        {
            if(dict.entry("writeType").value() == "columns")
                m_writer = std::auto_ptr< ValueWriter<T> >(new ColumnWrite<T>);
            else
                m_writer = std::auto_ptr< ValueWriter<T> >(new DefaultWrite<T>);
        }
        else
            m_writer = std::auto_ptr< ValueWriter<T> >(new DefaultWrite<T>);
        
        
        m_fileName = dict.entry("fileName").value();
        
        
        Vector start, end;
        
        dict.entry("start") >> start;
        dict.entry("end") >> end;
        
        size_t n_points;
        
        dict.entry("resolution") >> n_points;
        
        Vector dV = (end - start)/( (Scalar)(n_points-1) );
        
        m_points.resize(n_points);
        m_points[0]=start;
        
        for(size_t i=1; i<n_points; ++i)
            m_points[i] = m_points[i-1] + dV;
        
        m_interpolants.resize(n_points);
        for(size_t i=0; i<n_points; ++i)
        {
            size_t element=0;
            bool found=false;
            for(size_t e=0; e<mesh.size(); ++e)
            {
                if(mesh[e].containsPoint(m_points[i]))
                {
                    found = true;
                    element =e;
                    break;
                }
            }
            
            if(found)
            {
                m_interpolants[i] = std::make_pair( element,mesh[element].computeInterpCoeffs( m_points[i] ) );
            }
        }
    }
    
    void calculate()
    {
        numArray<T> values(m_points.size());
        size_t dim = CmpTraits<T>::dim();
        
        for(size_t i=0; i<m_interpolants.size(); ++i)
        {
            if(m_interpolants[i].second.size()>0)
            {
                auto map= m_mesh[m_interpolants[i].first].indexVectorMask();
                for(size_t d=0; d<dim; ++d)
                {
                    auto cmpField = CmpTraits<T>::cmpArray(m_field,d);
                    auto cmpFieldSlice = cmpField.slice(map);
                    CmpTraits<T>::component(values[i],d) = array::sum( m_interpolants[i].second * cmpFieldSlice );
                }
            }
            else
            {
                std::cout<<"point "<<m_points[i].x()<<", "<<m_points[i].y()<<" not found in element"<<std::endl;
            }
        }
        
        
        boost::filesystem::path filePath= m_time.localPath()->path()/(m_fileName);
        
        m_writer->write(m_points,values,filePath);
    }
    
};

}//SEM


#endif // VALUESONLINE_H
