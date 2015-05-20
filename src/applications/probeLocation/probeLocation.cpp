#include <vector>
#include <memory>
#include <fstream>

#include "iomanagment/case.h"
#include "time/Time.h"
#include "mesh/Mesh.h"
#include "postprocessing/PostObject.h"
#include "iomanagment/InfoStream.h"
#include "utilities/ImplementationFactory.h"
#include "fields/GeometricField.h"
#include "iomanagment/Dictionary2.h"
#include "iomanagment/InfoStream.h"
#include "utilities/ArrayFunctions.h"

namespace SEM {
class ProbeLocation
{
    DECLARE_IMPLEMENTATION_FACTORY(ProbeLocation, const iomanagment::Dictionary&, const mesh::Mesh& ,Time &)
public:
  virtual void writeValue()=0;
  virtual void finish()=0;
};


template<typename T>
class FieldProbeLocation : public ProbeLocation
{
    field::GeometricField<T> m_field;
    numArray<Scalar> m_interpolant;
    size_t m_element;
    Time & m_time;
    std::ofstream m_output;
public:
    FieldProbeLocation(const iomanagment::Dictionary & dict, const mesh::Mesh& m, Time &t)
    : m_field(dict.entry("field").value(),m,iomanagment::ALWAYS,iomanagment::NO_WRITE),
      m_element(0),
      m_time(t)
    {
        Vector p;
        dict.entry("point") >> p;
        
        for(size_t e=0; e<m.size(); ++e)
        {
            if( m[e].containsPoint(p) )
            {
                m_interpolant.resize(m[e].indexVectorMask().size());
                m_interpolant=m[e].computeInterpCoeffs(p);
                m_element= e;
            }
        }
        
        if(dict.hasEntry("file"))
        {
            m_output.open((Case::path()/dict.entry("file").getValue<std::string>()).string());
        }
        else
        {
            m_output.open( (Case::path()/dict.name()).string() );
        }
    }
    
    void writeValue() 
    {
        if(m_interpolant.size()==0)
        {
            iomanagment::Warning 
            << "location is not found, can't write location value of field "
            << m_field.objectName()<<iomanagment::endInfo;
            return;
        }
        
        m_output << m_time.timeName() << '\t';
        
        for(size_t d=0; d<CmpTraits<T>::dim(); ++d)
        {
            m_output << array::sum(m_interpolant*CmpTraits<T>::cmpArray(m_field.element(m_element),d) );
            m_output << '\t';
        }
        
        m_output <<'\n';
    }
    
    void finish()
    {
        m_output.close();
    }
    
};

}//SEM



int main(int argc, char* args[])
{
    using namespace SEM;
    using namespace SEM::iomanagment;
    using namespace SEM::mesh;
    
    //Setup case
    Case::setup(argc, args);
    
    //Create mesh
    Mesh mesh(Case::meshPath(),Case::elementsPath());
    
    // Get times in working dir
    std::vector<double> times =Case::timesInWorkingDir();
    
    if(times.size()==0)
    {
        ErrorInFunction<<"no time results in working directory"<<iomanagment::endProgram;
    }
    
    Case::time().setCurrentTime(times[0]);
    
    iomanagment::Dictionary probeDict;
    Case::probeLocationDefinition() >> probeDict;

    //Create probe locations for apropriate fields type
    typedef std::auto_ptr<ProbeLocation> ProbePtr;
    std::vector<ProbePtr > probes;
    for(auto objDictPair : probeDict.subDictionaries())
    {
        Dictionary dict;
        Case::time().localPath()->path()/(objDictPair.second->entry("field").value()+".sem") >> dict;
        
        probes.push_back( ProbePtr( ProbeLocation::Impl(dict.entry("fieldType").value() )( *objDictPair.second, mesh, Case::time() ) ) );
    }

    //write values
    for(double t : times)
    {
        Case::time().setCurrentTime(t);
        
        for( ProbePtr & probe : probes)
            probe->writeValue();
    }
    
    //finish writting
    for(ProbePtr &probe : probes)
        probe->finish();
    
}

namespace SEM {

DEFINE_IMPLEMENTATION_FACTORY(ProbeLocation,const iomanagment::Dictionary&,const mesh::Mesh&, Time&)
typedef FieldProbeLocation<Scalar> ScalarProbeLocation;
REGISTER_IMPLEMENATION(ProbeLocation, ScalarProbeLocation,"scalar")
typedef FieldProbeLocation<Vector> VectorProbeLocation;
REGISTER_IMPLEMENATION(ProbeLocation,VectorProbeLocation,"vector")

}//SEM