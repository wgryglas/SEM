#include "boost/assign/list_of.hpp"
#include "materials/Property.h"


namespace SEM { namespace materials {

typename Property::typesMap Property::Type = boost::assign::map_list_of("ScalarValue", ScalarValue)
                                                                       ("ScalarField", ScalarField)
                                                                       ("VectorValue", VectorValue)
                                                                       ("VectorFiled", VectorField)
                                                                       ("TensorValue", TensorValue);


Property::Property(const std::string &name, const PropertyType &type, const iomanagment::DictEntry& valueFromFile)
: m_name(name), m_type(type)
{
    using namespace iomanagment;
    
    switch(type)
    {
        case ScalarValue:
            valueFromFile>>m_scalarValue; 
            break;
        case ScalarField:
            valueFromFile>>m_scalarField;
            break;
        default:
            ErrorInFunction<<" Not yet supported Property type named "<<typeToString(type)
                           <<"for property with name "<<name<<"." 
                           <<"Change definition of property in apropriate file"<<endProgram;
    }
}




Property::Property(const Property& other)
: m_name(other.m_name), m_type(other.m_type), m_scalarValue(other.m_scalarValue)
{
}


Property::Property(const std::string &name,  const Scalar& valueToSet )
: m_name(name), m_type(ScalarValue), m_scalarValue(valueToSet)
{
}


Property::~Property()
{

}

Property& Property::operator=(const Property& other)
{
    m_type = other.m_type;
    m_name = other.m_name;
    
    using namespace iomanagment;
    switch(m_type)
    {
        case ScalarValue:
            m_scalarValue = other.m_scalarValue;
            break;
        case ScalarField:
            m_scalarField = other.m_scalarField;
            break;
        default:
            ErrorInFunction<<" Not yet supported Property type named "<<typeToString(m_type)<<"\n"
            <<"for property with name "<<m_name<<".\n" 
            <<"Error in assigment operator, possible in some place was used\n"
            <<"not supported property"<<endProgram;
    }
    
    return *this;
}


std::string Property::typeToString(const PropertyType &type)
{
    for(typename typesMap::const_iterator it=Property::Type.begin(); it!=Property::Type.end(); ++it)
    {
        if(it->second == type)
        {
            return it->first;
        }
    }
}


}//materials
}//SEM