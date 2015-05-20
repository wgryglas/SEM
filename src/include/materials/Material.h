
#ifndef MATERIAL_H
#define MATERIAL_H

#include <map>

#include "boost/filesystem/path.hpp"

#include "iomanagment/RegistryObject.h"
#include "components/Basic.h"

namespace SEM {namespace materials {

class Material : public iomanagment::RegistryObject
{
    boost::filesystem::path m_dataBasePath;
    
    typedef std::map<std::string,Scalar> propMap;
    propMap properties;
    
public:
    Material(iomanagment::RegistryFile::ref casFile, const boost::filesystem::path &dataBase);
    
    Material(const Material &other);
    
    Material& operator=(const Material &other);
    
    inline Scalar operator[](const std::string & propName) const;
    
    virtual ~Material();
    
    /// \brief implementation of RegistryObject interfec.
    /// below method called when file is registrated, and
    /// if RegistryFile would decide-->here it should 
    /// be called again only if local mateiral file would
    /// change.
    virtual void read(const iomanagment::Dictionary& dict);
    virtual void write(iomanagment::Dictionary& dict);

};


    Scalar Material::operator[](const std::string& propName) const
    {
        using namespace iomanagment;
        
        if( properties.find(propName) != properties.end() )
        {
            return properties.at(propName);
        }
        else
        {
            ErrorInFunction<<"Property "<<propName<<" was not declared \n"
                           <<"neither in data base material definition nor\n "
                           <<"in local case material definition file. "
                           <<"If this property was not found, then consider\n"
                           <<"updating your material in data base or changing\n"
                           <<"source code, because those kinde of files are \n"
                           <<"suppouse to lower user work and setup some commmon\n"
                           <<"properties automaticly so if data base would be \n"
                           <<"correctly defined then this error shall not appear."
                           <<endProgram;
                           
        }
        
    }




}//materials
}//SEM

#endif // MATERIAL_H
