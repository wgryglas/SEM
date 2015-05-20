#include "Material.h"

#include "iomanagment/Dictionary2.h"

namespace SEM {namespace materials {

Material::Material(iomanagment::RegistryFile::ref caseFile, const boost::filesystem::path& dataBase)
:   RegistryObject("Material"),
    m_dataBasePath(dataBase)
{
    setRegistryFile(caseFile);
}
    
Material& Material::operator=(const Material& other)
{
    m_dataBasePath = other.m_dataBasePath;
    setRegistryFile(other.file());
}


Material::Material(const Material& other)
: RegistryObject("Material"), m_dataBasePath(other.m_dataBasePath)
{
}

Material::~Material()
{
}

void Material::read(const iomanagment::Dictionary& dict)
{
    //Read dictionary from dataBase from name specified by dict
    boost::filesystem::path dbFile = m_dataBasePath/(dict.entry("material").value()+".sem");
    
    iomanagment::Dictionary dbDict;
    dbFile>>dbDict;
    
    //Override elements introduced in case local material definition
    typedef iomanagment::Dictionary::entryMap::value_type entyrPair;
//    typedef iomanagment::Dictionary::entryMap::const_iterator entryItr;
    
    if(dict.hasSubDictionary("override"))
    {
        for(const entyrPair &entry : dict.subDictionary("override").entries() )
        {
            //note that current impl. store all dictionary elements as raw ptrs 
            // it should be done in better way
            dbDict.add(new iomanagment::DictEntry(entry.first, entry.second->value()) );
        }
    }
    
    properties.clear();
    Scalar tmp;
    for(const entyrPair &entry : dbDict.entries())
    {
        *(entry.second)>>tmp;
        properties.insert(propMap::value_type(entry.first,tmp));
    }
    
}


void Material::write(iomanagment::Dictionary& dict)
{
}


}//materials 
}//SEM