#include "RegistryObject.h"

#include "Register.h"

namespace SEM { namespace iomanagment {

RegistryObject::RegistryObject(const std::string &name) 
: m_objectName(name)
{
}
    
void RegistryObject::setRegistryFile(RegistryFile::ref regFile)
{
    if(m_file)
    {
        m_file->unregisterObject(this);
    }

    m_file=regFile;
    m_file->registerObject(this);
}

RegistryObject::~RegistryObject()
{
    if(m_file)
        m_file->unregisterObject(this);
    
}

}//iomanagment
}//SEM
