#include "RegistryFile.h"

#include "RegistryObject.h"
#include "Register.h"

#include "iomanagment/InfoStream.h"
#include "iomanagment/case.h"

namespace SEM { namespace iomanagment {

RegistryFile::RegistryFile(Register &reg, const std::string &objectName, const boost::filesystem::path &local, const READ &readOpt, const WRITE &writeOpt)
    : m_register(reg),
      m_objectName(objectName),
      m_localPath(new LocalPath(local)),
      m_read(readOpt),
      m_write(writeOpt),
      currentDict(objectName)
{
    initFile();
}

RegistryFile::RegistryFile(Register &reg, const std::string &objectName, LocalPath::ref local, const READ &readOpt, const WRITE &writeOpt)
    : m_register(reg),
      m_objectName(objectName),
      m_localPath(local),
      m_read(readOpt),
      m_write(writeOpt),
      currentDict(objectName)
{
    initFile();
}

void RegistryFile::initFile()
{
    // Register this object
    m_register.registerFile(this);

    // Read file into dictionary if it's allowed
    boost::filesystem::path fullPath = filePath();
    bool fileExist=boost::filesystem::exists(fullPath);

    if(fileExist)
        lastModificationTime=boost::filesystem::last_write_time(fullPath);


    bool read=true;
    switch(m_read)
    {
    case NO_READ:
        read = false;
        break;
    case READ_IF_PRESENT:
        read = fileExist;
        break;
    default:
        if( !fileExist )
            ErrorInFunction<<"File "<<fullPath<<" don't exist, and this file must be read"<<endProgram;
        break;
    }

    if(read)
    {
        filePath()>>currentDict;
    }
}

// Disallowed assigment
RegistryFile &RegistryFile::operator =(const RegistryFile &other)
{
}

RegistryFile::~RegistryFile()
{
    m_register.unregisterFile(this);
}

void RegistryFile::registerObject(RegistryObject *regObj)
{
    m_objects.push_back(regObj);
    
    if(allowedRead() || m_read == READ_ONCE || m_read == READ_IF_MODIFIED || (exist() && m_read==READ_IF_PRESENT)) 
        regObj->read(currentDict);
}

void RegistryFile::unregisterObject(RegistryObject *regObj)
{
    std::vector<RegistryObject*>::iterator itr= std::find(m_objects.begin(),m_objects.end(),regObj);
    if(itr!=m_objects.end())
    {
        m_objects.erase(itr);
    }
}

boost::filesystem::path RegistryFile::filePath() const
{
    boost::filesystem::path p= Case::path()/m_localPath->path()/(m_objectName+".sem");
    return p;
}

void RegistryFile::fireReading()
{
    if(allowedRead())
    {
        currentDict.clear();
        filePath()>>currentDict;

        for(RegistryObject* obj: m_objects)
            obj->read(currentDict);
    }
}

void RegistryFile::fireWriting()
{
    if(allowdWrite())
    {
        currentDict.clear();

        for(RegistryObject* obj : m_objects)
            obj->write(currentDict);

        filePath()<<currentDict;
    }

}

bool RegistryFile::exist() const
{
    return boost::filesystem::exists(filePath());
}


bool RegistryFile::allowedRead()
{
    if(m_read == ALWAYS)
        return true;
    
    if(m_read==READ_IF_MODIFIED)
    {
        boost::filesystem::path fullPath = filePath();
        if( boost::filesystem::exists(fullPath) )
        {
            std::time_t modifTime = boost::filesystem::last_write_time(fullPath);

            if(modifTime!=lastModificationTime)
            {
                lastModificationTime=modifTime;
                return true;
            }
        }
        else
           ErrorInFunction<<"File "<<fullPath<<" don't exist, and this file must be read if was modified"<<endProgram;
    }

    return false;
}

bool RegistryFile::allowdWrite()
{
    return m_write == AUTO;
}
RegistryObject* RegistryFile::object(const std::string& objName) const 
    {
        for(auto itr=m_objects.begin(); itr!=m_objects.end(); ++itr )
        {
            if((*itr)->objectName()==objName)
                return *itr;
        }
        return nullptr;
    }




}//iomanagment
}//SEM
