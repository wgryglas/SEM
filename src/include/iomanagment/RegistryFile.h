#ifndef REGISTERYFILE_H
#define REGISTERYFILE_H

#include <string>
#include <ctime>

#include "boost/shared_ptr.hpp"
#include "boost/filesystem.hpp"

#include "utilities/Reference.h"
#include "iomanagment/LocalPath.h"
#include "iomanagment/Dictionary2.h"

namespace SEM { namespace iomanagment {

/** **********************
* @enum READ
* Define file read options
*************************/
enum READ
{
    NO_READ=0, // File never read
    READ_ONCE, // File read once-if don't exist program throw error
    READ_IF_PRESENT, //File read once if exist
    READ_IF_MODIFIED, //File read again if was modified at runtime
    ALWAYS
};

/** **********************
* @enum WRITE
* Define file write options
*************************/
enum WRITE
{
    NO_WRITE=0, // Never write
    AUTO  // Write each time when Register call reading method on registred object
};

class Register;
class RegistryObject;

/** **********************
* \class RegistryFile
* Class for controling reading/writing
* with assocciated file. Due to use of
* LocalPath varible it's possible to
* automaticly change file location during
* program runtime.
* --------------------------------------
* Any object assocciated with this RegFile
* will be automaticly updated or writed
* into the file defined here by path:
* CaseLocation/LocalPath/objectName
*************************/
class RegistryFile
{
    REFERENCE_TYPE(RegistryFile)

    /// \var m_register
    /// object which controls moment when files shall be read/write
    Register & m_register;

    /// \var m_localPath
    /// local path provider--> provieds local directory path for this file
    LocalPath::ref m_localPath;

    /// \var m_objectName
    /// variable equal to fileName
    std::string m_objectName;

    /// \var m_read, m_write
    /// read/write file options
    READ m_read;
    WRITE m_write;


    /// \var lastModificationTime
    /// file property storing file last modification time
    std::time_t lastModificationTime;

    /// \var objects
    /// collection of registred objectes which
    /// are linked with this file
    std::vector<RegistryObject*> m_objects;

    /// \var currnetDict
    /// dictionary read from file for last time
    Dictionary currentDict;

    // Disalow copy/assigment constructors
    RegistryFile( const RegistryFile & other);
    RegistryFile& operator =(const RegistryFile & other);

public:
    /// Construct RegistryFile from setup and LocalPath pointer is build
    /// from string as defualt implementation
    /// --> due to this fact this constructor shall be called for files
    /// which will not change location during program runtime
    RegistryFile( Register& reg, const std::string & objectName,
                  const boost::filesystem::path& local,
                  const READ& readOpt, const WRITE& writeOpt );


    /// Construct RegistryFile from setup and LocalPath pointer
    /// -->local path due to polymorphism can provide non-constant local path
    RegistryFile( Register& reg, const std::string & objectName,
                  LocalPath::ref local,
                  const READ& readOpt, const WRITE& writeOpt );

    
    
    /// unregister from Registry
    virtual ~RegistryFile();

    /// return object name-> equal to fileName
    std::string fileName() const {return m_objectName;}

    /// subscribe object which will read/write file
    void registerObject(RegistryObject* regObj);

    /// unsubscribe object which will read/write file
    void unregisterObject(RegistryObject* regObj);

    RegistryObject* object(const std::string & objName) const;
    
    
//    /// return set register (right know this method is not neccessery)
//    inline Register & getRegister() const {return m_register;}

    /// full current file path
    boost::filesystem::path filePath() const;

    /// Call read on registrated objects
    void fireReading();

    /// Call write on registrated objects
    void fireWriting();

private:
    /// add this file into Register
    /// Read file(if allowed) at creation
    /// to pass this dictionary for new objects
    /// registred in this file
    void initFile();

    /// Check flags and file state for reading
    bool allowedRead();

    bool exist() const;
    
    /// Check flags and file state for writing
    bool allowdWrite();
};



}//iomanagment
}//SEM
#endif // REGISTERYFILE_H


