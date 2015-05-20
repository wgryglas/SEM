#ifndef REGISTRYOBJECT_H
#define REGISTRYOBJECT_H

#include "iomanagment/RegistryFile.h"
#include "iomanagment/Dictionary2.h"

namespace SEM { namespace iomanagment {

///////////////////////////////////////
/// @class RegistryObject
/// ----------------------------------
/// Base class for any object that
/// shall be linked with file.
/// ----------------------------------
/// Note:: actual registration is done
/// in setRegistryFile method, so in
/// this manner user can easly chose
/// if want to allow automatic object
/// filling with file data or when
/// setRegistryFile method is not called
/// then such object will not be linked
/// with file, so can be used in normal
/// way.
////////////////////////////////////////
class RegistryObject
{
protected: 
    /// @variable m_file
    /// Linked RegistryFile reference
    /// note: used ref. because one instance
    /// of file may be linked with many objects
    /// --> for purpose when you want to control
    /// many objects by single file.
    RegistryFile::ref m_file;
    
    std::string m_objectName;
    
public:
    RegistryObject(const std::string &name);

    virtual ~RegistryObject();
    
    /// @arg refFile- file where obj. shall be link to
    /// method for actual object linking with file
    void setRegistryFile(RegistryFile::ref refFile);
    
    std::string objectName() const { return m_objectName;}
    
    RegistryFile::ref registryFile()
    {
        return m_file;
    }

    /// @arg dict - object read from file for furhter processing.
    /// method to be called when object shall be updated
    /// (always called at begining, when object is being registrated in regFile
    virtual void read(const Dictionary& dict)=0;

    /// @arg dict - dictionary where data from this object shall be palced in
    /// In this level dict is whole contest of file.
    /// Methid called when registry will decide.
    /// This method shall be implemented by derived object for
    /// handling writing object contest.
    virtual void write(Dictionary& dict)=0;
    
protected:
    RegistryFile::ref file() const {return m_file;}
};




}//iomanagment
}//SEM


#endif // REGISTRYOBJECT_H
