
#ifndef FLUID_H
#define FLUID_H

#include "iomanagment/RegistryObject.h"
#include "Property.h"

namespace SEM { namespace materials {

class Fluid : public iomanagment::RegistryObject
{
    
    
    
protected:
    using iomanagment::RegistryObject::m_file;
public:
    Fluid(iomanagment::RegistryFile::ref regFile);
    
    Fluid(iomanagment::RegistryFile* regFilePtr);
    
    /// \brief Copy constructor. Automaticly make 
    /// registration in file 
    Fluid(const Fluid& other);
    
    
    ~Fluid();
    Fluid& operator=(const Fluid& other);

};



}//materials
}//SEM

#endif // FLUID_H
