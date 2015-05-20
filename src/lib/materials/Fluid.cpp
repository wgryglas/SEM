#include "Fluid.h"

namespace SEM { namespace materials {
    
    Fluid::Fluid(iomanagment::RegistryFile::ref regFile)
    : RegistryObject("Fluid")
    {
       setRegistryFile(regFile); 
    }
    
    Fluid::Fluid(iomanagment::RegistryFile* regFilePtr)
    : RegistryObject("Fluid")
    {
        iomanagment::RegistryFile::ref regFile(regFilePtr);
        setRegistryFile(regFile);
    }
    
    Fluid::Fluid(const Fluid& other)
    : RegistryObject("Fluid")
    {
        setRegistryFile(other.m_file);
    }

    Fluid::~Fluid()
    {
    }

    Fluid & Fluid::operator=(const Fluid& other)
    {
    }

}//materials
}//SEM