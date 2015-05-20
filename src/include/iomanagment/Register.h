#ifndef REGISTER_H
#define REGISTER_H

#include <vector>
#include <string>

#include "RegistryObject.h"
#include "RegistryFile.h"

namespace SEM { namespace iomanagment {

class Register
{
    std::vector<RegistryFile*> files;
public:

    void registerFile(RegistryFile* entry);
    void unregisterFile(RegistryFile* entry);
    
    virtual void fireReading();
    virtual void fireWriting();
    
    
    template<typename T>
    T* object(const std::string &name) const
    {
        for(RegistryFile* file: files)
        {
            if(RegistryObject* obj=file->object(name))
            {
                if(T* typeObj = dynamic_cast<T*>(obj))
                    return typeObj;
                else
                {
                    using namespace SEM::iomanagment;
                    ErrorInFunction << "Object "<<name<<" found in register, but type of this object is diffrent then desired"<<endProgram;
                }
            }
        }
        
        ErrorInFunction << " object "<<name<<" was not found in register and it's required by program"<<endProgram;
        
        return nullptr;
    }
    

};








} //iomanagment
} //SEM


#endif // REGISTER_H
