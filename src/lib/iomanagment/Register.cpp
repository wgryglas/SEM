#include "iomanagment/Register.h"

#include <algorithm>

namespace SEM {namespace iomanagment {

void Register::registerFile(RegistryFile *obj)
{
    files.push_back(obj);
}

void Register::unregisterFile(RegistryFile *obj)
{
    std::vector<RegistryFile*>::iterator itr;
    itr=std::find(files.begin(),files.end(),obj);

    if(itr!=files.end())
        files.erase(itr);
}

void Register::fireReading()
{
    for(RegistryFile* obj: files)
        obj->fireReading();
}

void Register::fireWriting()
{
    for(RegistryFile* obj : files)
        obj->fireWriting();
}


} //iomanagment
} //SEM
