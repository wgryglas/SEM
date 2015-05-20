#include "DictEntry.h"

namespace SEM { namespace iomanagment {

DictEntry::DictEntry(const std::string & key, const std::string & value):
    m_key(key), m_value(value)
{
}

DictEntry::DictEntry(const DictEntry& d) 
: m_key(d.m_key), m_value(d.m_value)
{
    
}

DictEntry& DictEntry::operator=(const DictEntry& d) 
{
    m_key = d.m_key;
    m_value =d.m_value;
    return *this;
}

bool DictEntry::operator ==(const std::string &compareValue) const
{
    return m_value==compareValue;
}

bool DictEntry::operator !=(const std::string &compareValue) const
{
    return m_value!=compareValue;
}

void DictEntry::writeEntry(std::ostream & of) const
{
    of<<m_key<<"    "<<m_value<<";";
}
void DictEntry::setStrValue(const std::string& value) 
    {
        m_value = value;
    }




}//iomanagment
}//SEM
