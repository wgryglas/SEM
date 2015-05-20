#ifndef DICTENTRY_H
#define DICTENTRY_H

#include "iomanagment/ReadWriteStream.h"

#include <sstream>

namespace SEM { namespace iomanagment {

/// \brief The DictEntry class
/// Class storing pairs of key-value read
/// from file. Key and value are held here
/// as std::string. To get desired type
/// user shall call template method getValue/operator>>
/// which will call apropriate reader for desired value type
/// This way of storing data allows very flexible file definitions.
/// None token must be earlier passed to parser, the only
/// thing is the expected type of value.
/// -----------------------------------------------------
/// Disadventage: storing data as string copy of the file part.
/// Possible in future data shall be moved to some sort
/// file-stream, and key-value then shall not store data directly,
/// butonly begin-end stream position. The only problem wchich
/// may arise is that, the file will be opened(may not, but
/// it would slow down processing).
class DictEntry
{
private:
    std::string m_key;
    std::string m_value;

public:
    /// \brief DictEntry  constructor
    /// \param key       -entry name
    /// \param value     -entry final value (final means that exacly this string will be printed)
    ///                   setValue and operator<<  allows to convert any type to string automaticly
    ///                  note: setValue and operator<< converts std::string=text to value="text",
    ///                   so to set not converted value user shall use this constructor or
    ///                   setAsValue method.
    DictEntry(const std::string& key, const std::string& value="");
    
    DictEntry(const DictEntry& d);
    
    DictEntry & operator =(const DictEntry& d);
    

    /// \brief ~DictEntry -virtual destructor(empty, no heap value alocated here)
    virtual ~DictEntry(){}

    /// \brief key getter for entry name
    /// \return string identifying this entry in dictionary
    inline std::string key() const {return m_key;}

    /// \brief value getter for value.
    /// \return value held by this entry. Value is in form as it will be printed/read from file
    inline std::string value() const {return m_value;}

    /// \brief operator = - sets value string directly
    ///                     --> method required due to fact, that std::string
    ///                         value is automaticly wraped by chars " " in
    ///                         setValue and operator << methods.
    /// \param valueToBeSet
    /// \return
    inline DictEntry& operator=(const std::string &valueToBeSet)
    {
        m_value = valueToBeSet;
        return *this;
    }

    /// \brief operator == compare some string with value held by this entry
    /// \param compareValue - std::string to compare
    /// \return  true if compareValue is equal to value held by this entry.
    bool operator == (const std::string& compareValue) const;

    /// \brief operator == compare some string with value held by this entry
    /// \param compareValue - std::string to compare
    /// \return  true if compareValue is not equal to value held by this entry.
    bool operator != (const std::string& compareValue) const;

    /// method which takes any type value and writes it to string value
    /// -->type resolved at compile time
    template<typename T> void setValue(const T& val);

    /// method which set provided string as value. 
    /// \note that setValue method converts any type(string including)
    /// to output string form (in case of string adds " to provided string )
    /// Below method do not convert value, provided value must be 
    /// in correct form for output.
    void setStrValue(const std::string & value);
    
    /// method reads string value to any type
    /// -->type resolved at compile time
    template<typename T> void getValue(T& val);
    
    /// method reads string value to any type
    /// -->type resolved at compile time.
    /// \note desired type require default constructor and assigment operator
    template<typename T> T getValue();

    /// \brief operator << - method to set value from any type (the same as setValue)
    /// \param - value to be written to this entry
    /// \return - this object reference
    template<typename T> DictEntry& operator <<(const T& val);

    /// \brief operator << - method to get value from any type (the same as getValue)
    /// \param - value to be read from this entry
    /// \return - this object reference
    template<typename T> DictEntry& operator >>(T& val);

    /// \brief operator << (cosnt vers.) method to get value from any type (the same as getValue)
    /// \param - value to be read from this entry
    /// \return - this object reference
    template<typename T> const DictEntry& operator >>(T& val) const;
    
    /// \brief writeEntry - method which format entry how should be placed in file
    /// \param writFile   - stream where entry shall be placed
    virtual void writeEntry(std::ostream & writFile) const;
};

template<typename T>
inline void DictEntry::setValue(const T& val)
{
    std::ostringstream ss;
    write<T>(val, ss);
    m_value=ss.str();
}

template<typename T>
inline void DictEntry::getValue(T& val)
{
    std::istringstream ss(m_value);
    read<T>(ss, val);
}

template<typename T>
inline T DictEntry::getValue()
{
    std::istringstream ss(m_value);
    T val;
    read<T>(ss, val);
    return val;
}


template<typename T>
DictEntry& DictEntry::operator <<(const T& val)
{
    std::ostringstream ss;
    write<T>(val,ss);
    m_value=ss.str();
    return *this;
}

template<typename T>
DictEntry& DictEntry::operator >>(T& val) 
{
    std::istringstream ss(m_value);
    read<T>(ss,val);
    return *this;
}

template<typename T>
const DictEntry& DictEntry::operator >>(T& val) const
{
    std::istringstream ss(m_value);
    read<T>(ss,val);
    return *this;
}



} //iomanagment
} //SEM

#endif // DICTENTRY_H
