#ifndef DICTIONARY_H_
#define DICTIONARY_H_

#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <map>
#include <sstream>

#include "boost/filesystem/path.hpp"

#include "iomanagment/DictEntry.h"
#include "utilities/Utilities.h"

namespace SEM{ namespace iomanagment{

/// \brief The Dictionary class
///  base storage for read data from file.
///  -------------------------------------------------------------
///  File is being splited into separate
///  dictionaries by definition:  dicName { ... }.
///  Inside each dictionary (between {})
///  parser reads "key value;" pairs and inner dictionaries.
///  Key-value pair is stored here as DictEntry type.
///  -------------------------------------------------------------
///  Desired data read from file is stored
///  in DictEntry "value" member variable.
///  DictEntry have special methods, which
///  automaticly converts any string "value"
///  to desired type.
///  -------------------------------------------------------------
///  Note 1: While reading/writing file user have
///  to keep in mind that dictionary name is equal
///  to file-name. If you would set some name to new dictionary,
///  and then you will read it from file, then your dictionary
///  will change name to the file name.
///  -------------------------------------------------------------
///  Note 2: File is assumed to be the first level of dictionary.
///  So if you specify file, then you shouldn't write your
///  dictionary name and first { } tokens - eg. file :
///  ---- File test.sem ---------
///  |                          
///  | key somthing;            
///  | innerDict                
///  | {                        
///  |   key2 somthing2;        
///  | }                        
///  ------ EOF test.sem---------
///   from library user view above will look like:
/// ----------------------------
/// | test                      
/// |  {                        
/// |     key somthing;         
/// |                           
/// |    innerDict              
/// |    {                      
/// |        key2 somthing2;    
/// |    }                      
/// | }                         
/// -----------------------------
///  ---> name of top level dictionary will be equal to "test"
class Dictionary
{

public:
    typedef std::pair<std::string,Dictionary*> dictionaryPair;
    typedef std::pair<std::string,DictEntry*> entryPair;
    typedef std::map<std::string,Dictionary*> dictionarMap;
    typedef std::map<std::string,DictEntry*> entryMap;

private:
    std::string m_dicName;
    Dictionary* m_parent;
    dictionarMap m_subDictionaries;
    entryMap m_entries;

    friend std::ostream& SEM::iomanagment::operator <<(std::ostream & o, const Dictionary & dict);

    friend Dictionary& operator >>(const boost::filesystem::path& file, SEM::iomanagment::Dictionary& dict);

public:

    /// \brief Dictionary
    /// \param dicName -this di
    Dictionary(const std::string &dicName="");
    
    Dictionary(const Dictionary & other, Dictionary * parent=nullptr);
    
    Dictionary & operator=(const Dictionary & other);
    

    virtual ~Dictionary();

    inline std::string name() const { return m_dicName;}

    Dictionary* getParent();

    void setParent(Dictionary* parent);

	bool hasParent() const;

    void add(Dictionary* subDict);

    /// \brief add new entry to dictionary
    /// if entyr already exist then it's 
    /// removed and destroyied. 
    /// \param entry - pointer to entry 
    /// to be inserted.
	void add(DictEntry* entry);

    bool hasEntry(const std::string &key) const;

    DictEntry& entry(const std::string & key) const;

    bool hasSubDictionary(const std::string & key) const;

    Dictionary& subDictionary(const std::string & dictName) const;

    inline const dictionarMap & subDictionaries() const{ return m_subDictionaries; }

    inline const entryMap & entries() const { return m_entries; }

    void wrtiAsFile() const;
    void write(std::ostream & out ) const;

    void clear();

private:
    void writeDicsAndEntries(std::ostream & of) const;

    void writeFileDef(std::ostream & of) const;

};

std::ostream & operator << (std::ostream & o, const Dictionary & dict );

/////////////////////////////////////////////
/// operator:: istream >> Dictionary
/// -----------------------------------------
/// function for creating dictionary object
/// from input stream object
/////////////////////////////////////////////
std::istream & operator >>(std::istream & stream, Dictionary & dict);

/////////////////////////////////////////////
/// operator:: path >> Dictionary
/// -----------------------------------------
/// function for creating dictionary object
/// from specified file
/// -----------------------------------------
/// Declared as friend of Dictionary
/////////////////////////////////////////////
Dictionary& operator >>(const boost::filesystem::path& file, Dictionary& dict);

/////////////////////////////////////////////
/// operator:: path >> Dictionary
/// -----------------------------------------
/// function for writing dictionary object
/// to specified file
/////////////////////////////////////////////
Dictionary& operator <<(const boost::filesystem::path& file, Dictionary& dict);


}//IOMANAGMENT
}//SEM


#endif /* DICTIONARY_H_ */
