#include "Dictionary2.h"

#include "iomanagment/InfoStream.h"
#include "iomanagment/ReadWriteStream.h"
#include "iomanagment/Indenter.h"
#include "iomanagment/FileProcessing.h"

namespace SEM { namespace iomanagment {

Dictionary::Dictionary(const std::string& dicName)
 : m_dicName(dicName), m_parent(NULL)
{
}

Dictionary::Dictionary(const Dictionary& other, Dictionary* parent) 
: m_dicName(other.m_dicName), m_parent(parent) 
{
    for(auto e : other.m_entries)
    {
        m_entries.insert(entryPair(e.first,new DictEntry(*e.second)));
    }
    
    for(auto e : other.m_subDictionaries)
    {
        m_subDictionaries.insert(dictionaryPair(e.first,new Dictionary(*e.second,this)));
    }
}

Dictionary& Dictionary::operator=(const Dictionary& other) 
{
    clear();
    
    for(auto e : other.m_entries)
    {
        m_entries.insert(entryPair(e.first,new DictEntry(*e.second)));
    }
    
    for(auto e : other.m_subDictionaries)
    {
        m_subDictionaries.insert(dictionaryPair(e.first,new Dictionary(*e.second,this)));
    }
}

Dictionary::~Dictionary()
{
   clear();
}

Dictionary *Dictionary::getParent()
{
    return m_parent;
}


void Dictionary::add(Dictionary* subDict)
{
	subDict->setParent(this);
    m_subDictionaries.insert( std::pair<std::string,Dictionary*>(subDict->name(), subDict) );
}

void Dictionary::add(DictEntry* entry)
{
    if(hasEntry(entry->key()))
    {
        delete m_entries[entry->key()];
        m_entries.erase(entry->key());
    }
    
    m_entries.insert(std::pair<std::string,DictEntry*>(entry->key(),entry));
}

bool Dictionary::hasEntry(const std::string &key) const
{
    entryMap::const_iterator itr= m_entries.find(key);
    if(itr!=m_entries.end())
    {
        return true;
    }

    return false;
}

DictEntry &Dictionary::entry(const std::string &key) const
{
    entryMap::const_iterator itr= m_entries.find(key);
    if(itr!=m_entries.end())
    {
        return *(itr->second);
    }
    else
    {
        ErrorInFunction<<"Entry with key "<<key<<" was not found in dictionary "<<m_dicName<<endProgram;
    }
}

bool Dictionary::hasSubDictionary(const std::string &key) const
{
    dictionarMap::const_iterator itr= m_subDictionaries.find(key);
    if(itr == m_subDictionaries.end())
    {
        return false;
    }
    
    return true;
}

Dictionary &Dictionary::subDictionary(const std::string &dictName) const
{
    dictionarMap::const_iterator itr= m_subDictionaries.find(dictName);

    if(itr!=m_subDictionaries.end())
    {
        return *(itr->second);
    }
    else
    {
        ErrorInFunction<<"Dictionary with name "<<dictName<<" was not found in dictionary "<<m_dicName<<endProgram;
    }
}

void Dictionary::setParent(Dictionary* parent)
{
    m_parent = parent;
}

bool Dictionary::hasParent() const
{
    return m_parent!= NULL;
}

void Dictionary::writeDicsAndEntries(std::ostream & of) const
{
    using namespace std;

    if(m_entries.size()>0)
    {
        //DictEntry* last=entries.back();
        for(const entryMap::value_type& e: m_entries)
        {
            e.second->writeEntry(of);
//            if(e!=last)
            of<<endl;
        }
    }

    if(m_subDictionaries.size()>0)
    {
        of<<endl;
        //Dictionary* lastDict = subDictionaries.back();
        for(const dictionarMap::value_type &d :m_subDictionaries)
        {
            d.second->write(of);
//            if(d!=lastDict)
                of<<endl;
        }
    }
}

void Dictionary::wrtiAsFile() const
{
    using namespace std;

    if(hasParent())
    {
        ErrorInFunction<<"can't write dictionary where it is not top level dictionary"<<endProgram;
    }

    ofstream file(m_dicName.c_str(),ios::trunc);
    if(!file)
    {
        ErrorInFunction<<"coudn't write file"<<m_dicName.c_str()<<endProgram;
    }

    stringstream textStream;
    writeFileDef(textStream);
    writeDicsAndEntries(textStream);

    Indenter().indetStream<ofstream>(textStream, file);
    file.close();
}

void Dictionary::write(std::ostream & out) const
{
    using std::endl;

    if(hasParent())
    {
        out<<m_dicName<<endl
           <<"{"<<endl;
    }

    writeDicsAndEntries(out);

    if(hasParent())
    {
        out<<endl
           <<'}'<<endl;
    }
}

void Dictionary::clear()
{
    for(dictionarMap::value_type p : m_subDictionaries)
        delete p.second;

    for(entryMap::value_type p : m_entries)
        delete p.second;

    m_subDictionaries.clear();
    m_entries.clear();
}

void Dictionary::writeFileDef(std::ostream & of) const
{
    //file header shall be placed here
}

std::ostream & operator << (std::ostream & stream, const Dictionary & dict )
{
    std::stringstream ss;
    dict.write(ss);
    Indenter().indetStream<std::ostream>(ss,stream);

    return stream;
}

std::istream & operator >>(std::istream & stream, Dictionary & dict)
{
    using namespace std;

    stream >> noskipws;

    StreamProcessingIterator itr(stream);
    StreamProcessingIterator end(stream,true);

    std::vector<char> tokens = {'{','}',';'};
    SEMDictTokenizer tokenizedFile(itr,end,tokens);

    Dictionary* current = &dict;
    for(SEMDictTokenizer::iterator itr = tokenizedFile.begin(); itr!=tokenizedFile.end(); itr++ )
    {
        switch( itr->first )
        {
        case '{':
        {
            Dictionary* subDict=new Dictionary(itr->second);
            current->add(subDict);
            current = subDict;
        break;
        }
        case '}':
        {
            current = current->getParent();
        break;
        }
        case ';':
        {
            int wsPos  = itr->second.find_first_of(' ');
            string key = itr->second.substr(0,wsPos);
            string value = itr->second.substr(wsPos+1,itr->second.length()-wsPos-1);

            DictEntry* entry = new DictEntry(key,value);
            current->add(entry);
        break;
        }
        }
    }

    return stream;
}

Dictionary& operator >>(const boost::filesystem::path& file, SEM::iomanagment::Dictionary& dict)
{
    using namespace std;

    if( ! boost::filesystem::exists(file))
    {
        ErrorInFunction<<"Specified file "<<file<<" don't exist"<<endProgram;
    }
    else if( boost::filesystem::is_directory(file) )
    {
        ErrorInFunction<<"Specified file "<<file<<" is directory, where program expects here file"<<endProgram;
    }

    boost::filesystem::ifstream file_stream(file);

    dict.m_dicName = file.stem().string();

    file_stream >> dict;

    file_stream.close();

    return dict;
}

Dictionary& operator <<(const boost::filesystem::path& file, SEM::iomanagment::Dictionary& dict)
{
    using namespace std;

    if( boost::filesystem::is_directory(file) )
    {
        ErrorInFunction<<"Specified file is directory, where program expects here file"<<endProgram;
    }

    boost::filesystem::path parentDir = file.parent_path();
    if( ! boost::filesystem::exists( parentDir ) )
    {
         boost::filesystem::create_directory( parentDir );
    }

    boost::filesystem::ofstream file_stream(file);

    file_stream << dict;

    file_stream.close();

    return dict;
}


}//IOMANAGMENT
}//SEM


