#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <string>
#include <vector>
#include <map>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "word.h"

namespace SEM
{
namespace iomanagment
{
  class DictionaryOld
  {
    typedef boost::shared_ptr<DictionaryOld> dic_ptr;
    typedef boost::shared_ptr<Word> word_ptr;
    typedef std::pair<std::string, dic_ptr> dic_pair;
    typedef std::pair<std::string, word_ptr> word_pair;

    std::map<std::string, word_ptr > words_;
    std::map<std::string, dic_ptr > dictionaries_;

    std::string name_;

  public:

    DictionaryOld();
    DictionaryOld(std::string name);
    virtual ~DictionaryOld();

    void set_name(const std::string name);

    std::string get_name() const;

    void add_word(Word* const w_ptr);

    void add_dictionary(DictionaryOld* const d_ptr);

    std::vector< dic_ptr > get_dictionaries() const;

    std::vector< word_ptr > get_words() const;

    bool contain_dictionary(const std::string name) const;

    dic_ptr get_dictionary(const std::string name) const;

    bool contain_word(const std::string name) const;

    word_ptr get_word(const std::string name) const;

  };
} //IOMANAGMENT
} //SEM

#endif //DICTIONARY_H

