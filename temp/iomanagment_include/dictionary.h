
#ifndef DICTIONARY_H
#define DICTIONARY_H

#include <string>
#include <vector>
#include <map>

#include "boost/smart_ptr/shared_ptr.hpp"
#include "iomanagment/word.h"

namespace SEM
{
namespace IOMANAGMENT
{
  using std::vector;
  using std::string;
  using std::map;

  class Dictionary
  {
    typedef boost::shared_ptr<Dictionary> dic_ptr;
    typedef boost::shared_ptr<Word> word_ptr;
    typedef std::pair<string, dic_ptr> dic_pair;
    typedef std::pair<string, word_ptr> word_pair;

    map<string, word_ptr > words_;
    map<string, dic_ptr > dictionaries_;

    string name_;

  public:

    Dictionary();
    Dictionary(string name);
    virtual ~Dictionary();

    void set_name(const string name);

    string get_name() const;

    void add_word(Word* const w_ptr);

    void add_dictionary(Dictionary* const d_ptr);

    vector< dic_ptr > get_dictionaries() const;

    vector< word_ptr > get_words() const;

    bool contain_dictionary(const string name) const;

    dic_ptr get_dictionary(const string name) const;

    bool contain_word(const string name) const;

    word_ptr get_word(const string name) const;



  };
} //IOMANAGMENT
} //SEM

#endif //DICTIONARY_H

