
#include "dictionary.h"
#include <iostream>

namespace SEM
{
namespace iomanagment
{
  DictionaryOld::DictionaryOld()
  {
  }

  DictionaryOld::DictionaryOld(const std::string name): name_(name)
  {
  }

  DictionaryOld::~DictionaryOld()
  {
  }

  void DictionaryOld::set_name(const std::string name)
  {
    name_ = name;
  }

  std::string DictionaryOld::get_name() const
  {
    return name_;
  }

  void DictionaryOld::add_word(Word* const w_ptr)
  {
    words_[w_ptr->get_name()]= word_ptr(w_ptr);
  }

  void DictionaryOld::add_dictionary(DictionaryOld* const d_ptr)
  {
    dictionaries_[d_ptr->get_name()]= dic_ptr(d_ptr);
  }

  std::vector< DictionaryOld::dic_ptr > DictionaryOld::get_dictionaries() const
  {
    std::vector<dic_ptr> dic_vec( dictionaries_.size() );
    int i=0;
    for(std::map<const std::string, dic_ptr >::const_iterator it=dictionaries_.begin(); it!=dictionaries_.end(); it++)
    {
      dic_vec[i] = it->second;
      i++;
    }

    return dic_vec;
  }

  std::vector<DictionaryOld::word_ptr > DictionaryOld::get_words() const
  {
    std::vector<word_ptr> word_vec( words_.size() );

    int i=0;
    for(std::map< const std::string, word_ptr >::const_iterator it=words_.begin(); it!=words_.end(); it++)
    {
      word_vec[i] = it->second;
      i++;
    }
    return word_vec;
  }

  bool DictionaryOld::contain_dictionary(const std::string name) const
  {
      return dictionaries_.find(name)!=dictionaries_.end();
  }

  DictionaryOld::dic_ptr DictionaryOld::get_dictionary(const std::string name) const
  {
    return dictionaries_.find(name)->second;
  }

  bool DictionaryOld::contain_word(const std::string name) const
  {
    return words_.find(name)!=words_.end();
  }

  DictionaryOld::word_ptr DictionaryOld::get_word(const std::string name) const
  {
    return words_.find(name)->second;
  }

} //IOMANAGMENT
} //SEM
