
#ifndef WORD_H
#define WORD_H

#include <string>
#include <vector>

#include "boost/smart_ptr/shared_ptr.hpp"

#include "iomanagment/data_entry.h"
#include "iomanagment/scalar_entry.h"
#include "iomanagment/vector_entry.h"
#include "iomanagment/tensor_entry.h"
#include "iomanagment/list_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  using std::vector;
  using std::string;

  class Word
  {
    const string name_;
    const int line_id_;
    vector< boost::shared_ptr<DataEntry> > data_entries_;

    Word(Word& w): name_(w.get_name()), line_id_(w.get_line())
    {
    };

    void operator=(Word& w)
    {
    };

  public:

    Word(const string name, int line_id);

    string get_name() const;
    int get_line() const;
    void add_data_entry(DataEntry* const data_ptr);
    const vector<boost::shared_ptr<DataEntry> > get_data_entries() const;

    virtual ~Word();

  };
} //IOMANAGMENT
} //SEM

#endif //WORD_H
