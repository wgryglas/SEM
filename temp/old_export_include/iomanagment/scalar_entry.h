
#ifndef SCALAR_ENTRY_H
#define SCALAR_ENTRY_H

#include "data_entry.h"
#include "leaf_entry.h"
#include<iostream>

namespace SEM
{
namespace IOMANAGMENT
{
  using std::vector;
  using std::string;

  class ScalarEntry: public LeafEntry
  {
    const vector<string> value_;

    ScalarEntry(ScalarEntry& sE): value_(sE.get_data())
    {
    }

    void operator=(ScalarEntry& sE)
    {
    }

  public:

    ScalarEntry(const string value): value_(1, value)
    {
    }

    virtual ~ScalarEntry();

    virtual void add_data(const string data)
    {
      std::cout<<"WARNING:: can't add data to scalar entry"<<std::cout;
    }

    virtual DataType get_data_type()
    {
      return DataEntry::kScalar;
    }

    virtual const vector<string>& get_data() const
    {
      return value_;
    }

  };

} //IOMANAGMENT
} //SEM




#endif//SCALAR_ENTRY_H
