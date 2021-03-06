
#ifndef VECTOR_ENTRY_H
#define VECTOR_ENTRY_H

#include "iomanagment/data_entry.h"
#include "iomanagment/leaf_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  class VectorEntry: public LeafEntry
  {
    vector<string> value_;

    VectorEntry(VectorEntry& vE): value_(vE.get_data())
    {
    };

    void operator=(VectorEntry& vE)
    {
    };

  public:

    VectorEntry()
    {
    };

    VectorEntry(const string x, const string y): value_(2,x)
    {
      value_[1] = y;
    };

    virtual ~VectorEntry();

    virtual void add_data(const string data)
    {
      if(value_.size()<2)
      {
        value_.push_back(data);
      }
      else
      {
        std::cout<<"WARNING:: can't add more data to vector entry"<<std::cout;
      }
    };


    virtual DataType get_data_type()
    {
      return DataEntry::kVector;
    };

    virtual const vector<string>& get_data() const
    {
      return value_;
    };

  };

}//IOMANAGMENT
}//sem

#endif//VECTOR_ENTRY_H

