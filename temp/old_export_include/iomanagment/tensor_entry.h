

#ifndef TENSOR_ENTRY_H
#define TENSOR_ENTRY_H

#include "data_entry.h"
#include "leaf_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  class TensorEntry: public LeafEntry
  {
    vector<string> value_;

    TensorEntry(TensorEntry& vE): value_(vE.get_data())
    {
    };

    void operator=(TensorEntry& vE)
    {
    };

  public:

    TensorEntry(const string xx, const string xy, const string yx, const string yy): value_(4,xx)
    {
      value_[1] = xy;
      value_[2] = yx;
      value_[3] = yy;
    };

    virtual ~TensorEntry();

    virtual void add_data(const string data)
    {
      if(value_.size()<4)
      {
        value_.push_back(data);
      }
      else
      {
        std::cout<<"WARNING:: can't add more data to tensor entry"<<std::cout;
      }
    };

    virtual DataType get_data_type()
    {
        return DataEntry::kTensor;
    };

    virtual const vector<string>& get_data() const
    {
      return value_;
    };

  };

}
}




#endif//TENSOR_ENTRY_H

