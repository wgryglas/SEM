#ifndef LEAF_ENTRY_H
#define LEAF_ENTRY_H

#include <vector>
#include <string>
#include <iostream>

#include "iomanagment/data_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  using std::vector;
  using std::string;

  class LeafEntry : public DataEntry
  {

  public:

    virtual ~LeafEntry();

    bool is_container() const
    {
      return false;
    }

    virtual DataType get_data_type() const;

    void add_entry(const DataEntry& entry)
    {
      std::cout<<"Can't add entry to leaf entry, only data values are permited"<<std::endl;
    }

    virtual void add_data(const string data);

    const vector<DataEntry>& get_data_entries() const
    {
      return vector<DataEntry>(0);
    }

    virtual const vector<string>& get_data() const;

  };

} //IOMANAGMENT
} //SEM

#endif //LEAF_ENTRY_H
