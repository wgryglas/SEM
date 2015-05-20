#ifndef LEAF_ENTRY_H
#define LEAF_ENTRY_H

#include <vector>
#include <string>
#include <iostream>

#include "data_entry.h"

namespace SEM
{
namespace iomanagment
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

    const std::vector<DataEntry>& get_data_entries() const
    {
      return std::vector<DataEntry>(0);
    }

    virtual const std::vector<string>& get_data() const;

  };

} //IOMANAGMENT
} //SEM

#endif //LEAF_ENTRY_H
