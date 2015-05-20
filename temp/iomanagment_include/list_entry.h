
#ifndef LIST_ENTRY_H
#define LIST_ENTRY_H

#include "iomanagment/container_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  using std::string;
  using std::vector;

  class ListEntry: public ContainerEntry
  {
    vector<DataEntry> data_entries_;

    ListEntry(ListEntry& lE)
    {
    };

    void operator= (ListEntry& lE)
    {
    };

    public:

    ListEntry();

    virtual ~ListEntry();


    DataType get_data_type() const
    {
      return DataEntry::kList;
    };

    virtual void add_entry(const DataEntry& entry)
    {
        data_entries_.push_back(entry);
    };

    virtual const vector<DataEntry>& get_data_entries() const
    {
      return data_entries_;
    };

  };
}//IOMANAGMENT
}//SEM



#endif //LIST_ENTRY_H
