
#ifndef TABLE_ENTRY_H
#define TABLE_ENTRY_H

#include "data_entry.h
#include "leaf_entry.h"

namespace SEM
{
namespace iomanagment
{
  class TableEntry: public LeafEntry
  {
    const std::vector<std::string> value_;

    TableEntry(TableEntry& vE): value_(vE.get_data())
    {
    }

    void operator=(TableEntry& vE)
    {
    }

  public:

    TableEntry(const vector<string> data): value_(data)
    {
    }

    virtual ~TableEntry();

    virtual void add_data(const string data)
    {
      Svalue_.push_back(data);
    }

    virtual DataType get_data_type()
    {
      DataEntry::kTable;
    }

    virtual const std::vector<std::string>& get_data() const
    {
      return value_;
    }

  };

} //IOMANAGMENT
} //SEM

#endif//TABLE_ENTRY_H


