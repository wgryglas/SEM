
#ifndef DATA_ENTRY_H
#define DATA_ENTRY_H

#include <string>
#include <vector>

namespace SEM
{
namespace iomanagment
{
  using std::vector;
  using std::string;

  class DataEntry
  {
  public:
    enum DataType
    {
      kScalar,
      kVector,
      kTensor,
      kTable,
      kList
    };

    DataEntry()
    {
    }

    virtual bool is_container() const;
    virtual DataType get_data_type() const;
    virtual void add_entry(const DataEntry& entry);
    virtual void add_data(const string data);
    virtual const vector<DataEntry>& get_data_entries() const;
    virtual const vector<string>& get_data() const;

    virtual ~DataEntry()
    {
    }

  };
} //IOMANAGMENT
} //SEM

#endif //DATA_ENTRY_H

