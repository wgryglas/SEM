#ifndef CONTAINER_ENTRY_H
#define CONTAINER_ENTRY_H

#include <vector>
#include <string>

#include "data_entry.h"

namespace SEM
{
namespace iomanagment
{
  using std::vector;
  using std::string;

  class ContainerEntry : public DataEntry
  {

  public:


    virtual ~ContainerEntry();

    bool is_container() const
    {
      return true;
    }

    virtual DataType get_data_type() const;

    virtual void add_entry(const DataEntry& entry);

    void add_data(const string data)
    {
      std::cout<<"WARNING:: can't add data to container entry - only DataEntry is stored"<<std::cout;
    }

    const std::vector<std::string>& get_data() const
    {
      return std::vector<std::string>(0);
    }

  };

} //IOMANAGMENT
} //SEM

#endif //CONTAINER_ENTRY_H
