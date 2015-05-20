
#include "word.h"

namespace SEM
{
namespace iomanagment
{
  using std::vector;
  using std::string;

  Word::Word(const string name, int line_id): name_(name), line_id_(line_id_)
  {
  }

  Word::~Word()
  {
  }

  string Word::get_name() const
  {
      return name_;
  }

  int Word::get_line() const
  {
      return line_id_;
  }

  const vector<boost::shared_ptr<DataEntry> > Word::get_data_entries() const
  {
      return data_entries_;
  }

  void Word::add_data_entry(DataEntry* const data_ptr)
  {
      data_entries_.push_back( boost::shared_ptr<DataEntry>(data_ptr) );
  }


} //IOMANAGMENT
} //SEM
