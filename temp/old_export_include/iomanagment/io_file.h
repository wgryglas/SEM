
#ifndef IO_FILE_H
#define IO_FILE_H

#include<string>

#include "file_processor.h"
#include "dictionary.h"
#include "word.h"
#include "data_entry.h"

namespace SEM
{
namespace IOMANAGMENT
{
  using std::string;

  class IOFile : public Dictionary
  {
    IOFile(IOFile& file)
    {
    };

    void operator=(IOFile& file)
    {
    };

    public:

    IOFile(const string file_name);



  };

}//IOMANAGMENT
}//SEM







#endif //IO_FILE_H
