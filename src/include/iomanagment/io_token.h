#ifndef IO_TOKEN_H
#define IO_TOKEN_H

#include<string>
#include "boost/smart_ptr/shared_ptr.hpp"

namespace SEM
{
namespace iomanagment
{
  using std::string;

  class IOToken
  {
    IOToken(IOToken& file)
    {
    }

    void operator=(IOToken& file)
    {
    }

    public:

   // IOToken(const string file_name);
    virtual void manage(string s_token, boost::shared_ptr<Dictionary> current_dic);


  };

}//IOMANAGMENT
}//SEM


#endif //IO_TOKEN_H

