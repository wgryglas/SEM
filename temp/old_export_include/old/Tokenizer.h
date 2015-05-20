#ifndef Tokenizer_H
#define Tokenizer_H

#include <vector>
#include <string>


namespace SEM
{
namespace IOStream
{
using namespace std;

class Tokenizer
{
  const string comment_del_;
  const string split_del_;

  Tokenizer(const Tokenizer& tok){};
  void operator=(const Tokenizer& tok){};

public:
  Tokenizer();
  Tokenizer(const string comment_del, const string split_del );

  virtual ~Tokenizer();

  void parse(const char* fileName, vector<string>& result);

};

} //IOStream
} //SEM
#endif //Tokenizer_H
