
#include "iomanagment/Tokenizer.h"
#include <iostream>
#include <fstream>
//#include "boost/algorithm/string/split.hpp"
#include "boost/tokenizer.hpp"
#include "boost/regex.hpp"

namespace SEM
{
namespace IOStream
{
  using namespace std;

  Tokenizer::Tokenizer():comment_del_("/"), split_del_("{}()<>,;[]")
  {
  }
//  Tokenizer::Tokenizer():comment_del_("(//)"), split_del_("\\{|\\}|\\(|\\)|\\<|\\>|\\,|\\;|\\[|\\]")
//  {
//  }


  Tokenizer::Tokenizer(const string comment_del, const string split_del): comment_del_(comment_del), split_del_(split_del)
  {
  }

  Tokenizer::~Tokenizer()
  {
  }

  void Tokenizer::parse(const char* fileName, vector<string>& result)
  {
    ifstream file(fileName, ifstream::in);
    if(!file.is_open())
    {
      cout<<"Error file parsing. File <"<<fileName<<"> couldn't be opened!"<<endl;
      return;
    }

    boost::char_separator<char> comment_sep("", comment_del_.data());
    boost::char_separator<char> split_sep("\t ",split_del_.data());//\t\n
    typedef boost::tokenizer<boost::char_separator<char> > token;

    string line;
    int wordCount = 0;
    int lineCount = 0;

    while(getline(file, line))
    {
      lineCount++;
      //remove line comment:
      token tok_com(line, comment_sep);
      for(token::iterator tok_com_itr=tok_com.begin(); tok_com_itr!=tok_com.end(); ++tok_com_itr)
      {
          if(tok_com_itr!=tok_com.begin())
          {
            line = *(tok_com.begin());
            break;
          }
      }
      //split non-commented line part:
      token tok_split(line, split_sep);
      for(token::iterator it=tok_split.begin(); it!=tok_split.end(); it++)
      {
         result.push_back(*it);
      }

    }

//    boost::regex re_comment("//");
//    boost::regex re_split("\\{|\\}|\\(|\\)|\\<|\\>|\\,|\\;|\\[|\\]");
//    //boost::regex re_comm(comment_del_.data());
//    //boost::regbase re_split(split_del_.data());
//    while(getline(file, line))
//    {
//      //remove comment
//      {
//        boost::sregex_token_iterator i_c(line.begin(), line.end(), re_comment, -1);
//        boost::sregex_token_iterator j_c(i_c);
//        boost::sregex_token_iterator k;
//        if(++j_c!=k)
//        {
//          line = *i_c;
//        }
//
//        boost::sregex_token_iterator i_s(line.begin(), line.end(), re_split, -1);
//        while(i_s!=k)
//        {
//          result.push_back(*i_s);
//          i_s++;
//        }
//
//      }
//    }

//    for(vector<string>::iterator s_itr=result.begin(); s_itr!=result.end(); s_itr++)
//    {
//      cout<<*s_itr<<endl;
//    }






  }


}//IOStream
}//SEM
