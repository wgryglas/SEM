
#ifndef _FileProcessor_H
#define _FileProcessor_H

#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "boost/tokenizer.hpp"
#include "boost/regex.hpp"

namespace SEM
{
namespace iomanagment
{
  using std::string;
  using std::vector;
  using std::cout;
  using std::cin;

	class FileProcessor
	{
	  vector<string> tokens_;

    static const boost::char_separator<char> delimiters_;
    static const boost::regex line_comment_delimiter_;

    int itr_;

  public:

    FileProcessor( const string filePath );

		string get_token();

		string get_token(int offset);

		void go_foreward();

		void go_foreward( int step );

		void go_to_begin();

		bool go_to( int place );

		bool go_to( const string token_name );

		void go_to_end();

		int get_current_token_id();

		vector<string> get_tokens();

		int get_number_of_tokens();

		bool is_end();

		void operator++();

		void operator++( int t );

		void operator+=( int step );

	private:
    //split file specified by filePath arg into tokens, which are words and delimiters as well
		static void tokenize_file( const string filePath, vector<string>& tokens );

    //if rest of line is comment return true
		static bool split_by_delimiters( string s, vector<string>& tokens );

	};

}//IOMANAGMENT
}//SEM

#endif //_FileProcessor_H
