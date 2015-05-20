
#include "file_processor.h"

namespace SEM
{
namespace iomanagment
{
  const boost::char_separator<char> FileProcessor::delimiters_("","{}[];,()");
  const boost::regex FileProcessor::line_comment_delimiter_("(//)");
  typedef boost::tokenizer< boost::char_separator<char> > boost_tokenizer;

    FileProcessor::FileProcessor( const string filePath ): itr_(0)
	{
		FileProcessor::tokenize_file(filePath, tokens_);
    }

	void FileProcessor::go_foreward()
	{
		itr_++;
	}

    void FileProcessor::go_foreward( int step )
	{
		itr_++;
	}

	void FileProcessor::operator++()
	{
		++itr_;
	}

	void FileProcessor::operator++(int t)
	{
		itr_++;
	}

	void FileProcessor::operator+=( int step )
	{
		itr_+=step;
    }

    std::string FileProcessor::get_token()
	{
        if(itr_<tokens_.size())
          return tokens_[itr_];
        else
          return "";
    }

    std::string FileProcessor::get_token(int offset)
	{
        if(itr_+offset<tokens_.size())
        {
          return tokens_[itr_+offset];
        }
        else
        {
          return "";
        }
	}

	void FileProcessor::go_to_begin()
	{
		itr_ = 0;
    }

	bool FileProcessor::go_to( int place )
	{
		if(place<tokens_.size())
		{
			itr_ = place;
			return true;
		}
		else
		{
			return false;
		}
    }

	bool FileProcessor::go_to( const string token_name )
	{
		int temp=-1;
		bool search;
        std::vector<std::string>::iterator i=tokens_.begin();

		while(i < tokens_.end() && search)
		{
			if(*i==token_name)
			{
				itr_= temp;
				search=false;
			}

			i++;
			temp++;
		}

		if(search)
		{
			return false;
		}
		else
		{
			return true;
		}
    }

	void FileProcessor::go_to_end()
	{
		itr_ = tokens_.size()-1;
    }


	int FileProcessor::get_current_token_id()
	{
		return itr_;
	}

    std::vector<std::string> FileProcessor::get_tokens()
	{
		return tokens_;
	}

	int FileProcessor::get_number_of_tokens()
	{
		return tokens_.size();
	}

	bool FileProcessor::is_end()
	{
		if(itr_==tokens_.size())
			return true;
		else
			return false;
	}


	void FileProcessor::tokenize_file( const string filePath, vector<string>& tokens_ )
	{
        std::ifstream file;
        file.open(filePath.data());
        std::string temp;
		while(file>>temp)
		{
			if(FileProcessor::split_by_delimiters(temp,tokens_)) //when line comment delimiter occurs ignor rest input from this line
				 file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
		}

		file.close();
    }

    bool FileProcessor::split_by_delimiters( string s, vector<string>& tokens)
	{
        //Check if comment line delimiter exist in input string:
        bool comment_exist = false;
        boost::sregex_token_iterator i_c(s.begin(), s.end(), line_comment_delimiter_, -1);
        boost::sregex_token_iterator j_c(i_c);
        boost::sregex_token_iterator k;
        if(++j_c!=k)
        {
          s = *i_c; // parse only uncommented string part
          comment_exist=true;
        }

        //Split string with delimiters, keep delimiters as tokens_
        boost_tokenizer tok(s,delimiters_);
        for(boost_tokenizer::iterator itr=tok.begin(); itr!=tok.end(); ++itr)
        {
          tokens.push_back(*itr);
        }

        return comment_exist;
    }


}//IOMANAGMENT
}//SEM
