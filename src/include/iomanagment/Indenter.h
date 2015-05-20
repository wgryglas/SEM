#ifndef INDENTER_H_
#define INDENTER_H_

#include <sstream>
#include <iostream>
#include <vector>
#include <algorithm>


////////////////////////////////////////////
///  \class Indenter
/// ---------------------------------------
/// Utility for indenting output stream
/// This class as input take some text stream
/// and when find special character then
/// after it increments indenting next lines.
/// As a result retruns new indented stream
///////////////////////////////////////////

namespace SEM { namespace iomanagment {

class Indenter
{
private:
	int value_;
	int step_;

    std::vector<char> indentingChars;
    std::vector<char> unindentingChars;

public:
	Indenter(int val=0, int step=4);

	template<typename S>
	void addIndentation(S& stream);

	template<typename sOut>
    void indetStream(std::stringstream& stream, sOut& out);

	Indenter& operator++();
	Indenter& operator--();
	Indenter& operator++(int unused);
	Indenter& operator--(int unused);

private:
	inline bool isIndenting(const char& sign) { return containsChar(indentingChars,sign); }
	inline bool isUnindenting(const char& sign) { return containsChar(unindentingChars,sign); }

    static bool containsChar(const std::vector<char>& v, const char& val);
};

inline Indenter::Indenter(int val, int step): value_(val), step_(step)
{
	indentingChars.push_back('{');
	indentingChars.push_back('(');

	unindentingChars.push_back('}');
	unindentingChars.push_back(')');
}

inline bool Indenter::containsChar(const std::vector<char>& v, const char& val)
{
    return std::find(v.begin(),v.end(),val) != v.end();
}

template<typename sOut>
inline void Indenter::indetStream(std::stringstream& stream, sOut& out)
{
    using namespace std;
	value_ = 0;
	stream.seekp(0,iostream::beg);

    if(!stream.good())
        return;

    char curr,next;
    stream.get(next);

    while(stream.good())
	{
		curr = next;
        out<<curr;
        next = stream.get();

		if(isIndenting(curr))
		{
			(*this)++;
			continue;
		}
		else if(isUnindenting(curr))
		{
			(*this)--;
			continue;
		}
		else if(curr=='\n')
		{
			if(isUnindenting(next))//If next one is deindenting - indent earelier
			{
				(*this)--;
				addIndentation<sOut>(out);
				(*this)++;
			}
			else
			{
				addIndentation<sOut>(out);
			}
		}
	}

}

template<typename S>
inline void Indenter::addIndentation(S& stream)
{
    stream<<std::setw(value_*step_);
}

inline Indenter& Indenter::operator ++()
{
	value_++;
	return *this;
}

inline Indenter& Indenter::operator --()
{
	value_--;
	if(value_<0)
		value_=0;

	return *this;
}

inline Indenter& Indenter::operator++(int unused)
{
	value_++;
	return *this;
}

inline Indenter& Indenter::operator--(int unused)
{
	value_--;
	if(value_<0)
		value_=0;

	return *this;
}


} //iomanagment
} //SEM




#endif /* INDENTER_H_ */
