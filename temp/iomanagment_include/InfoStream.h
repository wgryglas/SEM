

#ifndef _InfoStream_H_
#define _InfoStream_H_

#include<iostream>


using namespace std;
namespace SEM
{
namespace IOMANAGMENT
{

	/*class InfoStream
	{

	private:
		ostream& os;

	public:
		InfoStream
		(
		);

		template <typename T>
		ostream& operator<<
		(
			const T &input
		)
		{

			os<<"------------------ SEM ----------------"<<endl<<"Info: "<<input;
			return this->os;
		 }

		ostream& operator<<
		(
			std::ostream& (*pf) (std::ostream&)
		)
		{
			os<<pf;
			return this->os;
		}
	};*/

	////object Info
	//extern InfoStream Info;

	//manipulator:
	ostream& endInfo(ostream&);

	ostream& info(ostream& o);


}//IOMANAGMENT
}//SEM

#endif _InfoStream_H_
