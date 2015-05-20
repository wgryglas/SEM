
#ifndef _IOObject_H_
#define _IOObject_H_

#include<string>

#include"IOSubObject.h"

using namespace std;

namespace SEM
{
namespace IOStream
{

	class IOObject: public IOSubObject
	{
		public:
			IOObject
			(
			);

			IOObject
			(
			string name
			);

			void changeName(string name);

			bool isSubObject
			(
			);

			~IOObject
			(
			);

	};


	ostream& operator<<
	(
		ostream& o,
		IOObject io
	);



}
}






#endif //_IOObject_H_
