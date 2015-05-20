
#ifndef _IOTokens_H_
#define _IOTokens_H_

#include <string>

#include"IOObject.h"
#include"FileProcessor.h"

using namespace std;
namespace SEM
{
namespace IOStream
{
	class IOToken
	{
	public:
		string token;

		IOToken
		(
			string name
		);

		virtual bool act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		)=0;
	};

	class IOOpenObject: public IOToken
	{
	public:
		IOOpenObject
		(
		);

		bool act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		);
	};

	class IOCloseObject: public IOToken
	{
	public:
		IOCloseObject
		(
		);

		bool act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		);
	};


	class IOEndOfFiled: public IOToken
	{
		public:
		IOEndOfFiled
		(
		);

		bool act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		);

	};

	class IOList: public IOToken
	{
	private:
		enum listType
		{
			Scalar,
			Vector,
			Tensor,
			Table,
			None,
		};

		listType str2ListType(string str);

		public:
		IOList
		(
		);

		bool act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		);

	};

}
}
#endif //_IOTokens_H_
