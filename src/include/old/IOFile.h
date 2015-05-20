
#ifndef _IOFile_H_
#define _IOFile_H_

#include<string>

#include "IOObject.h"
#include "FileProcessor.h"
#include "IOTokens.h"

namespace SEM
{
namespace IOStream
{
	class IOFile
	{
	private:
		IOSubObject* currIOObj;
		FileProcessor fProcessor;


		IOObject result = *(new IOObject("SEMFile"));//("SEMFile");

		///vector

	public:
		IOFile
		(
			string filePath
		);

		IOObject& getReadObject
		(
		);


		~IOFile
		(
		);




	};

}
}

#endif //_IOFile_H_
