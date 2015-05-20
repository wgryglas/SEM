
#include "IOFile.h"

namespace SEM
{
namespace IOStream
{
	IOFile::IOFile
	(
		string filePath
	): fProcessor(filePath), result(IOObject("SEMFile"))
	{
		currIOObj = &result;

		vector<IOToken*> tokens;
		tokens.push_back(new IOOpenObject());
		tokens.push_back(new IOCloseObject());
		tokens.push_back(new IOEndOfFiled());
		tokens.push_back(new IOList());

		vector<IOToken*>:: iterator i;
		bool searchTokens;

		while(!fProcessor.end())
		{
			searchTokens = false;
			for(i=tokens.begin();i!=tokens.end() && !searchTokens ;i++)
			{
				searchTokens = (*i)->act(
										  fProcessor,
										  &currIOObj
									    );
			}
			fProcessor.goForeward();
		}
	};



	IOObject& IOFile:: getReadObject
	(
	)
	{
		return result;
	}


	IOFile::~IOFile
	(
	)
	{
	};


}
}
