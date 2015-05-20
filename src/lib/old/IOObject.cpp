
#include "IOObject.h"

namespace SEM
{
namespace IOStream
{

	IOObject::IOObject
	(
	): IOSubObject("", 0)
	{
	};

	IOObject::IOObject
	(
		string name
	): IOSubObject(name, 0)
	{
	};


	void IOObject::changeName
	(
		string name
	)
	{
		this->name = name;
	}

	bool IOObject::isSubObject
	(
	)
	{
		return false;
	};

	IOObject::~IOObject
	(
	)
	{
	};

	ostream& operator<<
	(
		ostream& o,
		IOObject io
	)
	{
	    /*
		o<<"============ "o<<io.getName()<<" ================="<<endl;
		o<<endl<<"fields::"<<endl;
		for(map<string,string>::iterator i=io.fields.begin();i!=io.fields.end();i++)
		{
			o<<"\t"<<i->first<<"="<<i->second<<endl;
		}
		o<<endl<<"scalar fields::"<<endl;
		for(map<string,vector<string>>::iterator i=io.scalarFields.begin();i!=io.scalarFields.end();i++)
		{
			o<<"\t"<<i->first<<"=";
			for(vector<string>::iterator k=i->second.begin();k!=i->second.end(); k++)
					o<<*k<<",";
			o<<endl;
		}
		o<<endl<<"vector fields::"<<endl;
		for(map<string,vector<array<string,2>>>::iterator i=io.vectorFields.begin();i!=io.vectorFields.end();i++)
		{
			o<<"\t"<<i->first<<"=";
			for(vector<array<string,2>>::iterator k=i->second.begin();k!=i->second.end(); k++)
					o<<"("<<(*k)[0]<<" "<<(*k)[1]<<"), ";
			o<<endl;
		}

		o<<endl<<"tensor fields::"<<endl;
		for(map<string,vector<array<string,4>>>::iterator i=io.tensorFields.begin();i!=io.tensorFields.end();i++)
		{
			o<<"\t"<<i->first<<"=";
			for(vector<array<string,4>>::iterator k=i->second.begin();k!=i->second.end(); k++)
					o<<"("<<(*k)[0]<<" "<<(*k)[1]<<" "<<(*k)[2]<<" "<<(*k)[3]<<"), ";
			o<<endl;
		}

		o<<endl<<"table fields::"<<endl;
		for(map<string,vector<vector<string>>>::iterator i=io.tableFields.begin();i!=io.tableFields.end();i++)
		{
			o<<"\t"<<i->first<<"=";
			for(vector<vector<string>>::iterator k=i->second.begin();k!=i->second.end(); k++)
			{
					o<<"(";
					for(vector<string>:: iterator j=k->begin();j!=k->end();j++)
						o<<*j<<" ";

					o<<"), ";
			}
			o<<endl;
		}


		o<<endl<<"sub-objects::"<<endl;
		for(map<string,IOSubObject*>::iterator i=io.iOSubObjects.begin();i!=io.iOSubObjects.end(); i++)
		{
			o<<*(i->second)<<endl;
			o<<endl;
		}
		o<<"============ END OF "<<io.getName()<<" ================="<<endl;

		*/

		return o;
	}


}
}
