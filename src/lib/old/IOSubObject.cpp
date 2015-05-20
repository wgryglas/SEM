#include "IOSubObject.h"

#include<iostream>

namespace SEM
{
namespace IOStream
{
	IOSubObject:: IOSubObject
	(
		string name,
		IOSubObject* parent
	): name(name), parent(parent)
	{
	};

	string IOSubObject::getName
	(
	)
	{
		return name;
	}

	void IOSubObject::addField
	(
		string name, 
		string value
	)
	{
		fields.insert(pair<string,string>(name,value));
	};

	void IOSubObject::addScalarFiled
	(
		string name,
		vector<string> value
	)
	{
		scalarFields.insert(pair<string,vector<string>>(name, value));	
	};

	void IOSubObject::addVectorFiled
	(
		string name,
		vector<array<string,2>> value
	)
	{
		vectorFields.insert(pair<string,vector<array<string,2> > >(name, value));
	};
		
	void IOSubObject::addTensorFiled
	(
		string name,
		vector<array<string,4>> value
	)
	{
		tensorFields.insert(pair<string,vector<array<string,4> > >(name, value));
	};

	void IOSubObject::addTableFiled
	(
		string name,
		vector<vector<string>> value
	)
	{
		tableFields.insert(pair< string, vector< vector<string> > >(name, value));
	};


	void IOSubObject::addIOObject
	(
		IOSubObject* value
	)
	{
		iOSubObjects.insert(pair<string,IOSubObject*>(value->getName(),value));
	};

	bool IOSubObject::isSubObject
	(
	)
	{
		return true;
	};

	IOSubObject* IOSubObject::getParentAdress
	(
	)
	{
		return parent;
	}

	bool IOSubObject::haveObj
	(
		string name
	)
	{
		if(fields.find(name)==fields.end())
		{
			return false;
		}
		else
		{
			return true;
		}	
	};

	bool IOSubObject::haveField
	(
		string name
	)
	{
		if(fields.find(name)==fields.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	bool IOSubObject::haveScalarField
	(
		string name
	)
	{
		if(scalarFields.find(name)==scalarFields.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	bool IOSubObject::haveVectorField
	(
		string name
	)
	{
		if(vectorFields.find(name)==vectorFields.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	bool IOSubObject::haveTensorField
	(
		string name
	)
	{
		if(tensorFields.find(name)==tensorFields.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	};

	bool IOSubObject::haveTableField
	(
		string name
	)
	{
		if(tableFields.find(name)==tableFields.end())
		{
			return false;
		}
		else
		{
			return true;
		}
	};



	string& IOSubObject::getFiledAt
	(
		string name
	)
	{
		return fields.at(name);
	}

	vector<string>& IOSubObject::getScalarFieldAt
	(
		string name
	)
	{
		return scalarFields.at(name);
	};

	vector<array<string,2> >& IOSubObject::getVectorFieldAt
	(
		string name
	)
	{
		return vectorFields.at(name);
	};

	vector<array<string,4> >& IOSubObject::getTensorFieldAt
	(
		string name
	)
	{
		return tensorFields.at(name);
	};

	vector<vector<string> >& IOSubObject::getTableFieldAt
	(
		string name
	)
	{
		return tableFields.at(name);
	};


	IOSubObject* IOSubObject::getObjAt
	(
		string name
	)
	{
		return iOSubObjects.at(name);
	}

	map<string,IOSubObject*>& IOSubObject::getIOSubObjects
	(
	)
	{
		return iOSubObjects;
	}

	IOSubObject::~IOSubObject
	(
	)
	{
	};


	ostream& operator<<
	(
		ostream& o,
		IOSubObject& io
	)
	{
		o<<"---------- "<<io.getName()<<" ------------------"<<endl;
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
		o<<"------------ END OF "<<io.getName()<<" -----------------"<<endl;
		return o;
	}


}
}