
#ifndef _IOSubOBject_H_
#define _IOSubOBject_H_

#include<string>
#include<map>
#include<vector>
#include<array>


using namespace std;

namespace SEM
{
namespace IOStream
{
	class IOSubObject
	{

	public:
		IOSubObject* parent;

		map<string, string> fields;
		map<string, vector<string> > scalarFields;
		map<string, vector< string[2] > > vectorFields;
		map<string, vector< string[4] > > tensorFields;
		map<string, vector< vector<string> > >  tableFields;

		map<string,IOSubObject*> iOSubObjects;


	protected:
		string name;

	public:

		IOSubObject
		(
		string name,
		IOSubObject* parent
		);

		~IOSubObject
		(
		);

		void addField
		(
			string name,
			string value
		);

		void addScalarFiled
		(
			string name,
			vector<string> value
		);

		void addVectorFiled
		(
			string name,
			vector<array<string,2>> value
		);

		void addTensorFiled
		(
			string name,
			vector<array<string,4>> value
		);

		void addTableFiled
		(
			string name,
			vector<vector<string>> value
		);

		void addIOObject
		(
			IOSubObject* value
		);

		string& getFiledAt
		(
			string name
		);

		vector<string>& getScalarFieldAt
		(
			string name
		);

		vector<array<string,2> >&  getVectorFieldAt
		(
			string name
		);

		vector<array<string,4> >& getTensorFieldAt
		(
			string name
		);

		vector<vector<string> >& getTableFieldAt
		(
			string name
		);

		IOSubObject* getObjAt
		(
			string name
		);

		map<string,IOSubObject*>& getIOSubObjects
		(
		);

		bool haveField
		(
			string name
		);

		bool haveScalarField
		(
			string name
		);

		bool haveVectorField
		(
			string name
		);

		bool haveTensorField
		(
			string name
		);

		bool haveTableField
		(
			string name
		);


		bool haveObj
		(
			string name
		);




		string getName
		(
		);

		IOSubObject* getParentAdress
		(
		);


		virtual bool isSubObject
		(
		);


	};

	ostream& operator<<(ostream& o, IOSubObject& io);


}
}


#endif //_IOSubOBject_H_
