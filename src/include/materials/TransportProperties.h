
// #ifndef _TransportProperties_H_
// #define _TransportProperties_H_
// 
// #include <string>
// #include <map>
// 
// #include "IOFile.h"
// #include "FileProcessor.h"
// #include "IOObject.h"
// 
// 
// using namespace std;
// 
// namespace SEM
// {
// namespace FIELDS
// {
// 
// 	enum propertyType
// 	{
// 		constantScalar,
// 		constantTensor,
// 		scalarField,
// 		tensorField
// 	};
// 
// 	propertyType str2propertyType(string s);
// 
// 	class Property
// 	{
// 	private:
// 		string name;
// 		propertyType type;
// 		double constValue;
// 
// 	public :
// 		Property(string name, double value);
// 		Property();
// 
// 		double getValue();
// 		propertyType getType();
// 
// 	};
// 
// 	class TransportProperties
// 	{
// 	private:
// 		map<string,Property> properites;
// 
// 	public:
// 		TransportProperties(SEM::IOStream::IOFile& file);
// 		TransportProperties(){};
// 
// 		Property& getProperty(string name);
// 
// 
// 	};
// 
// 
// 
// }
// }
// 
// 
// #endif _TransportProperties_H_
