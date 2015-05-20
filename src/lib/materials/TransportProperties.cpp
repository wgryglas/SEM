/*
#include "TransportProperties.h"

namespace SEM
{
namespace FIELDS
{

	propertyType str2propertyType(string s)
	{
		if(s == "constScalar")
			return constantScalar;

		if(s == "constTensor")
			return constantTensor;

		if(s == "scalarField")
			return scalarField;

		if(s == "tensorField")
			return tensorField;

		return constantScalar;
	}

	Property::Property(string name, double value):name(name), constValue(value), type(constantScalar)
	{
	};

	Property::Property(): name(""), constValue(1), type(constantScalar)
	{
	};

	double Property::getValue()
	{
		return constValue;
	}

	propertyType Property::getType()
	{
		return type;
	};

	TransportProperties::TransportProperties(SEM::IOStream::IOFile& file)
	{
		SEM::IOStream::IOObject& inputData = file.getReadObject();

		double value;

		for(map<string,SEM::IOStream::IOSubObject*>:: iterator obj=inputData.getIOSubObjects().begin(); obj!=inputData.getIOSubObjects().end(); obj++)
		{

			switch( str2propertyType(obj->second->getFiledAt("type")) )
			{
			case constantScalar:
				value = stod(obj->second->getFiledAt("value"));
				properites.insert(pair<string,Property>(obj->second->getName(), Property(obj->second->getName(),value) ) );
				break;

			case constantTensor:

				break;

			case scalarField:

				break;

			case tensorField:

				break;
			}
		}
	};

	 Property& TransportProperties::getProperty(string name)
	 {
		 return properites[name];
	 };
}
}*/
