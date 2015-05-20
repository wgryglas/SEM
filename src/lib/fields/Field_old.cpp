/*
#include "Field.h"


namespace SEM
{
namespace FIELDS
{

	Field::Field
	(
		string name,
		Mesh& mesh
	): mesh(mesh), name(name)
	{
	};

	inputType str2InputType
	(
		string s
	)
	{
		if(s=="uniform")
			return inputType::uniform;

		if(s=="nodalValue")
			return inputType::nodalValue;

		return inputType::nodalValue;
	}

	boundaryType str2BoundaryType
	(
		string s
	)
	{
		if(s=="dirichlet")
			return dirichlet;

		if(s=="neumann")
			return neumann;

		return boundaryType::dirichlet;
	}


	Field::~Field()
	{
	};

}
}
*/
