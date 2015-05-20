
#ifndef _Field_H_
#define _Field_H_

#include "Mesh.h"

using namespace SEM::MESH;

namespace SEM
{
namespace FIELDS
{
	class Field
	{
	public:
		
		Mesh& mesh;
		string name;

		Field
		(
			string name,
			Mesh& mesh
		);
		

		virtual ~Field
		(
		);

	protected:
		
	};


	enum inputType
	{
		uniform,
		nodalValue
	};

	inputType str2InputType
	(
		string s
	);

	enum boundaryType
	{
		dirichlet,
		neumann
	};

	boundaryType str2BoundaryType
	(
		string s
	);


}
}






#endif _Field_H_