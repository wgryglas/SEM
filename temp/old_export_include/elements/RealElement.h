
#ifndef _RealElement_H_
#define _RealElement_H_

#include<iostream>
#include <array>

#include "SE_QLGL.h"



namespace SEM
{
namespace MESH
{
	class RealElement
	{
	public:
		int id;
		
		vector<array<double,2>> nodes;
		SpectralElement& sEnt;
		vector<vector<int>> mask;
		vector<vector<array<double,6>>> H_Matrix;
		vector<vector<double>> M_Matrix;
		

		/*RealElement
		(
		);*/

		/*RealElement
		(
			RealElement& e
		);*/
		
		void operator=
		(
			const RealElement& a
		);

		RealElement
		(
			vector<array<double,2>> nodes,
			SpectralElement& sEnt,
			vector<vector<int>> mask
		);

	};
}
}


#endif _RealElement_H_