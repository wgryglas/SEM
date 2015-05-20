
#ifndef _laplacian_H_
#define _laplacian_H_

#include"Equation.h"
#include"ScalarField.h"
#include "TransportProperties.h"

using namespace SEM::FIELDS;

namespace SEM
{
namespace FIELDMATH
{


	Equation& laplacian(SEM::FIELDS::Property& prop,SEM::FIELDS::ScalarField& T);

	void laplaceConstScalarCoeff(double a, SEM::FIELDS::ScalarField& T,Equation& eq);



	/*
	void assignElementStiffMatrix
	(
		vector<vector<int>>& mask,
		vector<vector<vector<vector<double>>>>& eStiff,
		vector<vector<vector<double>>>& esStiff
	)
	{
		for(int alpha=0;alpha<eStiff.size();alpha++)
		{
			for(int betha=0;betha<eStiff[alpha].size();betha++)
			{
				if(esStiff[mask[alpha][betha]].size()==0)
				{
					esStiff[mask[alpha][betha]] = vector<vector<double>>(mask.size());
					for(int k=0;k<mask.size();k++)
						esStiff[mask[alpha][betha]][k] = vector<double>(mask[k].size());
				}
				else
				{
					cout<<esStiff[mask[alpha][betha]].size()<<"\t";
				}

				for(int i=0;i<eStiff[alpha][betha].size();i++)
				{
					for(int j=0;j<eStiff[alpha][betha][i].size();j++)
					{

						esStiff[mask[alpha][betha]][i][j] += eStiff[alpha][betha][i][j];

					}
				}
			}
		}

	};
	*/


}
}



#endif _laplacian_H_
