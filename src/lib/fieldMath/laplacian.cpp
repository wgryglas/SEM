
#include "laplacian.h"

#include <fstream>


//namespace SEM
//{
//namespace FIELDMATH
//{

	/*
	Equation laplacian
	(
		SEM::FIELDS::ScalarField& T
	)
	{
		cout<<info<<"building stiffnes matrix ...."<<endl;

		Equation equation = Equation(T);

		//equation.matrix.resize(T.Tacc.size(),T.Tacc.size());//resize and set to 0
		equation.rhsVector.setZero(T.Tacc.size());//resize and set to 0

		//temp element matrix
		vector<vector<vector<vector<double>>>> eStiff;


		//build elements of stiff matrix, and initiali allocate it to global coefficients
		int row;
		int col;
		for( vector<RealElement*>:: iterator e=T.mesh.rEnts.begin(); e!=T.mesh.rEnts.end();e++ )
		{
			(*e)->sEnt.generateStiff_Matrix(eStiff,1,(*e)->H_Matrix);

			// Apply stiffnes element matrix coeffs to global matrix.
			//Take also to evaluation dirichlet BC information(don't write) coresponding equations,
			//evaluate column*dirValue and place it to rhs vector of nondirichlet BC rows
			for(int p=0;p<eStiff.size();p++)
			{
				for(int q=0;q<eStiff[p].size();q++)
				{
					for(int i=0;i<eStiff[p][q].size();i++)
					{
						for(int j=0;j<eStiff[p][q][i].size();j++)
						{
							row = (*e)->mask[p][q];
							col = (*e)->mask[i][j];

							if(!T.dirichletMask[row])
							{
								//Write only those equation(rows), which don't corespond to dirichlet BC in row and col
								if(T.dirichletMask[col])//<-- column coresponds to dirichlet node
								{
									// Place product of multiplication coeff*bound_value to rhs
									equation.rhsVector[row]-=eStiff[p][q][i][j]*T.Tacc[col];
								}
								else //<--- neither row nor column corespond to dirichlet
								{
									// place matrix coefficient. Matrix coeffs can corespond to row,col many times. Sp. matrix builder will add coffs if row,col is def. again
									equation.coeffs.push_back( SMCoeff(row, col, eStiff[p][q][i][j]));
								}
							}
						}
					}
				}
			}

		}

		//Apply integral from Neumann BC
		vector<double>* nodeNeuman;
		int aa=0;
		for(map<string, boundaryType>:: iterator bound=T.bType.begin(); bound!=T.bType.end(); bound++)
		{
			if(bound->second == neumann)
			{

				for(vector<boundaryEdge>::iterator bEdge = T.mesh.boundaryInfo[bound->first].begin(); bEdge!= T.mesh.boundaryInfo[bound->first].end(); bEdge++)
				{
					aa++;
					cout<<aa<<endl;
					nodeNeuman= &(T.neumanValues[&(*bEdge)]);
					for(int i=0;i<bEdge->nodesId.size();i++)
					{
						equation.rhsVector[bEdge->nodesId[i]]-=bEdge->neumanBCPreValue[i] * (*nodeNeuman)[i];
					}
				}
			}
		}

		return equation;
	};
	*/



//	Equation& laplacian(SEM::FIELDS::Property& prop,SEM::FIELDS::ScalarField& T)
//	{
//		Equation equation = Equation(T);

//		switch(prop.getType())
//		{
//		case constantScalar:

//			laplaceConstScalarCoeff(prop.getValue(),T,equation);

//			break;

//		case constantTensor:

//			break;

//		case scalarField:

//			break;

//		case tensorField:

//			break;
//		}


//		return equation;
//	}

//	void laplaceConstScalarCoeff(double a, SEM::FIELDS::ScalarField& T, Equation& equation)
//	{
//		cout<<info<<"building stiffnes matrix ...."<<endl;

//		equation.rhsVector.setZero(T.Tacc.size());//resize and set to 0
//	 	vector<vector<vector<vector<double>>>> eStiff;

//		//build elements of stiff matrix, and initiali allocate it to global coefficients
//		int row;
//		int col;
//		for( vector<RealElement*>:: iterator e=T.mesh.rEnts.begin(); e!=T.mesh.rEnts.end();e++ )
//		{
//			(*e)->sEnt.generateStiff_Matrix(eStiff,a,(*e)->H_Matrix);

//			// Apply stiffnes element matrix coeffs to global matrix.
//			//Take also to evaluation dirichlet BC information(don't write) coresponding equations,
//			//evaluate column*dirValue and place it to rhs vector of nondirichlet BC rows
//			for(int p=0;p<eStiff.size();p++)
//			{
//				for(int q=0;q<eStiff[p].size();q++)
//				{
//					for(int i=0;i<eStiff[p][q].size();i++)
//					{
//						for(int j=0;j<eStiff[p][q][i].size();j++)
//						{
//							row = (*e)->mask[p][q];
//							col = (*e)->mask[i][j];

//							if(!T.dirichletMask[row])
//							{
//								//Write only those equation(rows), which don't corespond to dirichlet BC in row and col
//								if(T.dirichletMask[col])//<-- column coresponds to dirichlet node
//								{
//									// Place product of multiplication coeff*bound_value to rhs
//									equation.rhsVector[row]-=eStiff[p][q][i][j]*T.Tacc[col];
//								}
//								else //<--- neither row nor column corespond to dirichlet
//								{
//									// place matrix coefficient. Matrix coeffs can corespond to row,col many times. Sp. matrix builder will add coffs if row,col is def. again
//                                    equation.coeffs.push_back( MCoeff(row, col, eStiff[p][q][i][j]));
//								}
//							}
//						}
//					}
//				}
//			}

//		}

//		//Apply integral from Neumann BC
//		vector<double>* nodeNeuman;
//		for(map<string, boundaryType>:: iterator bound=T.bType.begin(); bound!=T.bType.end(); bound++)
//		{
//			if(bound->second == neumann)
//			{

//				for(vector<boundaryEdge>::iterator bEdge = T.mesh.boundaryInfo[bound->first].begin(); bEdge!= T.mesh.boundaryInfo[bound->first].end(); bEdge++)
//				{
//					nodeNeuman= &(T.neumanValues[&(*bEdge)]);
//					for(int i=0;i<bEdge->nodesId.size();i++)
//					{
//						equation.rhsVector[bEdge->nodesId[i]]-=bEdge->neumanBCPreValue[i] *a*(*nodeNeuman)[i];
//					}
//				}
//			}
//		}

//	};



//}
//}
