//#include "ddt.h"


//namespace SEM
//{
//namespace FIELDMATH
//{
//	Equation& ddt
//	(
//		ScalarField& U
//	)
//	{
//		//Init equation:
//		Equation equation(U);

//		switch(SolutionControl::timeDiscretization)
//		{
//			case BDF:
//				BDF2(equation);
//				break;

//			default:
//				Euler(equation);
//				break;
//		}

//		return equation;
//	};

//	void Euler
//	(
//		Equation& eq
//	)
//	{
//		//Init rhsVector
//		eq.rhsVector.setZero(eq.field.Tacc.size());

//		//init helper values

//		int row_col;
//		double invDt = 1./SolutionControl::timeStep;
//		double val;

//		// Evaluate and allocate ddt(U) matrix and rhsVector. Take in consideration dirichletBC
//		for( vector<RealElement*>:: iterator e=eq.field.mesh.rEnts.begin(); e!=eq.field.mesh.rEnts.end();e++ )
//		{
//			//Allocate eMass in global matrix:
//			for(int p=0;p<(*e)->M_Matrix.size();p++ )
//			{
//				for(int q=0;q<(*e)->M_Matrix[p].size();q++)
//				{
//					row_col = (*e)->mask[p][q];
//					//Only if node coresponding to equation is not dirichlet type(there whole eq. must have 1 on diag)
//					if(!eq.field.dirichletMask[row_col])
//					{
//						val=invDt*( (*e)->M_Matrix[p][q] );
//						eq.coeffs.push_back( MCoeff(row_col,row_col, val) );
//						eq.rhsVector[row_col] +=val*eq.field.Tacc[row_col];
//					}

//				}
//			}
//		}
//	};


//	void BDF2
//	(
//		Equation& eq
//	)
//	{
//		//Init rhsVector
//		eq.rhsVector.setZero(eq.field.Tacc.size());

//		//init helper values

//		int row_col;
//		double invDt = 1./SolutionControl::timeStep;
//		double val;

//		// Evaluate and allocate ddt(U) matrix and rhsVector. Take in consideration dirichletBC
//		for( vector<RealElement*>:: iterator e=eq.field.mesh.rEnts.begin(); e!=eq.field.mesh.rEnts.end();e++ )
//		{
//			//Allocate eMass in global matrix:
//			for(int p=0;p<(*e)->M_Matrix.size();p++ )
//			{
//				for(int q=0;q<(*e)->M_Matrix[p].size();q++)
//				{
//					row_col = (*e)->mask[p][q];
//					//Only if node coresponding to equation is not dirichlet type(there whole eq. must have 1 on diag)
//					if(!eq.field.dirichletMask[row_col])
//					{
//						val=invDt*( (*e)->M_Matrix[p][q] );

//						eq.coeffs.push_back( MCoeff(row_col,row_col, 1.5*val) );
//						eq.rhsVector[row_col] +=( 2*val*eq.field.Tacc[row_col]-0.5*val*eq.field.Tprev[row_col] );

//					}

//				}
//			}
//		}
//	};



//}
//}




