//#include "f.h"

//using namespace SEM::FIELDS;

//namespace SEM
//{
//namespace FIELDMATH
//{
//	Eigen::VectorXd& f(ScalarField U)
//	{
//		//Init equation:
//		Equation equation(U);
//		equation.matrix.resize(U.Tacc.size(),U.Tacc.size());
//		equation.rhsVector.setZero(U.Tacc.size());

//		Vector* vec = new Vector(U.Tacc.size());
//		vec->setZero();

//		int row_col;
		
//		for( vector<RealElement*>:: iterator e=U.mesh.rEnts.begin(); e!=U.mesh.rEnts.end();e++ )
//		{
//			//Allocate eMass in global matrix:
//			for(int p=0;p<(*e)->M_Matrix.size();p++ )
//			{
//				for(int q=0;q<(*e)->M_Matrix[p].size();q++)
//				{
//					row_col = (*e)->mask[p][q];
					
//					if(!U.dirichletMask[row_col])
//						(*vec)[row_col] += U.Tacc[row_col] * ((*e)->M_Matrix[p][q]);
//				}
//			}
//		}

//		return *vec;

//	};
//}
//}

