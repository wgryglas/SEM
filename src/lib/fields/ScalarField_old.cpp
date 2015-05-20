/*
#include "ScalarField.h"
#include "time/Time.h"
#include "iomanagment/case.h"

namespace SEM
{
namespace FIELDS
{
	ScalarField::ScalarField
	(
		IOFile& iof,
		Mesh& mesh
	)	:Field(iof.getReadObject().getFiledAt("name"), mesh), Tacc(vector<double>()), Tprev(vector<double>())
	{

		IOObject& io = iof.getReadObject();

		for(map<string, IOSubObject*>::iterator obj= io.getIOSubObjects().begin(); obj!= io.getIOSubObjects().end(); obj++ )
		{
			if(obj->second->getFiledAt("object")=="internalField")
				applyInteriorInintialCondition(*(obj->second));
		}

		dirichletMask = vector<bool>(mesh.nSNodes,false);
		for(map<string, IOSubObject*>::iterator obj= io.getIOSubObjects().begin(); obj!= io.getIOSubObjects().end(); obj++ )
		{
			if(obj->second->getFiledAt("object")=="boundaryField")
				applyBoundaryCondition(*(obj->second));
		}

	};


	void ScalarField::applyInteriorInintialCondition
	(
		IOSubObject& io
	)
	{
		switch ( str2InputType(io.getFiledAt("valueType")) )
		{
		case uniform:
			Tacc = vector<double>(mesh.nSNodes,stod(io.getFiledAt("value")) );
			Tprev = vector<double>(Tacc);

			break;
		case nodalValue:
			
			break;
		}
    }

	void ScalarField::applyBoundaryCondition
	(
		IOSubObject& io
	)
	{
		inputType inType;

		inType = str2InputType( io.getFiledAt("valueType") );
		switch( str2BoundaryType(io.getFiledAt("boundaryType")) ) 
		{
		case dirichlet:
			{	
				string bName = io.getName();
				bType.insert( pair<string,boundaryType>(bName,dirichlet) );

				switch (inType)
				{
				case uniform:
					{
						double val = stod( io.getFiledAt("value") );
						for(vector<boundaryEdge>::iterator bEdge=mesh.boundaryInfo[bName].begin();bEdge!=mesh.boundaryInfo[bName].end();bEdge++ )
							for(int node=0; node< bEdge->nodesId.size(); node++)
							{
								Tacc[bEdge->nodesId[node]] = val;
								Tprev[bEdge->nodesId[node]] = val;
								dirichletMask[bEdge->nodesId[node]] = true;
							}
					break;
					}
				case nodalValue:
					{
						int itr =0;
						for(vector<boundaryEdge>::iterator bEdge=mesh.boundaryInfo[bName].begin();bEdge!=mesh.boundaryInfo[bName].end();bEdge++ )
							for(int node=0; node< bEdge->nodesId.size(); node++)
							{
								Tacc[bEdge->nodesId[node]] = stod( io.getScalarFieldAt("value")[itr] );
								Tprev[bEdge->nodesId[node]]= Tacc[bEdge->nodesId[node]];
								dirichletMask[bEdge->nodesId[node]] = true;
								itr++;
							}
						break;
					}
				default:
					{
						break;
					}

				break;
				}
			break;
			}
		case neumann:
			{
				string bName = io.getName();
				bType.insert( pair<string,boundaryType>(io.getName(),neumann) );
				vector<double> temp;
				switch (inType)
				{
				case uniform:
					{
						double val = stod( io.getFiledAt("value") );
						for(vector<boundaryEdge>::iterator bEdge=mesh.boundaryInfo[bName].begin();bEdge!=mesh.boundaryInfo[bName].end();bEdge++ )
						{
							temp.resize(bEdge->nodesId.size());
							for(int node=0; node< bEdge->nodesId.size(); node++)
							{
								temp[node]=val;	
							}
							neumanValues.insert(pair< boundaryEdge*, vector<double> >(&(*bEdge),temp));
						}
						break;
					}
				case nodalValue:
					{
						double val;
						int itr =0;
						for(vector<boundaryEdge>::iterator bEdge=mesh.boundaryInfo[bName].begin();bEdge!=mesh.boundaryInfo[bName].end();bEdge++ )
						{
							temp.resize(bEdge->nodesId.size());
							for(int node=0; node< bEdge->nodesId.size(); node++)
							{
								val = stod( io.getScalarFieldAt("value")[itr] );
								itr++;
								temp[node]=val;	
							}
							neumanValues.insert(pair< boundaryEdge*, vector<double>>(&(*bEdge),temp));
						}
					
					break;
					}
				}
			}
		}
		
    }


	ScalarField::~ScalarField
	(
	)
	{
    }


	ostream& ScalarField::displayNodal
	(
		ostream& o
	)
	{
		o<<endl;
		for(vector<double>::iterator i=Tacc.begin();i!=Tacc.end();i++)
			o<<*i<<endl;

		o<<endl;

		return o;
    }
		
	ostream& ScalarField::displayNodalWithCoords
	(
		ostream& o
	)
	{
		o<<endl;
		for(int i=0;i<Tacc.size();i++)
		{
			o<<mesh.sCoords[i][0]<<"\t"<<mesh.sCoords[i][1]<<"\t"<< Tacc[i]<<endl;
		}

		o<<endl;

		return o;
    }

	Eigen::VectorXd ScalarField::getAccSolution
	(
	)
	{
		Eigen::VectorXd vec(Tacc.size());
		for(int i=0;i<Tacc.size();i++)
			vec[i]=Tacc[i];

		return vec;
	};

	Eigen::VectorXd ScalarField::getPrevSolution
	(
	)
	{
		Eigen::VectorXd vec(Tprev.size());
		for(int i=0;i<Tprev.size();i++)
			vec[i]=Tprev[i];

		return vec;
    }


	void ScalarField::writeListWithNodes
	(
	)
	{
		ofstream file;
		file.open(Case::casePath+name+"_"+Time::timeName()+".dat");

		file<<"TITLE = \""<<name<<"_"<<Time::timeName()<<"\""<<endl;
		file<<"Variables = \"X\", \"Y\",\""+name+"\""<<endl;
		//file<<"Zone I="<<Tacc.size()<<"\t J=2 \t K=1"<<endl;
		file<<endl;
		for(int i=0;i<Tacc.size();i++)
		{
			file<<mesh.sCoords[i][0]<<"\t"<<mesh.sCoords[i][1]<<"\t"<<Tacc[i]<<endl;
		}

		file<<endl;

		file.close();
	}


	void ScalarField::writeList
	(
	)
	{
		ofstream file;
		file.open(Case::casePath+name+"_"+Time::timeName()+".semsol");

		file<<endl;
		for(int i=0;i<Tacc.size();i++)
		{
			file<<Tacc[i]<<endl;
		}

		file<<endl;

		file.close();
    }


	ScalarField& ScalarField::operator=
	(
		Eigen::VectorXd& vector
	)
	{
		Tprev = Tacc;

		for(int i=0;i<Tacc.size(); i++)
		{
			Tacc[i] = vector[i];
		}

		return *this;
    }

	ScalarField ScalarField::operator+
	(
		ScalarField& f
	)
	{
		ScalarField sf(*this);

		for(int i=0; i<Tacc.size();i++)
		{
			sf.Tacc[i]+=f.Tacc[i];
			sf.Tprev[i]+=f.Tprev[i];
		}

		return sf;
    }
		
	ScalarField ScalarField::operator-
	(
		ScalarField& f
	)
	{
		ScalarField sf(*this);

		for(int i=0; i<Tacc.size();i++)
		{
			sf.Tacc[i]-=f.Tacc[i];
			sf.Tprev[i]-=f.Tprev[i];
		}

		return sf;
    }



	ostream& operator<<
	(
		ostream& o,
		ScalarField field
	)
	{
		o<<"scalar filed::"<<field.name<<endl;
		o<<"with nodal values:"<<endl;
		for(int n=0;n<field.Tacc.size();n++)
		{
			o<<"\t"<<n<<"\t"<<field.Tacc[n]<<endl;//<<"\t dirichlet mask::"<<field.dirichletMask[n]<<endl;
		}

		return o;
    }
}
}
*/
