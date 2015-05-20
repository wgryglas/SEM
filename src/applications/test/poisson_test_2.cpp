#include <iostream>
#include <cstdlib>
#include <array>

#include "SE_QLGL.h"
#include "Quad.h"
#include "Mesh.h"
#include<fstream>
#include<regex>
#include"FileProcessor.h"
#include"IOFile.h"
#include "InfoStream.h"
#include "ScalarField.h"
#include "laplacian.h"
#include "Solver.h"
#include "Time.h"
#include "ddt.h"
#include "Case.h"
#include "f.h"
#include "Equation.h"
#include "TransportProperties.h"

using namespace SEM::FIELDMATH;
using namespace SEM::FIELDS;
using namespace SEM::IOStream;
using namespace SEM::MESH;
using namespace SEM;
using namespace std;

int main(int argc, char* args[])
{
	std::cout<<"=========== POISSON SEM APPLICATION =============="<<endl;

	#include "setRootCase.h"
	#include "createTime.h"
	#include "createMesh.h"
	#include "createFileds.h"
	#include "readProperties.h"
	
	//Test rhs field -> equal to ddt(sol)-laplacian(sol), where sol-analitical solution;
	ScalarField rhsField(IOFile(rootCase+"rhsPoisson.sem"), mesh);
	
	//Time Loop
	double x,y;
	while(!runTime.end())
	{
		cout<<"run time: "<<runTime<<endl;

		//Evaluate rhs of known solution
		for(int i=0;i<mesh.nSNodes; i++)
		{
			x = mesh.sCoords[i][0];
			y = mesh.sCoords[i][1];  
			rhsField.Tacc[i] =( x*y*(x-1)*(y-1) - 2*a.getValue()*(y*(y - 1) + x*(x -1))*Time::time() );
		}
		
		//Apply neumann BC values on apropriate edge(du/dx):
		vector<double> temp;
		T.neumanValues.clear();
		for(map<string,boundaryType>::iterator b=T.bType.begin(); b != T.bType.end();b++)
		{
			if(b->second==neumann)
			{
				for(vector<boundaryEdge>::iterator bEdge=mesh.boundaryInfo[b->first].begin();bEdge!=mesh.boundaryInfo[b->first].end();bEdge++ )
				{
					temp.resize(bEdge->nodesId.size());
					for(int node=0; node< bEdge->nodesId.size(); node++)
					{
						x = mesh.sCoords[bEdge->nodesId[node]][0];
						y = mesh.sCoords[bEdge->nodesId[node]][1];
						temp[node]= ( (2*x-1)*y*(y-1)*Time::time() );
					}
					
					T.neumanValues.insert(pair< boundaryEdge*, vector<double> >(&(*bEdge),temp));
				}
			}
		}

		T = solve( ddt(T) - laplacian(a,T)=f(rhsField) ); 
		
		T.writeList();

		runTime++;
	}

	std::cout<<info<<"SOLUTION IS DONE. POISSON APPLICATION EXITING..."<<endInfo;
	
	return 0;
}




