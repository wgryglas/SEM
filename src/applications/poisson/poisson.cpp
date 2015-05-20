// #include <iostream>
// #include <cstdlib>
// #include <array>
// 
// #include "SE_QLGL.h"
// #include "Quad.h"
// #include "Mesh.h"
// #include<fstream>
// #include<regex>
// #include"FileProcessor.h"
// #include"IOFile.h"
// #include "InfoStream.h"
// #include "ScalarField.h"
// #include "laplacian.h"
// #include "Solver.h"
// #include "Time.h"
// #include "ddt.h"
// #include "Case.h"
// #include "TransportProperties.h"
// 
// 
// using namespace SEM::FIELDMATH;
// using namespace SEM::FIELDS;
// using namespace SEM::IOStream;
// using namespace SEM::MESH;
// using namespace SEM;
// using namespace std;
// 
// int main(int argc, char* args[])
// {
// 	std::cout<<"=========== POISSON SEM APPLICATION =============="<<endl;
// 
// 	#include "setRootCase.h"
// 	#include "createTime.h"
// 	#include "createMesh.h"
// 	#include "createFileds.h"
// 	#include "readProperties.h"
// 
// 	//Time Loop
// 	while(!runTime.end())
// 	{
// 		cout<<"run time: "<<runTime<<endl;
// 
// 		T = solve( ddt(T) - laplacian(a,T) );
// 		T.writeList();
// 
// 		runTime++;
// 	}
// 
// 	std::cout<<info<<"SOLUTION IS DONE. POISSON APPLICATION EXITING..."<<endInfo;
// 
// 
// 
// 	return 0;
// }
// 
// 
// 
// 
