
#include "IOTokens.h"
namespace SEM
{
	namespace IOStream
	{

		//------- IOTOKEN ------------------------
		IOToken::IOToken
			(
			string name
			):token(name)
		{
		};

		//------- IOOpenObjectToken ---------------
		IOOpenObject::IOOpenObject
			(
			): IOToken("{")
		{
		};

		bool IOOpenObject::act
			(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
			)
		{
			if(token==fProcessor.getToken())
			{
				if( fProcessor.getPreviousToken()!="SEMFile" )
				{
					if( ! ( (*accIOObject)->haveObj(fProcessor.getPreviousToken()) ) )
					{
						IOSubObject* newObj = new IOSubObject(fProcessor.getPreviousToken(), *accIOObject);
						(*accIOObject)->addIOObject( newObj );
						(*accIOObject)=newObj;
					}
					else
					{
						(*accIOObject)=( (*accIOObject)->getObjAt(fProcessor.getPreviousToken()) );
					}
				}

				return true;
			}
			else
			{
				return false;
			}
		};

		//------- IOCloseObjectToken --------------
		IOCloseObject::IOCloseObject
			(
			): IOToken("}")
		{
		};

		bool IOCloseObject::act
			(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
			)
		{
			if(token==fProcessor.getToken())
			{
				if( (*accIOObject)->getName()!="SEMFile" )
					(*accIOObject) = (*accIOObject)->getParentAdress();

				return true;
			}
			else
			{
				return false;
			}
		};

		//------- IOEndOfFiled ---------------------
		IOEndOfFiled::IOEndOfFiled
			(
			): IOToken(";")
		{
		};

		bool IOEndOfFiled::act
			(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
			)
		{

			if( token == fProcessor.getToken() && fProcessor.getPreviousToken()!=")")
			{
				if( !( (*accIOObject)->haveField(fProcessor.getPreviousPreviousToken() ) ) )
				{	
					(*accIOObject)->addField( fProcessor.getPreviousPreviousToken(), fProcessor.getPreviousToken() );
				}
				else
				{
					cout<<"Warrning! field already defined! this value= "<<fProcessor.getPreviousToken()<<" will not be assigned to field "<<fProcessor.getPreviousPreviousToken()<<endl;
					//(*accIOObject)->getFiledAt(fProcessor.getPreviousPreviousToken()).push_back(fProcessor.getPreviousToken());
				}

				return true;
			}
			else
			{
				return false;
			}

		};

		//------- IOList ---------------------
		IOList::IOList
			(
			): IOToken("List")//"List"
		{

		};

		IOList::listType IOList:: str2ListType
			(
			string str
			)
		{
			if(str=="vector")
				return listType::Vector;

			if(str=="scalar")
				return listType::Scalar;

			if(str=="tensor")
				return listType::Tensor;

			if(str=="table")
				return listType::Table;

			return listType::None;
		};

		bool IOList::act
		(
			FileProcessor& fProcessor,
			IOSubObject** accIOObject
		)
		{
			if( token == fProcessor.getToken() )
			{
				string fieldName = fProcessor.getPreviousToken();

				//def. of list: List<type,size>( 1 1, 2 2,...); --> here for vector 
				//fProcessor is at word "List"

				fProcessor+=2; // go from ("List","<") to "type"


				listType listReadType= str2ListType(fProcessor.getToken());


				fProcessor+=2;// go to 2nd arg from ("type",",") to "size"
				int size =  std::stoi(fProcessor.getToken());// lenght of list

				fProcessor+=3; //go to first digit

				switch (listReadType)
				{
				case Scalar:
					{
						vector<string> listScalar(size);

						for(int i=0;i<size;i++)
						{
							listScalar[i]=fProcessor.getToken();

							//go to next token end ommit "," (at end go on ";")
							fProcessor+=2;
						}

						if( !( (*accIOObject)->haveScalarField(fieldName ) ) ) //add new list to new token
						{	
							(*accIOObject)->addScalarFiled(fieldName,listScalar);
						}
						else // fill existing list
						{
							(*accIOObject)->getScalarFieldAt(fieldName).insert((*accIOObject)->getScalarFieldAt(fieldName).end(),listScalar.begin(),listScalar.end());
						}

						break;
					}
				case Vector:
					{
						vector< array<string,2> > listVector(size);

						for(int i=0;i<size;i++)
						{
							listVector[i][0]=fProcessor.getToken();
							listVector[i][1]=fProcessor.getNextToken();
							//go to next token end ommit "," (at end go on ";")
							fProcessor+=3;
						}

						if( !( (*accIOObject)->haveVectorField(fieldName ) ) ) //add new list to new token
						{	
							(*accIOObject)->addVectorFiled(fieldName,listVector);
						}
						else // fill existing list
						{
							(*accIOObject)->getVectorFieldAt(fieldName).insert((*accIOObject)->getVectorFieldAt(fieldName).end(),listVector.begin(),listVector.end());
						}

						break;
					}
				case Tensor:
					{
						vector< array<string,4> > listTensor(size);

						for(int i=0;i<size;i++)
						{
							listTensor[i][0]=fProcessor.getToken();
							listTensor[i][1]=fProcessor.getNextToken();
							fProcessor+=2; //go to next 2
							listTensor[i][2]=fProcessor.getToken();
							listTensor[i][3]=fProcessor.getNextToken();
							//go to next token end ommit "," (at end go on ";")
							fProcessor+=3;
						}

						if( !( (*accIOObject)->haveTensorField(fieldName ) ) ) //add new list to new token
						{	
							(*accIOObject)->addTensorFiled(fieldName,listTensor);
						}
						else // fill existing list
						{
							(*accIOObject)->getTensorFieldAt(fieldName).insert((*accIOObject)->getTensorFieldAt(fieldName).end(),listTensor.begin(),listTensor.end());
						}

						break;
					}

				case Table:
					{
						vector< vector<string> > listTable(size);
						vector<string> table;
						bool endTable;
						for(int i=0;i<size;i++)
						{
							endTable=false;

							while(!endTable)
							{
								table.push_back(fProcessor.getToken());
								fProcessor++;

								if(size!=1)
								{
									if(fProcessor.getToken()=="," || fProcessor.getToken()==")")
									{
										endTable=true;
										fProcessor++;//go to next digit, ommit ","
									}
								}
								else
								{
									if(fProcessor.getToken()==")")
									{
										endTable=true;
									}
								}

							}

							listTable[i]=table;
							table.clear();
						}

						if( !( (*accIOObject)->haveTableField(fieldName ) ) ) //add new list to new token
						{	
							(*accIOObject)->addTableFiled(fieldName,listTable);
						}
						else // fill existing list
						{
							(*accIOObject)->getTableFieldAt(fieldName).insert((*accIOObject)->getTableFieldAt(fieldName).end(),listTable.begin(),listTable.end());
						}

					break;
					}
				default:
					cout<<"Error! not acceptable list type"<<endl;
					break;
				}


				return true;
			}
			else
			{
				return false;
			}

		};

		////------- IOList ---------------------
		//	IOList::IOList
		//	(
		//	): IOToken("List")
		//	{
		//	};
		//
		//	bool IOList::act
		//	(
		//		FileProcessor& fProcessor,
		//		string** currField,
		//		string** currIOSubObject,
		//		IOSubObject** accIOObject
		//	)
		//	{
		//		if( token == fProcessor.getToken() )
		//		{
		//			return true;
		//		}
		//		else
		//		{
		//			return false;
		//		}
		//
		//	};


	}
}