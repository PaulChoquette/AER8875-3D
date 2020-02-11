
#include <iostream>
#include <metis.h>

int nElements=26;
int nNodes=36;

int eprt_dummy[37]={0,6,10,14,18,22,28,32,36,40,44,48,52,56,59,62,66,70,74,78,82,86,90,94,98,102,106,110,114,118,122,126,130,134,138,142,146};
int eind_dummy[150]={1,2,3,12,27,34,11,12,33,34,10,11,32,33,9,10,31,32,8,9,30,21,5,6,7,8,29,30,4,5,28,29,3,4,15,16,4,5,16,17,5,6,17,18,6,18,19,1
,13,26,1,12,25,26,11,12,24,25,10,11,23,24,9,10,22,23,8,9,21,22,7,8,20,21,6,7,19,20,19,20,41,42,20,21,42,43,21,22,43,44,22,23,44,45,23,
24,45,46,24,25,46,47,25,26,47,48,13,26,35,48,13,14,36,37,15,36,37,15,16,37,38,16,17,38,39,17,18,39,40,18,19,40,41};


int* epart;
int* npart;
int* eprt;
int* eind;


int main()
{

	int nPart=2;

	int objval;
	int ncommon=2;

	epart = new int [nElements];
	npart = new int [nNodes];
	eprt = new int [37];
	eind = new int [150];

	for (int i =0; i<37; i++)
	{
		eprt[i]=eprt_dummy[i];
	}
	for (int i=0; i<150; i++)
	{
		eind[i]=eind_dummy[i];
	}

	int success = METIS_PartMeshDual(&nElements, &nNodes, &eprt[0], &eind[0], NULL, NULL,
																	 &ncommon, &nPart, NULL, NULL, &objval,
																	 &epart[0], &npart[0]);

		for (int i=0;i<nElements; i++)
		{
			std::cout<<epart[i]<<std::endl ;
		}
		std:: cout<< std::endl;
		for (int i=0; i<nNodes; i++)
		{
			std::cout<<npart[i]<<std::endl;
		}
};
