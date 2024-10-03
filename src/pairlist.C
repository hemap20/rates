//***********************************************
//function to generate pair list!! 
//***********************************************

#include "libinclude.h" 
#include "mathdec.h"

pair <int, int> compare(int n_atomtype,string *atomtype,  string atom1, string atom2);

void pairlist(double **PosIons, int natoms, float **boxcell, string *atomtype, int *natoms_type, int n_atomtype, string *type)
{

	//open the file to get cutoffs for bonds in angstrom 

	ifstream CutoffIn("cutoff.in",ios::in);
	if(!CutoffIn)
	{
		cerr << "File  CutoffIn could not be opened" <<endl;
		exit(1);
	}


	int ncutoff;    //number of cuttoff parameters
	string garbage;
	string *parameters;  //potential parameters
	float *cutoff;

	ncutoff = n_atomtype*(n_atomtype+1)/2;

	cutoff=new float [ncutoff];
	parameters=new string [ncutoff];

	for(int i=0;i<ncutoff;i++)
	{
		CutoffIn>>garbage>>garbage>>cutoff[i];
		getline(CutoffIn, parameters[i]);

	}

	int tmp = 0;
	CutoffIn.close();


	//generate pair list


	int num1, num2;
	double distatoms;

	FILE *pbondlist;
	FILE *pnonbondlist;
	FILE *pmolecules;

	pbondlist = fopen("bondlist.out", "w");
	pnonbondlist= fopen("nonbondlist.out", "w");
	pmolecules = fopen("molecule.out", "w");


	for(int i=0; i<natoms; i++)
	{
		for(int j=i+1;j<natoms; j++)
		{
			///calculate distance
			distatoms=dist(PosIons, i, j, boxcell);
			pair <int, int> comp = compare(n_atomtype,atomtype, type[i] , type[j]);
			num1 =comp.first;
			num2 = comp.second;

			//total number of pairs = nC2 + n


			tmp = (n_atomtype*(n_atomtype+1) - (n_atomtype-num1+1)*(n_atomtype-num1+2))/2 + (num2-num1+1);		

			//this will work only if the numbering is same in the input POSCAR file and cutoff.in file...be careful while making your inputs

			if(distatoms<=cutoff[tmp-1])
			{
				//	cout<<type[i]<<"-"<<type[j]<<"   bonded   "<<distatoms<<"  "<<cutoff[tmp-1]<<endl;

				fprintf(pbondlist, "%d - %d  %11.6f %s\n",i,j,distatoms,parameters[tmp-1].c_str());
				if(type[i].compare("C") ==0 || type[j].compare("C") ==0)
				{
					fprintf(pmolecules, "%d  -  %d  \n",i,j);

				}


			}
			else
			{
				if(type[i].compare("Cu") ==0 && type[j].compare("Cu") ==0)
				{
					if(distatoms<=6)
					 fprintf(pnonbondlist, "%d - %d  %11.6f %s\n",i,j,distatoms,parameters[tmp-1].c_str());

				}
				else
					fprintf(pnonbondlist, "%d - %d  %11.6f %s\n",i,j,distatoms,parameters[tmp-1].c_str());
			}

		}

	}

	cout<<"Bonded atoms with the cuttoff specified by the user are dumped in file bondlist.out\n";
	cout<<"Non-bonded atoms are dumped in the file nonbondlist.out\n\n";
	fclose(pbondlist);
	fclose(pnonbondlist);
	fclose(pmolecules);
	system("wc -l molecule.out > nmoles");
	system("cat nmoles molecule.out > molecules.out");
	system("rm nmoles molecule.out");
	delete [] cutoff;
	delete [] parameters;
}


// function compares atom1 and atom2 to atomtypes and return  a pair of two numbers for comparison
pair <int, int> compare(int n_atomtype, string *atomtype,  string atom1, string atom2 )
{
	int num1, num2;

	for(int i=1;i<=n_atomtype;i++)
	{
		if(atomtype[i-1]==atom1)
		{
			num1=i;

			break;
		}

	}		


	for(int i=1;i<=n_atomtype;i++)
	{
		if(atomtype[i-1]==atom2)
		{
			num2=i;

			break;
		}

	}


	return make_pair(num1, num2);
}



void readbondpar(int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int **pairs, int natoms)
{

	string garbage;
	ifstream InputIn("bondlist.out",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}
	for(int i=0;i<nbondatoms; i++)
	{

		InputIn>>batom1[i]>>garbage>>batom2[i]>>garbage;

		pairs[batom1[i]][batom2[i]] = -(i+1);
		pairs[batom2[i]][batom1[i]] = -(i+1);


		InputIn>>bondpot[i]>>garbage;

		if(bondpot[i].compare("lj") ==0)
		{
			InputIn>>bondpar[i][0]>>bondpar[i][1];
		}
		else if(bondpot[i].compare("ms") ==0)
		{
			InputIn>>bondpar[i][0]>>bondpar[i][1]>>bondpar[i][2];
		}
		else if(bondpot[i].compare("ha") ==0)
		{
			InputIn>>bondpar[i][0]>>bondpar[i][1];
		}
		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}


		getline(InputIn, garbage);
//		printf("%d - %d  %s %11.5f %11.5f %11.5f\n", batom1[i], batom2[i],bondpot[i].c_str(), bondpar[i][0], bondpar[i][1], bondpar[i][2]);  

	}


	InputIn.close();

}


void readnonbondpar(int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **pairs, int natoms)
{

	string garbage;
	ifstream InputIn("nonbondlist.out",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}
	for(int i=0;i<nnonbondatoms; i++)
	{

		InputIn>>nbatom1[i]>>garbage>>nbatom2[i]>>garbage;
		InputIn>>garbage>>nonbondpot[i];

		pairs[nbatom1[i]][nbatom2[i]] = (i+1);
		pairs[nbatom2[i]][nbatom1[i]] = (i+1);

		if(garbage.compare("lj") ==0)
		{
			InputIn>>garbage>>garbage;
		}
		else if(garbage.compare("ms") ==0)
		{
			InputIn>>garbage>>garbage>>garbage;
		}
		else if(garbage.compare("ha") ==0)
		{
			InputIn>>garbage>>garbage;
		}
		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}
	
		if(nonbondpot[i].compare("lj") ==0)
		{
			InputIn>>nonbondpar[i][0]>>nonbondpar[i][1];
		}
		else if(nonbondpot[i].compare("ms") ==0)
		{
			InputIn>>nonbondpar[i][0]>>nonbondpar[i][1]>>nonbondpar[i][2];
		}
		else if(nonbondpot[i].compare("ha") ==0)
		{
			InputIn>>nonbondpar[i][0]>>nonbondpar[i][1];
		}
		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}



//		printf("%d - %d  %s %11.5f %11.5f %11.5f\n", nbatom1[i], nbatom2[i],nonbondpot[i].c_str(), nonbondpar[i][0], nonbondpar[i][1], nonbondpar[i][2]);  

	}


	InputIn.close();

}

