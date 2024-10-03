//***************************************
// Calculate energy and forces 
//***************************************

#include "libinclude.h"
#include "mathdec.h"
#include "potdec.h"
#include "const.h"


// The following function calculates energy for a single atom selected 
double CalEner1(double **PosIons, int natoms, float **boxcell, int *batom1, int *batom2, string *bondpot, int *nbatom1, int *nbatom2, string *nonbondpot, double **bondpar, double **nonbondpar, int **pairs, int randatom, int intdoub, int randatom2)
{

	double Pot = 0;
	double Dist;
	double *Par;
	int tmp, natom;
	string pot;

	Par=new double [10]; //upto 10 empirical parameters

//	cout<<"randatom="<<randatom<<endl;

	double (*potfunc)(double, double *);

	for(int i=0;i<natoms;i++)
	{
		if(intdoub==1 && i==randatom2)
		{
			continue;
		
		}
		//calculating distance between atoms
		if(pairs[randatom][i] != 0)
		{
			tmp = pairs[randatom][i];
//			cout<<tmp<<endl;

			if(tmp<0)
			{
				natom = -tmp -1;
				pot= bondpot[natom];

//				cout<<endl<<natom<<"\t"<<pot<<"\t"<<batom1[natom]<<"\t"<<batom2[natom]<<endl;

				Dist=dist(PosIons, batom1[natom], batom2[natom], boxcell);

//				cout<<"Distance between atoms: "<<Dist;
				

				if(pot.compare("lj") ==0)
				{
					potfunc = lj_pot;
					Par[0] = bondpar[natom][0];
					Par[1] = bondpar[natom][1];
				}
				else if(pot.compare("ms") ==0)
				{
					potfunc = morse_pot;
					Par[0] = bondpar[natom][0];
					Par[1] = bondpar[natom][1];
					Par[2] = bondpar[natom][2];
				}

				else
				{
					cerr<<"Check your input, function not available in this version of the Program";
				}

			}
			else
			{
				natom = tmp -1;
				pot= nonbondpot[natom];

//				cout<<endl<<natom<<"\t"<<pot<<"\t"<<nbatom1[natom]<<"\t"<<nbatom2[natom]<<endl;

				Dist=dist(PosIons, nbatom1[natom], nbatom2[natom], boxcell);

//				cout<<"Distance between atoms: "<<Dist;
				
				if(pot.compare("lj") ==0)
				{
					potfunc = lj_pot;
					Par[0] = nonbondpar[natom][0];
					Par[1] = nonbondpar[natom][1];
				}
				else if(pot.compare("ms") ==0)
				{
					potfunc = morse_pot;
					Par[0] = nonbondpar[natom][0];
					Par[1] = nonbondpar[natom][1];
					Par[2] = nonbondpar[natom][2];
				}

				else
				{
					cerr<<"Check your input, function not available in this version of the Program";
				}


			}

//			cout<<"Potential at every step: "<< Pot<<endl;

			Pot = Pot + potfunc(Dist, Par);

		}
	}


	delete [] Par;	
	return Pot;

}



double CalEner(double **PosIons, int natoms, float **boxcell, int nbatoms, int *atom1, int *atom2, string *npot, double **npar)
{

	double Pot = 0;
	double Dist;
	double *Par;

	Par=new double [10]; //upto 10 empirical parameters


	double (*potfunc)(double, double *);

	for(int i=0;i<nbatoms;i++)
	{
		//calculating distance between atoms
		Dist=dist(PosIons, atom1[i], atom2[i], boxcell);

		if(npot[i].compare("lj") ==0)
		{
			potfunc = lj_pot;
			Par[0] = npar[i][0];
			Par[1] = npar[i][1];
		}
		else if(npot[i].compare("ms") ==0)
		{
			potfunc = morse_pot;
			Par[0] = npar[i][0];
			Par[1] = npar[i][1];
			Par[2] = npar[i][2];
		}

		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}


		Pot = Pot + potfunc(Dist, Par);
	}


	delete [] Par;	
	return Pot;

}

void CalFor(double **PosIons, double **ForceIons,  int natoms, float **boxcell, int nbatoms, int *atom1, int *atom2, string *npot, double **npar)
{

	double Pot = 0;
	double Dist;
	double *Par;

	Par=new double [10]; //upto 10 empirical parameters


	void (*forcefunc)(int, int, double, double **, double **, double *, float **);	


	for(int i=0;i<nbatoms;i++)
	{
		//calculating distance between atoms
		Dist=dist(PosIons, atom1[i], atom2[i], boxcell);

		if(npot[i].compare("lj") ==0)
		{

			forcefunc = lj_force;
			Par[0] = npar[i][0];
			Par[1] = npar[i][1];
		}
		else if(npot[i].compare("ms") ==0)
		{
			forcefunc = morse_force;
			Par[0] = npar[i][0];
			Par[1] = npar[i][1];
			Par[2] = npar[i][2];
		}

		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}

		forcefunc(atom1[i], atom2[i], Dist, PosIons, ForceIons, Par, boxcell);

	}

	delete [] Par;

}


double CalEnerFor(double **PosIons, int natoms, float **boxcell, double **ForceIons, string filename, bool bonded, bool forces)
{
	string garbage, garbage1;
	char bondpot[2], nonbondpot[2];
	double Pot = 0;
	double Dist;
	int atom1, atom2;
	int npar;
	double *Par;

	Par=new double [10]; //upto 10 empirical parameters




	double (*potfunc)(double, double *);
	void (*forcefunc)(int, int, double, double **, double **, double *, float **);	



	ifstream FileIn(filename.c_str(),ios::in);
	if(!FileIn)
	{
		cerr << "File BondIn could not be opened" <<endl;
		exit(1);
	}

	// Get energy and forces from the bonded atoms	

	while(!FileIn.eof())
	{
		FileIn>>atom1>>garbage>>atom2>>Dist;

		if(FileIn.eof())  break;

		//calculating distance between atoms
		Dist=dist(PosIons, atom1, atom2, boxcell);

		FileIn>>bondpot>>nonbondpot;

		if(strcmp(bondpot,"lj") ==0)
		{
			potfunc = lj_pot;
			forcefunc = lj_force;
			FileIn>>Par[0]>>Par[1];
		}
		else if(strcmp(bondpot,"ms") ==0)
		{
			potfunc = morse_pot;
			forcefunc = morse_force;
			FileIn>>Par[0]>>Par[1]>>Par[2];
		}

		else
		{
			cerr<<"Check your input, function not available in this version of the Program";
		}


		if(bonded)
		{
			getline(FileIn, garbage);
			Pot = Pot + potfunc(Dist, Par);

			if(forces)
				forcefunc(atom1, atom2, Dist, PosIons, ForceIons, Par, boxcell);
		}
		else
		{
			if(strcmp(nonbondpot,"lj") ==0)
			{
				potfunc = lj_pot;
				forcefunc = lj_force;
				FileIn>>Par[0]>>Par[1];
			}
			else if(strcmp(nonbondpot,"ms") ==0)
			{
				potfunc = morse_pot;
				forcefunc = morse_force;
				FileIn>>Par[0]>>Par[1]>>Par[2];
			}
			else
			{
				cerr<<"Check you input, function not available in this version of the Program";
				exit (EXIT_FAILURE);
			}


			Pot = Pot + potfunc(Dist, Par); 


			if(forces)
				forcefunc(atom1, atom2, Dist, PosIons, ForceIons, Par, boxcell);
		}

	}

	cout <<"Energy from this routine: "<<Pot;
	delete [] Par;
	return Pot;
}


double CalKinEner(double **Vel, int natoms, float *mass)
{

	double KEner=0;
	for(int i=0;i<natoms;i++)
	{		
		for(int j=0;j<3;j++)
		{
			KEner = KEner + mass[i]*Vel[i][j]*Vel[i][j]/2;
		}

	}

	//KEner obtained is in g-angs2/mol-fs2, converting it into kJ/mol
	return KEner*1e4;


}

float CalTemp(double KEner, int natoms)
{
	float Temp;


	Temp = 2*KEner/(3*natoms-1)*1e6/KB/AVG;

	return Temp;


}
