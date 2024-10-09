///***************************************
//potentials used in this work
//***********************************

#include "libinclude.h"
#include "fundec.h"
#include "mathdec.h"
#include "const.h"
//Energies are obtained in Kj/mol



//Lennard Jones Potential
double lj_pot(double Dist, double *Par)
{
	double pot=0.000;

//	cout<<endl<<"Entering the Lennard Jones routine"<<endl;

//	cout<<Par[0]<<"\t"<<Par[1]<<endl;	

	if(Dist < 12)
		pot = 4*Par[0]*(pow(Par[1]/Dist,12)-pow(Par[1]/Dist,6));

//	cout<<pot<<endl;

	return pot;
}


//Morse Potential
double morse_pot(double Dist, double *Par)
{

//	cout<<endl<<"Entering the morse potential routine"<<endl;
	double pot=0.000;

//	cout<<Par[0]<<"\t"<<Par[1]<<"\t"<<Par[2]<<endl;	
	
	pot = Par[0]*(exp(-2*Par[1]*(Dist-Par[2]))-2*exp(-Par[1]*(Dist-Par[2])));

//	cout<<pot<<endl;

	return pot;

}



//Forces are obtained in Kj/mol/Angs
//remember while calculating forces one needs to apply minimum image convention 

void lj_force(int atom1, int atom2, double Dist, double **PosIons, double **ForceIons, double *Par, float **box)
{
	double Dx[3];
	double Dx1;
	double tmp;

	if(Dist < 12)
	{
		Dx[0] = PosIons[atom1][0] - PosIons[atom2][0];
		Dx[1] = PosIons[atom1][1] - PosIons[atom2][1];
		Dx[2] = PosIons[atom1][2] - PosIons[atom2][2];

		for(int i=0;i<3;i++)
		{

			Dx1 = Dx[i] - box[0][i]*ceil(Dx[0]/box[0][0]-0.5) - box[1][i]*ceil(Dx[1]/box[1][1]-0.5) - box[2][i]*ceil(Dx[2]/box[2][2]-0.5);
			tmp = 24*Par[0]*Dx1*(2*pow(Par[1]/Dist,12) - pow(Par[1]/Dist,6))/Dist/Dist;

			ForceIons[atom1][i] = ForceIons[atom1][i] + tmp;

			ForceIons[atom2][i] = ForceIons[atom2][i] - tmp;
		}

	}	

}

void morse_force(int atom1, int atom2, double Dist, double **PosIons, double **ForceIons, double *Par, float **box)
{
	double Dx[3];
	double Dx1;
	double tmp;

	Dx[0] = PosIons[atom1][0] - PosIons[atom2][0];
	Dx[1] = PosIons[atom1][1] - PosIons[atom2][1];
	Dx[2] = PosIons[atom1][2] - PosIons[atom2][2];


	for(int i=0;i<3;i++)
	{
		Dx1 = Dx[i] - box[0][i]*ceil(Dx[0]/box[0][0]-0.5) - box[1][i]*ceil(Dx[1]/box[1][1]-0.5) - box[2][i]*ceil(Dx[2]/box[2][2]-0.5);	
		tmp = Dx1*Par[0]*2*Par[1]*(exp(-2*Par[1]*(Dist-Par[2])) -  exp(-Par[1]*(Dist-Par[2])))/Dist;

		ForceIons[atom1][i] = ForceIons[atom1][i] + tmp;
		ForceIons[atom2][i] = ForceIons[atom2][i] - tmp;
	}

}



void fixatoms_forces(int natoms, double **ForceIons, int **fixatoms)
{
	for(int i = 0; i< natoms;i++)*********
	{
		for(int j=0;j<3;j++)
		{
			if(fixatoms[i][j] == 1)
				ForceIons[i][j] = 0.000;
		}
	}

}


