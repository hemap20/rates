//***********************************************
// function to for doing simulated annealing 
// 
// //***********************************************

#include "libinclude.h" 
#include "fundec.h"
#include "mathdec.h"
#include "const.h"
#include "potdec.h"



void simuAnn(double **PosIons, int natoms, float **boxcell, double **vel, double **ForceIons, float *mass, float dt, float TempB, float TempE, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms, int rampsize, int rampstep)
{

	int sasteps = ((TempE-TempB)/rampsize);
	float Temp;
	for(int i=0;i<abs(sasteps);i++)
	{
		Temp  = TempB + i*sasteps/abs(sasteps)*rampsize;

		for(int j=0;j<rampstep;j++)
		{

			res_vel(vel, natoms, Temp, mass);
			zeros(natoms, 3 ,ForceIons);
			CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
			CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);
			fixatoms_forces(natoms, ForceIons, fixatoms);

			md(PosIons, natoms, boxcell, vel, ForceIons, mass, dt, Temp, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar,fixatoms);


		}

		cout<<"*****Wrapping Coordinates"<<endl;
		wrap_coor(PosIons, natoms, boxcell);

	}
}
