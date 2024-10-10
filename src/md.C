//***************************************
// Molecular Dynamics step 
//***************************************

#include "libinclude.h"
#include "fundec.h"
#include "mathdec.h"
#include "const.h"
#include "potdec.h"


///velocity-verlet algorithm
// Positions are in Angs, Forces are in Kj/mol/angs, velocities are ang/fs, time in fs 

void md(double **PosIons, int natoms, float **boxcell, double **vel, double **ForceIons, float *mass, float dt, float Temp, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms )
{

	//convert forces in kJ/mol-angs to gm-angs^2/fs^2-mol
		
	double forc_conv = 1e-4;


	for(int i=0;i<natoms;i++)
	{
		for(int j=0;j<3;j++)
		{

				vel[i][j] = vel[i][j]  + dt*ForceIons[i][j]*forc_conv/2/mass[i];

				PosIons[i][j] = PosIons[i][j] + dt*vel[i][j];  			
		}
	
	}


		zeros(natoms, 3 ,ForceIons);
		CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
		CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);
		fixatoms_forces(natoms, ForceIons, fixatoms);



	for(int i=0;i<natoms;i++)
	{
		for(int j=0;j<3;j++)
		{
				vel[i][j] = vel[i][j]  + dt*ForceIons[i][j]*forc_conv/2/mass[i];
			
		}
	
	}


}



