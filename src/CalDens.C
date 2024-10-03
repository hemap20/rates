//***************************************************************
// Code for computing densities 
//***************************************************************

#include "libinclude.h"
#include "mathdec.h"
#include "potdec.h"
#include "fundec.h"
#include "const.h"

void CalDens(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, double lim1, double dx, double *density)
{

	int index;
	int j;

		
	
	for(j=natoms_type[0]; j<natoms; j++)
	{
		index = floor((PosIons[j][0]-lim1)/dx);
		density[index] = density[index] + 1;
	}


}
