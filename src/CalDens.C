//***************************************************************
// Code for computing densities 
//***************************************************************

#include "libinclude.h"
#include "mathdec.h"
#include "potdec.h"
#include "fundec.h"
#include "const.h"

void CalDens(double **PosIons, int natoms, float **boxcell, int n_atomtype, int *natoms_type, double lim1, double dx, double *density, int density_size)
{
    int index, j;

    // Each thread gets a private copy of a local density array
    #pragma omp parallel
    {
        // Create a local density array initialized to 0
        vector<double> local_density(density_size, 0);

        #pragma omp for private(index)
        for (j = natoms_type[0]; j < natoms; j++)
        {
            index = floor((PosIons[j][0] - lim1) / dx);
            local_density[index] += 1;
        }

        // Combine the local densities into the shared array
        #pragma omp critical
        {
            for (int i = 0; i < density_size; i++)
            {
                density[i] += local_density[i];
            }
        }
    }
}
