//----------------------------------------------------------------
// V. Agarwal 
// Starting date: 21.06.19
// Last Modified:
// code for computing rates in a liquid
//-----------------------------------------------------------------

#include "libinclude.h" 
#include "potdec.h"
#include "mathdec.h"
#include "fundec.h"
#include "const.h"

void vec_zero(double *vec, int len);
int main(int argc, char **argv)
{
	cout<<"******************************************************************"<<endl;
	cout<<"******************************************************************"<<endl;
	cout<<"      Code for adsorption/desorption near a sold surface          "<<endl;
	cout<<"                    Written by Vishal Agarwal                     "<<endl;
	cout<<"******************************************************************"<<endl;
	cout<<"******************************************************************"<<endl<<endl<<endl;

	//Start Clock;
	time_t tim;
	time(&tim); //pass variable tim to time function 
	cout << "Simulations started on: "<< ctime(&tim)<<endl; // this translates what was returned from time() into a readable format

	struct timeb starttime, endtime;
	ftime(&starttime);
	/////////define variables

	string filename;
	string posfile;
	string garbage, garbage1;
	char garbage2[200];
	int tmp;
	int i, j, MC, MD;
	int  nMD, md_print;
	float Temp;
	int  MDtime, MDeq;
	float dt;
	int genpairlist;
	int randnum;
	int hyb, NMCMD, rdvel;
	float tmp2;

	///////////////////////////////////////////////
	//////initialize variables////////////////////
	filename = "input.in";
	ifstream InputIn(filename.c_str(),ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	getline(InputIn, garbage);

	cout<<"Project Name: "<<garbage<<endl;

	InputIn>>garbage>>garbage;
	InputIn>>posfile;


	InputIn>>garbage>>garbage;
	InputIn>>genpairlist;


	InputIn>>garbage>>garbage;
	InputIn>>Temp;

	InputIn>>garbage>>garbage;
	InputIn>>hyb;

	InputIn>>garbage>>garbage;
	InputIn>>NMCMD;

	InputIn>>garbage>>garbage;
	InputIn>>rdvel;

	bool comdens;

	InputIn>>garbage>>garbage;
	InputIn>>comdens;

	cout<<"Compute Density:  "<<comdens;

	double lim1, lim2, dx;

	InputIn>>garbage>>garbage;
	InputIn>>lim1>>lim2>>dx;

	cout<<"Range of x-values considered for computing densities: "<<lim1<<"\t"<<lim2<<endl;

	//Parameters for simulated annealing		

	bool simann;

	InputIn>>garbage>>garbage;
	InputIn>>simann;

	float TempB, TempE;

	InputIn>>garbage>>garbage;
	InputIn>>TempB>>TempE;

	int rampsize, rampstep;

	InputIn>>garbage>>garbage;
	InputIn>>rampsize>>rampstep;

	int ncycle;

	InputIn>>garbage>>garbage;
	InputIn>>ncycle;
	

	if(simann==1)
	{	
		cout<<"Perform Simulated Annealing"<<endl;

		cout<<"Start Temperature: "<<TempB<<" K"<<endl;
		cout<<"End Temperature: "<<TempE<<" K"<<endl;
		cout<<"Delta T:"<<rampsize<<" K"<<endl;
		cout<<"Increase temperature after: "<<rampstep<<" MD steps"<<endl;
		cout<<"ncycles: "<<ncycle<<endl;
	}

	bool minimize;

	InputIn>>garbage>>garbage;
        InputIn>>minimize;                      //whether to minimize or not

	if(minimize)
		cout<<"****Perform minimization****"<<endl;

	InputIn.close();

	ifstream MDIn("md.in",ios::in);
	if(!MDIn)
	{
		cerr << "File MDIn could not be opened" <<endl;
		exit(1);
	}

	MDIn>>garbage>>garbage;
	MDIn>>MDtime;

	MDIn>>garbage>>garbage;
	MDIn>>dt;                                    

	MDIn>>garbage>>garbage;
	MDIn>>md_print;                             

	MDIn>>garbage>>garbage;
	MDIn>>MDeq;                             

	bool trj;

	MDIn>>garbage>>garbage;
	MDIn>>trj;                             

	MDIn.close();

	int nMDeq = MDeq/dt;

	nMD = MDtime/dt;

	////////////////////////////////////////////////////////////

	cout<<"Initial positions read from the file: "<<posfile<<endl;

	//////////////////////get positions from the POSCAR file//////////////////////////////////////////////////////////

	ifstream PosIn(posfile.c_str(),ios::in);
	if(!PosIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	cout<<"Temperature of the system:\t\t"<<Temp<<" K"<<endl;

	getline(PosIn, garbage);
	istringstream StrStream(garbage);

	int n_atomtype=0;

	while(StrStream)
	{
		getline(StrStream, garbage1, ' ');
		if(garbage1.compare("") != 0)
			n_atomtype = n_atomtype + 1;
	}

	string *atomtype;
	atomtype=new string [n_atomtype];
	istringstream strstream(garbage);
	tmp = 0;

	while(strstream)
	{
		getline(strstream, garbage1, ' ');
		if(garbage1.compare("") != 0)
		{
			atomtype[tmp]=garbage1;
			tmp = tmp + 1;
		}
	}

	int *natoms_type, natoms;  //number of atoms for each type
	natoms_type=new int [n_atomtype];
	natoms=0;

	//get box cell vectors
	float **boxcell;
	boxcell=new float * [3];

	for(i = 0; i<3;i++)
	{
		boxcell[i]=new float [3];
	}

	getline(PosIn, garbage);

	for(i = 0; i<3;i++)
	{
		for(j=0;j<3;j++)
		{
			PosIn>>boxcell[i][j];
		}
	}
	
	getline(PosIn, garbage);
	getline(PosIn, garbage);

	for(i=0;i<n_atomtype;i++)
	{
		PosIn>>natoms_type[i];
		natoms=natoms+natoms_type[i];
	}

	getline(PosIn, garbage);
	getline(PosIn, garbage);

	double **PosIons, **ForceIons, **vel;
	int **fixatoms;
	float *mass;

	mass=new float [natoms];
	PosIons=new double * [natoms];
	ForceIons=new double *[natoms];
	vel=new double *[natoms];
	fixatoms=new int *[natoms];

	for(i=0;i<natoms;i++)
	{
		PosIons[i]=new double [3];
		ForceIons[i]=new double [3];
		vel[i]=new double [3];		
		fixatoms[i]=new int [3];
	}

	//getting positions of each atom
	for(i=0;i<natoms;i++)
	{
		PosIn>>PosIons[i][0]>>PosIons[i][1]>>PosIons[i][2];
	}


	PosIn.close();

	//////////////////////////////////////////////////////////////////////////////////////////////////////////////

	// Generate pair list
	// initialize masses

	string *type;
	type=new string [natoms];

	tmp = 0;

	//get the vector of symbols for atoms
	for(i=0;i<n_atomtype;i++)
	{ 
		for(int j=0;j<natoms_type[i];j++)
		{
			type[tmp] = atomtype[i];
			tmp=tmp+1;
		}		
	}

	//get masses
	ini_mass(mass, natoms, type);

	//read which atoms to fix /////////////////////////////////////

	InputIn.open("fixatoms.in");

	cout<<"***Reading fixatoms file***"<<endl;

	for(i = 0; i< natoms;i++)
	{
		cout<<type[i]<<"\t";
		for(int j=0;j<3;j++)
		{
			InputIn>>garbage;

			if(garbage.compare("T") !=0)
				fixatoms[i][j] = 1;
			else
				fixatoms[i][j] = 0;

			cout<<fixatoms[i][j]<<"\t";
		}
		cout<<endl;
	}

	InputIn.close();


	///////////////////////////////////////////////////////////////////////////////////////


	/////////////////////////////////////////////////////////////////////////////////////////
	if(genpairlist==1)
	{
		cout<<"\n***********Generating pair list**************\n";
		pairlist(PosIons, natoms, boxcell, atomtype,natoms_type, n_atomtype,type);
	}
	////////////////////////////////////////////////////////////////////////////////////////


	//store bonded and nonbonded pairs in a file and store their parameters///////////////////

	int nbondatoms, nnonbondatoms;

	system("wc -l bondlist.out > tmp.txt");
	InputIn.open("tmp.txt");
	InputIn>>nbondatoms;
	InputIn.close();

	cout<<"No. of bonded interactions: "<<nbondatoms<<endl;

	system("wc -l nonbondlist.out > tmp.txt ");
	InputIn.open("tmp.txt");
	InputIn>>nnonbondatoms;
	InputIn.close();

	system("rm tmp.txt ");
	cout<<"No. of nonbonded interactions: "<<nnonbondatoms<<endl;

	int *batom1, *batom2, *nbatom1, *nbatom2, **pairs;
	string *bondpot, *nonbondpot;
	double **bondpar, **nonbondpar;

	batom1=new int [nbondatoms];
	batom2=new int [nbondatoms];
	nbatom1=new int [nnonbondatoms];
	nbatom2=new int [nnonbondatoms];
	pairs=new int *  [natoms];

	bondpot=new string [nbondatoms];
	nonbondpot=new string [nnonbondatoms];

	bondpar=new double * [nbondatoms];
	nonbondpar=new double * [nnonbondatoms];

	for(i=0;i<nbondatoms;i++)
	{
		bondpar[i]=new double [3];
	}

	for(i=0;i<nnonbondatoms;i++)
	{
		nonbondpar[i]=new double [3];
	}

	for(i=0;i<natoms;i++)
	{
		pairs[i]=new int [natoms];
	}


	//reading pairs of bonded and nonbonded atoms and storing these parameters in pointers
	zeros_int(natoms, natoms,pairs);
	readbondpar(nbondatoms, batom1,batom2, bondpot, bondpar, pairs, natoms);
	readnonbondpar(nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar, pairs, natoms);

	//print the pairs list
	//////////////////////////////////////////////////////////////////////////////////////////////////////
	//	print_mat_int(pairs, natoms, natoms, "pairs");
	// Get energy and forces///////////////////////////////////////////////////////////////////////////////

	double Pot = 0;
	Pot = CalEner(PosIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
	Pot = Pot + CalEner(PosIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);

	zeros(natoms, 3 ,ForceIons);
	CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
	CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);

	fixatoms_forces(natoms, ForceIons, fixatoms);

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	//Code for density calculations

	double *density;
	int len=(lim2-lim1)/dx;
	density=new double [len];
	vec_zero(density, len);
	int index;
	int nconf = 0;

	//print Coordinates, Forces//////////////////////////////////////////////////////////////////
	printCoor(PosIons, natoms, type);
	printFor(ForceIons,natoms,type);

	cout<<"\nTotal Potential Energy of the System:  "<<Pot/KJEV<<"\t\t eV"<<endl;


	//minimization code
	if(minimize)
        {       
        	minimization(PosIons, natoms, ForceIons, boxcell, Pot,type, fixatoms);
            print_coor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, 0,'w', "CONTCAR");

			Pot = CalEner(PosIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
			Pot = Pot + CalEner(PosIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);

			zeros(natoms, 3 ,ForceIons);
			CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
			CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);

			fixatoms_forces(natoms, ForceIons, fixatoms);

            printCoor(PosIons, natoms, type);
            printFor(ForceIons,natoms,type);
			print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");

            cout<<"\nTotal Potential Energy of the System:  "<<Pot/KJEV<<"\t\t eV"<<endl;
        }       

	//////////////////////////////////////////////////////////////////////////////////////////////////////////
	cout<<"******Intitializing Velocities******"<<endl<<endl;

	if(rdvel ==1)
	{	
		readvel(vel, natoms);
	}
	else
	{
		randnum =  rand() % 100000;
		ini_vel(vel,  natoms, Temp, mass, randnum, fixatoms);
	}

	if(simann==1)
	{
		for(j=0;j<ncycle;j++)
		{
			simuAnn(PosIons, natoms, boxcell, vel, ForceIons, mass, dt, TempB, TempE, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms,nbatom1, nbatom2, nonbondpot, nonbondpar, fixatoms, rampsize,rampstep);
			print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR.SA");
			tmp2 = TempB;
			TempB=TempE;
			TempE=tmp2;
		}
	}

	/////////////////
	if(hyb==1)
	{
		for(j=0;j<NMCMD;j++)
		{
			cout<<"*****hybrid MCMD step No. "<<j<<endl;
			if(j==0)
				print_coor(PosIons, natoms, boxcell, n_atomtype, natoms_type, atomtype, 0, j,'w', "XDATCAR");

			montecarlo(PosIons, natoms, boxcell, mass, type, n_atomtype, natoms_type, atomtype, Temp, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar, fixatoms, pairs);

			//////////////////////////////////////////////////////////////////////////////////////////////////////////

			cout<<"*****Entering the Molecular dynamics code*******"<<endl<<endl;
			cout<<"Starting Equilibration runs for:   "<<nMDeq<<endl<<endl;

			#pragma omp parallel
			{
				for (i = 0; i < nMDeq; i++)  
				{
					// Perform the velocity reset every 50 steps (handled outside the parallel loop)
					if (i % 50 == 0 && i != 0)
					{
						#pragma omp single // Ensure only one thread resets velocities
						{
							res_vel(vel, natoms, Temp, mass);
						}
					}

					// Reset forces for the current iteration
					zeros(natoms, 3, ForceIons);

					// Calculate forces
					CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
					CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);
					fixatoms_forces(natoms, ForceIons, fixatoms);

					// Perform molecular dynamics step
					md(PosIons, natoms, boxcell, vel, ForceIons, mass, dt, Temp, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar, fixatoms);

					// Print status every 10000 iterations (handled outside the parallel region)
					if (i % 10000 == 0)
					{
						#pragma omp critical // Ensure only one thread prints at a time
						{
							cout << "MD Eq step no.: " << i << endl;
							printVel(vel, natoms);
							print_carcoor(PosIons, natoms, boxcell, n_atomtype, natoms_type, atomtype, 0, i, 'w', "CONTCAR");
						}
					}

					// Wrap coordinates every 1000 iterations (handled outside the parallel region)
					if (i % 1000 == 0)
					{
						#pragma omp critical // Ensure only one thread prints at a time
						{
							cout << "MD Eq. Step No:  " << i << endl;
							cout << "*****Wrapping Coordinates" << endl;
							wrap_coor(PosIons, natoms, boxcell);
						}
					}
				}
			}

			cout<<"*****Wrapping Coordinates"<<endl;
			wrap_coor(PosIons, natoms, boxcell);

			print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR.eq");
			cout<<"****************Finished Equilibration runs"<<endl<<endl;

			#pragma omp parallel
			{
				// Declare an index variable for the thread
				int n_threads = omp_get_num_threads();
				int thread_id = omp_get_thread_num();

				// Parallelize the loop
				#pragma omp for
				for (int i = 0; i < nMD; i++)
				{
					if (i % md_print == 0)
					{
						// Use critical section to ensure only one thread prints at a time
						#pragma omp critical
						{
							cout << "MD step no.: " << i << endl;
							wrap_coor(PosIons, natoms, boxcell);
							if (trj)
								print_coor(PosIons, natoms, boxcell, n_atomtype, natoms_type, atomtype, 1, i, 'a', "XDATCAR");
							
							print_carcoor(PosIons, natoms, boxcell, n_atomtype, natoms_type, atomtype, 0, i, 'w', "CONTCAR");
							printMD(PosIons, natoms, boxcell, vel, i, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar, fixatoms, pairs, mass);
							printVel(vel, natoms);
							
							if (comdens)
							{
								nconf = nconf + 1;
								CalDens(PosIons, natoms, boxcell, n_atomtype, natoms_type, lim1, dx, density, len);
							}
						}
					}

					// Reset velocities every 10 iterations
					if (i % 10 == 0 && i != 0)
					{
						#pragma omp single // Ensure that only one thread resets velocities
						{
							res_vel(vel, natoms, Temp, mass);
						}
					}

					// Reset forces for the current iteration
					zeros(natoms, 3, ForceIons);

					// Calculate forces
					CalFor(PosIons, ForceIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
					CalFor(PosIons, ForceIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);
					fixatoms_forces(natoms, ForceIons, fixatoms);

					// Perform molecular dynamics step
					md(PosIons, natoms, boxcell, vel, ForceIons, mass, dt, Temp, nbondatoms, batom1, batom2, bondpot, bondpar, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar, fixatoms);
				}
			}
			print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");
		}	
	}


	if(comdens)
	{
		double vol = boxcell[1][1]*boxcell[2][2]*dx*AVG/1e24;      ///cm3/mol

		cout<<"**************Printing Densities"<<endl;
		for(i=0; i<len; i++)
		{
			cout<<(lim1+i*dx+dx/2)<<"\t"<<density[i]<<"\t"<<density[i]/vol/nconf<<endl; //calculates density in mol/cm3 
		}
	}


	// delete dynamic variables 

	for(i=0;i<3;i++)
	{
		delete [] boxcell[i]; 
	}

	for(i=0;i<natoms;i++)
	{
		delete [] PosIons[i];
		delete [] ForceIons[i];
		delete [] vel[i];
		delete [] fixatoms[i];
	}


	delete [] PosIons;
	delete [] boxcell;
	delete [] atomtype;
	delete [] natoms_type;
	delete [] ForceIons;
	delete [] vel;
	delete [] fixatoms;

	delete [] type;
	delete [] pairs;

	delete [] batom1;
	delete [] batom2;
	delete [] nbatom1;
	delete [] nbatom2;

	delete [] bondpot;
	delete [] nonbondpot; 
	delete [] density;


	for(i=0;i<nbondatoms;i++)
	{
		delete [] bondpar[i];
	}

	for(i=0;i<nnonbondatoms;i++)
	{

		delete [] nonbondpar[i];
	}


	delete [] bondpar;
	delete [] nonbondpar;







	///End time and output time
	ftime(&endtime);

	cout<<endl<<"Execution Time:  "<<duration(starttime,endtime)<<" seconds"<<endl;


	return 0;
}//end of main loop


double duration(struct timeb start, struct timeb end)
{ return(end.time+(end.millitm/1000.0)) - (start.time+(start.millitm/1000.0));
}
