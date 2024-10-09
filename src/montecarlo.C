//***************************************************************
// Monte Carlo Algorithm for getting ratio of partition functions
//***************************************************************

#include "libinclude.h"
#include "mathdec.h"
#include "potdec.h"
#include "fundec.h"
#include "const.h"
#include "dSFMT.h"



void print_MCdetail(int MCstep, string move, int atom1, int atom2, double dx, double dy, double dz, double Pot, double Vbias, double delPot, int accp_rej);

void montecarlo(double **PosIons, int natoms, float **boxcell, float *mass, string *type, int n_atomtype, int *natoms_type, string *atomtype,float Temp, int nbondatoms, int *batom1, int *batom2, string *bondpot, double **bondpar, int nnonbondatoms, int *nbatom1, int *nbatom2, string *nonbondpot, double **nonbondpar, int **fixatoms, int **pairs)
{

	cout<<"*****Entering the Monte Carlo code*******"<<endl<<endl;

	int seed;		//seed for random number generator
	dsfmt_t dsfmt;
	seed = time(NULL);
	dsfmt_init_gen_rand(&dsfmt, seed);
	int n_accp, n_rej, n_swap_accp, n_swap_rej;  //no. of moves which are accepted and rejected
	double delPot;    //change in the potential
	double prob;      //boltzman probability
	double Pot;       //Potential Energy

	float dq[3];                     //move in the three directions
	float randnum;
	float delx, delxMol, probdelx;   //delx is maximum atomic moves, delxMol is maximum molecule move, probdelx is the probability of atomic move
	char garbage[100];
	int printmc, niter;        //at what sweep to print MC output and no. of iterations
	int atom1, atom2;          //atom1 and atom2
	int fixmolecule;           //whether to fix molecule or not
	int fixatom1=1 , fixatom2 = 1;	 //which atom is fixed  

	int nrun;

	int tmp;	 
	int MCiter, MCeq;          //MC iteration and equilibration steps 

	int nseed;		//run random number nseed times before drawing numbers for monte carlo simulaitons

	int accp_rej=0;	

	ifstream InputIn("mc.in",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	InputIn>>MCiter;

	InputIn>>garbage>>garbage;
	InputIn>>MCeq;


	InputIn>>garbage>>garbage;
	InputIn>>delx;                                    // maximum movement of an atom

	InputIn>>garbage>>garbage;
	InputIn>>probdelx;                                    // maximum movement of an atom


	InputIn>>garbage>>garbage;
	InputIn>>printmc;


	InputIn>>garbage>>garbage;
	InputIn>>nseed;			   		    



	InputIn.close();


	/////random number generator
	for(int i=0;i<nseed; i++)
		dsfmt_genrand_close_open(&dsfmt);



	//calculate energy
	Pot = CalEner(PosIons, natoms, boxcell, nbondatoms, batom1, batom2, bondpot, bondpar);
	Pot = Pot + CalEner(PosIons, natoms, boxcell, nnonbondatoms, nbatom1, nbatom2, nonbondpot, nonbondpar);

	int randatom; 

	n_accp = 0;
	n_rej = 0;
	n_swap_accp = 0;
	n_swap_rej = 0;

	nrun = natoms;

	for(int r=0;r<2;r++)
	{
		if(r==0)
			niter=MCeq;
		else
		{
			niter=MCiter;
			n_accp = 0;
			n_rej = 0;
		}




		for(unsigned int i=0;i<niter*nrun;i++)************************
		{

			//pick a random number and decide whether to swap the atom or not

			randnum=dsfmt_genrand_close_open(&dsfmt);

			if(randnum<=probdelx)
			{
				//Monte Carlo move to just move the atoms
				//pick a random atom

				delPot = 0;
				randatom = dsfmt_genrand_close_open(&dsfmt)*natoms;

				//		cout<<"Random atom chosen is: "<<randatom<<endl;

				if(fixatoms[randatom][0] != 0 && fixatoms[randatom][1] != 0 && fixatoms[randatom][2] != 0)
				{
					//			cout<<"Atom "<<randatom<<"  is fixed"<<endl;
					continue;
				}


				delPot = CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, randatom, 0, 3);


				for(int j=0;j<3;j++)
				{
					//randomly move atom
					if(fixatoms[randatom][j] == 0)
						dq[j] = delx*(dsfmt_genrand_close_open(&dsfmt)-0.5);
					else
						dq[j] = 0.000;

					PosIons[randatom][j] = PosIons[randatom][j] + dq[j] ;

				}

				//calculate energy

				//fprintf(pFile, "%d   atom   %d  %d  %11.5f %11.5f %11.5f", i, randatom, randatom, dq[0], dq[1], dq[2]);

				if(dq[0] !=0.000 || dq[1] !=0.000 || dq[2] !=0.000)
				{

					//Calculate energy to due to this particular atom

					delPot = CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, randatom, 0, 3) - delPot;


					//fprintf(pFile, "  %11.5f ", delPot);


					if(delPot<=0) //accept move
					{
						Pot = Pot+delPot;
						n_accp = n_accp +1;
						accp_rej = 1;
					}
					else
					{
						//calculate boltzman probablity
						prob = exp(-delPot*1e6/KB/AVG/Temp);
						randnum=dsfmt_genrand_close_open(&dsfmt);

						//accept if probality is greater than random number
						if(prob>=randnum)
						{

							Pot = Pot+delPot;
							n_accp = n_accp +1;		
							accp_rej = 1;
						}	
						else
						{
							n_rej = n_rej +1;
							accp_rej = 0;
							for(int j=0;j<3;j++)
							{
								PosIons[randatom][j] = PosIons[randatom][j] - dq[j] ;
							}
						}
					}

				}


			}
			else
			{
				//pick two random atoms to swap


				atom1 = natoms_type[0]  + dsfmt_genrand_close_open(&dsfmt)*(natoms_type[1]);
				atom2 = natoms_type[0] + natoms_type[1] + dsfmt_genrand_close_open(&dsfmt)*(natoms_type[2]);

				if(fixatoms[atom1][0] != 0 && fixatoms[atom1][1] != 0 && fixatoms[atom1][2] != 0)
				{
					//	cout<<"Atom "<<atom1<<"  is fixed"<<endl;
					continue;
				}

				if(fixatoms[atom2][0] != 0 && fixatoms[atom2][1] != 0 && fixatoms[atom2][2] != 0)
				{
					//	cout<<"Atom "<<atom2<<"  is fixed"<<endl;
					continue;
				}




				delPot = 0;

				delPot = CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, atom1, 0, 3);

				delPot = delPot + CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, atom2, 1, atom1);


				for(int j=0;j<3;j++)
				{
					//swap positions of the two atoms
					dq[j] = PosIons[atom1][j];
					PosIons[atom1][j] = PosIons[atom2][j];
					PosIons[atom2][j] = dq[j] ;
				}


				delPot = CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, atom1, 0, 3) - delPot;

				delPot = delPot + CalEner1(PosIons, natoms, boxcell, batom1, batom2, bondpot, nbatom1, nbatom2, nonbondpot, bondpar, nonbondpar, pairs, atom2, 1, atom1);

				if(delPot<=0) //accept move
				{
					Pot = Pot + delPot;
					n_accp = n_accp +1;
					n_swap_accp = n_swap_accp + 1;
					accp_rej = 1;
				}
				else
				{
					//calculate boltzman probablity
					prob = exp(-delPot*1e6/KB/AVG/Temp);
					randnum=dsfmt_genrand_close_open(&dsfmt);

					//accept if probality is greater than random number
					if(prob>=randnum)
					{
						Pot = Pot + delPot;
						n_accp = n_accp +1;
						n_swap_accp = n_swap_accp + 1;	
						accp_rej=1;
					}	
					else
					{
						n_rej = n_rej +1;
						n_swap_rej = n_swap_rej + 1;
						accp_rej = 0;

						for(int j=0;j<3;j++)
						{
							//swap positions of the two atoms
							PosIons[atom2][j] = PosIons[atom1][j];
							PosIons[atom1][j] = dq[j] ;

						}
					}
				}

			}



			if(i%(int(printmc*nrun))==0 && r !=0)
			{
					print_coor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 1, i,'a', "XDATCAR");

					cout<<"Finished MC sweep no.: "<<int(float(i)/nrun)<<endl;
					cout<<"Percentage moves accepted: "<< 100*float(n_accp)/(n_accp+n_rej)<<endl;

				wrap_coor(PosIons, natoms, boxcell);

				if(i !=0 && i%(printmc*nrun*10)==0)
					print_carcoor(PosIons, natoms, boxcell,  n_atomtype, natoms_type, atomtype, 0, i,'w', "CONTCAR");

			}
			else if(i%int((printmc*nrun))==0 && r ==0 )
			{

				cout<<"Finished equibration MC step no.: "<<int(float(i)/(nrun))<<endl;
				cout<<"Percentage moves accepted: "<< 100*float(n_accp)/(n_accp+n_rej)<<endl;
			}


		}

	}

	cout<<"Total no. of accepted MC moves: "<<n_accp<<endl;
	cout<<"Total no. of rejected MC moves: "<<n_rej<<endl;
	cout<<"Total no. of swap move accepted: "<<n_swap_accp<<endl;
	cout<<"Total no. of swap move rejected: "<<n_swap_rej<<endl;
	cout<<"Total percentage of accepted moves: "<<100*float(n_accp)/(n_accp+n_rej)<<endl;
	cout<<"Total percentage of accepted atomic moves: "<<100*float(n_accp-n_swap_accp)/(n_accp+n_rej-n_swap_accp-n_swap_rej)<<endl;
	cout<<"Total percentage of accepted swap moves: "<<100*float(n_swap_accp)/(n_swap_accp+n_swap_rej)<<endl;


}

void wrap_coor(double **PosIons, int natoms, float **boxcell)
{

	for(int i=0;i<natoms;i++)
	{

		for(int j=0;j<3;j++)
		{

			PosIons[i][j] = PosIons[i][j] - floor(PosIons[i][j]/boxcell[j][j])*boxcell[j][j];

		}



	}


}
