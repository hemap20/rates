//***********************************************
//function to implement steepest descent method!! 
//***********************************************

#include "libinclude.h" 
#include "mathdec.h"
#include "potdec.h"
#include "fundec.h"
#include "const.h"

void minimization(double **PosIons, int natoms, double **ForceIons, float **boxcell, double Pot, string *type, int **fixatoms)
{

	cout<<"*****Entering the minimization code*******"<<endl<<endl;

	char Prec;
	float c1_wolfe, alpha, brac_par;
	int itermax;
	int wolfe_cond;
	int algo;
	int hesupdt;

	char garbage[100];

	ifstream InputIn("min.in",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	InputIn>>garbage>>garbage;
	InputIn>>Prec;

	InputIn>>garbage>>garbage;
	InputIn>>c1_wolfe;

	InputIn>>garbage>>garbage;
	InputIn>>itermax;

	InputIn>>garbage>>garbage;
	InputIn>>alpha;
	
	InputIn>>garbage>>garbage;
	InputIn>>brac_par;

	InputIn>>garbage>>garbage;
	InputIn>>wolfe_cond;

	InputIn>>garbage>>garbage;
	InputIn>>algo;

	InputIn>>garbage>>garbage;
	InputIn>>hesupdt;

	InputIn.close();



	float maxF, delPos, delV;   // parameters for stopping the code

	///////////////////////////////////////////////////

	switch (Prec)
	{
		case 'L': maxF=0.1*KJEV;              //KJ/angs 
			  delV =0.1*KJEV;          //KJ
			  delPos = 0.01;	//angs
			  break;
		case 'M': maxF=0.05*KJEV; 
			  delV=0.01*KJEV;
			  delPos = 0.001;
			  break;
		case 'H': maxF=0.01*KJEV; 
			  delV=0.001*KJEV; 
			  delPos = 0.0001;
			  break;
		default:  maxF=0.01*KJEV; 
			  delV=0.001*KJEV; 
			  delPos = 0.0001;
			  break;

	}
	//////////////////////////////////////////////////


	double alpha1=alpha;
	double normf = 0;
	double maxf = 1;
	double normx = 0;
	double delPot = 0;
	double **PrevPos;	

	maxf = getMatabsmax(ForceIons, natoms, 3);

	//print the progress 
	print_prog(maxf, normx, delPot, maxF, delPos, delV);
	int iter = 0;


	PrevPos=new double * [natoms];

	for(int i=0;i<natoms;i++)
	{
		PrevPos[i]=new double [3];
	}


	//set hessian to identity matrix!!! 




	while( maxf > maxF ||  abs(delPot) > delV)
	{
		cout<<"Starting Minimization Iteration Number: "<<iter +1<<endl<<endl;

		mat_equate(PosIons, PrevPos, natoms, 3);
	
		alpha1 = alpha;

		//get direction to udptate

		if(algo==1 || iter==0)
		{
			cout<<"Using the Steepest descent step for minimization"<<endl;
		}


		for(int i=0;i<natoms;i++)
		{
			for(int j=0;j<3;j++)
				PosIons[i][j]= PosIons[i][j] + alpha1*ForceIons[i][j]; //updating the postions
		}


		if(wolfe_cond==1 || iter==0)	
		{
			pair<double, float> answer  = first_wolfe_brac(PosIons, PrevPos, ForceIons, natoms,boxcell, Pot, c1_wolfe, alpha1, brac_par,fixatoms);
			delPot = answer.first;
			alpha1 = answer.second;
			Pot = Pot + delPot;
		}

		iter=iter+1;


		cout<<"Using Step length as: "<<alpha1<<endl;

		printCoor(PosIons, natoms, type);
		printFor(ForceIons,natoms,type);

		maxf = getMatabsmax(ForceIons, natoms, 3);


		print_prog(maxf, normx, delPot, maxF, delPos, delV);

		if(iter>itermax)
		{
			cout<<endl<<"Maximum iteration steps reached"<<endl;
			break;}


	}

	if(iter<itermax)
	{
		cout<<endl<<"Reached required accuracy, stopping minimization algorithm"<<endl;
	}
	for(int i=0;i<natoms;i++)
	{
		delete [] PrevPos[i]; 
	}


	delete [] PrevPos;

}

pair <double, float>  first_wolfe_brac(double **PosIons, double **PrevPos,double **ForceIons, int natoms,  float **boxcell, double Pot, float c1_wolfe, float alpha, float brac_par, int **fixatoms)
{

	cout<<"Checking if first wolfe's condition is satisfied"<<endl;

	double prevPot, prev_dot, forc_dot;
	double tmp;
	double **PrevForce;
	int iter =0;
	prevPot = Pot;
	string garbage;

	PrevForce=new double * [natoms];
	for(int i=0; i<natoms;i++)
	{
		PrevForce[i]=new double [3];
	}

	mat_equate(ForceIons, PrevForce, natoms, 3);

	ifstream InputIn("molecules.out",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}
	
	int n_molecules;	

	InputIn>>n_molecules>>garbage;

	int *molecule_atom1, *molecule_atom2;
	molecule_atom1=new int [n_molecules];
	molecule_atom2=new int [n_molecules];


	for(int i=0;i<n_molecules;i++)
	{
		InputIn>>molecule_atom1[i]>>garbage>>molecule_atom2[i];
	}

	InputIn.close();


	forc_dot = sum2mat(ForceIons, natoms, 3);
	tmp = prevPot - c1_wolfe*alpha*forc_dot;
	prev_dot = forc_dot;

	zeros(natoms, 3 ,ForceIons);
	Pot = CalEnerFor(PosIons, natoms, boxcell, ForceIons, "bondlist.out", true, true);
	Pot = Pot + CalEnerFor(PosIons, natoms, boxcell, ForceIons, "nonbondlist.out", false, true);


	fixatoms_forces(natoms, ForceIons, fixatoms);

	while(Pot> tmp || isinf(Pot)==1 || Pot != Pot)
	{


		if(iter==0)
			cout<<"Starting bracktracing iteration"<<endl;

		alpha=alpha*brac_par;

		cout<<"Iteration "<<iter<<": Changing step length to  "<<alpha<<endl;

		for(int i=0;i<natoms;i++)
		{
			for(int j=0;j<3;j++)
				PosIons[i][j]= PrevPos[i][j] + alpha*PrevForce[i][j];
		}


		zeros(natoms, 3 ,ForceIons);
		Pot = CalEnerFor(PosIons, natoms, boxcell, ForceIons, "bondlist.out", true, true);
		Pot = Pot + CalEnerFor(PosIons, natoms, boxcell, ForceIons, "nonbondlist.out", false, true);

		fixatoms_forces(natoms, ForceIons, fixatoms);

		tmp = prevPot - c1_wolfe*alpha*forc_dot;


		if(iter==20)
			break;

		iter=iter+1;
	}

	if(iter<20)
		cout<<"******First Wolfe's condition satisfied*****"<<endl;
	else
		cout<<"Maximum iterations in bracktracing method reached"<<endl;

	cout<<endl<<"Total Potential Energy of the system:  "<<Pot/KJEV<<" eV "<<endl;

	tmp = Pot - prevPot;

	return make_pair(tmp, alpha);

}
