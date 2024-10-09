//****************************************
// various math utlities
// ***************************************

#include "libinclude.h" 
#include "const.h"


void res_vel(double **vel, int natoms, float Temp, float *mass);
double gasdev(double randnum);


//function to get maximum value from an array
double getmax(double *val, int len)
{

	double a =val[0];
	for(int i=0;i<len;i++)
	{
		if(val[i]>= a)
			a =val[i];



	}


	return a;

}

double getabsmax(double *val, int len)
{

	double a =abs(val[0]);
	for(int i=0;i<len;i++)
	{
		if(abs(val[i])>= a)
			a =abs(val[i]);



	}


	return a;

}

//function to get 2-norm
double getnorm(double *val, int len)
{
	double norm = 0;

	for(int i=0;i<len;i++)
	{
		norm = norm + pow(val[i],2);

	}

	norm = sqrt(norm);

	return norm;

}

//get the maximum value in a matrix

double getMatmax(double **mat, int rows, int columns)
{

	double a =mat[0][0];

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns;j++)
		{
			if(mat[i][j]>= a)
				a =mat[i][j];
		}

	}

	return a;

}

//get the maximum absolute value in a matrix

double getMatabsmax(double **mat, int rows, int columns)
{
	double a =mat[0][0];

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns;j++)
		{
			if(abs(mat[i][j])>= a)
				a = abs(mat[i][j]);
		}

	}

	return a;
}

//function to get 2-norm
double getMatRMS(double **mat, int rows, int columns)
{
	double norm = 0;

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns; j++)
		{
			norm = norm + pow(mat[i][j],2);
		}

	}

	norm = sqrt(norm)/rows/columns;

	return norm;

}

double sum2mat(double **mat, int rows, int columns)
{
	double sum2= 0;

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns; j++)
		{
			sum2 = sum2 + pow(mat[i][j],2);
		}

	}

	return sum2;

}

//////////////////////matrix and vector operations/////////////////

void vec_equate(double *vec1, double *vec2, int len)
{
	////////equate two vectors

	for(int i=0;i<len;i++)
	{
		vec2[i]=vec1[i];
	}

}


void mat_equate(double **mat1, double **mat2, int rows, int columns)
{
	//equate two matrices

	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns;j++)
		{
			mat2[i][j]=mat1[i][j];		

		}	

	}

}


void invmat(double **mat, int rows, double **matinv)
{

	MatrixXd matrix(rows, rows);
	for(int i=0;i<rows;i++)
	{
		for(int j=0; j<rows;j++)
		{
			matrix(i,j) = mat[i][j];

		}

	}

	MatrixXd inv = matrix.inverse();

	for(int i=0;i<rows;i++)
	{
		for(int j=0; j<rows;j++)
		{
			matinv[i][j] = inv(i,j);

		}

	}

}

double vec_dot(double *vec1, double *vec2, int len)
{
	double vecdot=0;

	for(int i=0;i<len;i++)
		vecdot+=  vec1[i]*vec2[i];

	return vecdot;

}

///vector multiplication with its transpose!!!
void  vec_X_vecT(double *vec1, double *vec2, int len1, int len2, double **matres )
{
	for(int i=0;i<len1;i++)
	{	
		for(int j=0;j<len2;j++)
		{
			matres[i][j] = vec1[i]*vec2[j];
		}
	}
}


void eye(int rows, int columns, double **mat)
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns; j++)
		{
			if(i==j)
				mat[i][j] = 1;
			else
				mat[i][j] = 0;
		}
	}
}

void zeros(int rows, int columns, double **mat)
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns; j++)
		{
			if(i==j)
				mat[i][j] = 0.000;
			else
				mat[i][j] = 0.000;
		}
	}
}

void zeros_int(int rows, int columns, int **mat)
{
	for(int i=0;i<rows;i++)
	{
		for(int j=0;j<columns; j++)
		{
			if(i==j)
				mat[i][j] = 0.000;
			else
				mat[i][j] = 0.000;
		}
	}
}

void vec_zero(double *vec, int len)
{
	for(int i=0;i<len;i++)
		vec[i]=0;

}

void vec_zero_int(int *vec, int len)
{
	for(int i=0;i<len;i++)
		vec[i]=0;

}



int vec_dot_mat(double *vec, int veclen, double **mat, int rows, int columns, double *vecres)
{

	if(veclen != rows)
	{
		cerr<<"Error in vector.dot.matrix, matrix dimensions must agree";
		return 0;
	}

	for(int i=0;i<columns; i++)
	{
		vecres[i] = 0;
		for(int j=0;j<rows;j++) 
			vecres[i] += vec[j]*mat[j][i];


	}

	return columns;

}

int mat_X_vec(double *vec, int veclen, double **mat, int rows, int columns, double *vecres)
{

	if(veclen != columns)
	{
		cerr<<"Error in matrix.X.vector, matrix dimensions must agree";
		return 0;
	}

	for(int i=0;i<rows; i++)
	{
		vecres[i] = vec_dot(mat[i], vec, veclen);


	}

	return rows;

}

void mat_transpose(double **mat, int rows, int columns, double **mat_trans)
{

	for(int i=0;i<rows;i++)
	{
		for(int j=0; j<columns;j++)
		{

			mat_trans[j][i]	 = mat[i][j];	
		}
	}
}





void mat_X_mat(double **mat1, int rows1, int columns1, double **mat2, int rows2, int columns2, double **matres)
{

	for(int i=0;i<rows1;i++)
	{

		for(int k=0;k<columns2;k++)
		{
			matres[i][k] = 0;

			for(int j=0;j<rows2;j++)
			{
				matres[i][k]+= mat1[i][j]*mat2[j][k];

			}		

		}

	}


}

void matT_X_mat(double **mat1, int rows1, int columns1, double **mat2, int rows2, int columns2, double **matres)
{

	for(int i=0;i<columns1;i++)
	{
		for(int k=0;k<columns2;k++)
		{
			matres[i][k] = 0;

			for(int j=0;j<rows2;j++)
			{
				matres[i][k]+= mat1[j][i]*mat2[j][k];

			}		

		}

	}


}

/////printing functions
void print_vec(double *vec, int len, string str)
{
	cout<<endl;
	cout<<str<<"=  [";


	for(int i =0;i < len;i++)
	{
		if(i!=len-1)
			cout<<vec[i]<<endl;		
		else if(i==len-1)
			cout<<vec[i]<<"]"<<endl;

	}
}

void print_mat(double **mat, int rows, int columns, string str)
{
	cout<<endl;

	cout<<str<<"=\n  [";

	for(int i=0;i<rows;i++)
	{

		for(int j = 0; j<columns; j++)
		{
			if(j ==columns-1)
				cout<<mat[i][j]<<",";
			else 	
				cout<<mat[i][j]<<"\t";
		}

		if(i!=rows-1)
			cout<<"\n";


	}	

	cout<<"]"<<endl;

}


void print_mat_int(int **mat, int rows, int columns, string str)
{
	cout<<endl;

	cout<<str<<"=\n  [";

	for(int i=0;i<rows;i++)
	{

		for(int j = 0; j<columns; j++)
		{
			if(j ==columns-1)
				cout<<mat[i][j]<<",";
			else 	
				cout<<mat[i][j]<<"\t";
		}

		if(i!=rows-1)
			cout<<"\n";


	}	

	cout<<"]"<<endl;

}


void print_pos_forces(double *Pos, double *Jac, double len)
{

	cout<<endl<<left<<setw(40)<<"Postions"<< setw(40)<<"Forces"<<endl;
	cout<<setfill('-')<<setw(80)<<"-"<<endl;
	cout<<setfill(' ');
	cout<<left<<setw(20)<<"x"<<setw(20)<<"y"<<setw(20)<<"dx"<<setw(20)<<"dy"<<endl;
	cout<<left<<setw(20)<<Pos[0]<<setw(20)<<Pos[1]<<setw(20)<<Jac[0]<<setw(20)<<Jac[1]<<endl;
	cout<<setfill('-')<<setw(80)<<"-"<<endl;
	cout<<setfill(' ')<<endl;


}


void print_prog(double maxf, double normx, double delPot, float maxF, float delPos, float delV )
{
	cout<<endl<<left<<setw(20)<<" "<<setw(20)<<"Current"<< setw(20)<<"Stopping Criteria"<<endl;
	cout<<setfill('-')<<setw(80)<<"-"<<endl;
	cout<<setfill(' ');
	cout<<left<<setw(20)<<"Max. Force"<<setw(20)<<maxf<<setw(20)<<maxF<<endl;
	cout<<left<<setw(20)<<"Del Position"<<setw(20)<<normx<<setw(20)<<delPos<<endl;
	cout<<left<<setw(20)<<"Del Potential"<<setw(20)<<delPot<<setw(20)<<delV<<endl;
	cout<<setfill('-')<<setw(80)<<"-"<<endl;
	cout<<setfill(' ')<<endl;
}




int factorial(int num)
{

	int fact=1;

	for(int i=1;i<=num;i++)
	{
		fact = fact*i;


	}


	return fact;

}

void ini_mass(float *mass, int natoms, string *type)
{
	for(int i=0;i<natoms; i++)
	{
		if(type[i].compare("Cu") ==0)
			mass[i]	= 63.546;
		else if(type[i].compare("C") ==0)
			mass[i] = 18.0;
		else if(type[i].compare("O") ==0)
			mass[i] = 44.0;
		else if(type[i].compare("H") ==0)
			mass[i] = 1.00794;

	}

}


void ini_vel(double **vel, int natoms, float Temp, float *mass, int numrand, int **fixatoms)
{
	double sumV[3];
	sumV[0] = 0.00; sumV[1] = 0.0000 ; sumV[2] = 0.0000 ;
	double sumMass=0.000;
	double scal_fac;
	srand(time(NULL)*numrand);
	double randnum;

	for(int i=0; i<natoms; i++)
	{

		scal_fac = sqrt(KB*AVG*Temp/mass[i]);
		for(int j=0;j<3;j++)
		{
			randnum = rand();

			if(fixatoms[i][j] == 0)
			{
				vel[i][j] = gasdev(randnum);

				vel[i][j] = vel[i][j]*scal_fac;

			}
			else
				vel[i][j] == 0.000;

			
			sumV[j] = sumV[j] + mass[i]*vel[i][j];

		}


			sumMass = sumMass + mass[i];
	}

	for(int j=0;j<3;j++)
	{
		sumV[j] = sumV[j]/sumMass ;
	}


	for(int i=0; i<natoms; i++)
	{
		for(int j=0;j<3;j++)
		{
			if(fixatoms[i][j] == 0)
				vel[i][j] = (vel[i][j] - sumV[j]);  
		}
	}


	res_vel(vel, natoms, Temp, mass);

}

//Box muller method to get normal distribution
double gasdev(double randnum)
{
	srand(randnum);
	static int  iset=0;  //static variable once defined remains there till the end of the program
	static float gset;
	double fac, rsq, v1, v2;

	if (iset==0) {
		do {
			v1 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			v2 = 2.0 * rand() / double(RAND_MAX) - 1.0;
			rsq = v1 * v1 + v2 * v2;
		} while (rsq >= 1.0 || rsq == 0.0);
		fac = sqrt(-2.0 * log(rsq)/rsq);
		gset = v1 * fac;
		iset = 1;
		return v2*fac;
	} else {
		iset = 0;
		return gset;
	}
}


void res_vel(double **vel, int natoms, float Temp, float *mass)
{

	double sumV2 = 0.0000;
	double scal_fac = 0.0000;

	for(int i=0; i<natoms; i++)
	{
		for(int j=0;j<3;j++)
		{
			sumV2 = sumV2 + mass[i]*vel[i][j]*vel[i][j];
		}
	}


	scal_fac = sqrt((3*natoms-3)*KB*AVG*Temp/sumV2);

	for(int i=0; i<natoms; i++)
	{
		for(int j=0;j<3;j++)
		{
			vel[i][j] = vel[i][j]*scal_fac*1e-5;     //velocities are obtained in angs/fs
		}
	}


}

void readvel(double **vel, int natoms)
{
	ifstream InputIn("velocities.in",ios::in);
	if(!InputIn)
	{
		cerr << "File InputIn could not be opened" <<endl;
		exit(1);
	}

	for(int i=0; i<natoms; i++)
	{
		for(int j=0;j<3;j++)
		{
			InputIn>>vel[i][j];
		}
	}

	InputIn.close();

}





//get center of mass

void getCOM(double **PosIons, int natoms, int natom1, int natom2, float **box, float *mass, double *RCOM)
{

	double R1[3], R2[3];
	double dx[3];

	for(int i=0;i<3;i++)
	{
		R1[i] = PosIons[natom1][i];
		R2[i] = PosIons[natom2][i];
		dx[i] = R1[i]-R2[i];	
	}

	R2[0] = R2[0] + box[0][0]*ceil(dx[0]/box[0][0]-0.5) + box[1][0]*ceil(dx[1]/box[1][1]-0.5) + box[2][0]*ceil(dx[2]/box[2][2]-0.5);
	R2[1] = R2[1] + box[0][1]*ceil(dx[0]/box[0][0]-0.5) + box[1][1]*ceil(dx[1]/box[1][1]-0.5) + box[2][1]*ceil(dx[2]/box[2][2]-0.5);
	R2[2] = R2[2] + box[0][2]*ceil(dx[0]/box[0][0]-0.5) + box[1][2]*ceil(dx[1]/box[1][1]-0.5) + box[2][2]*ceil(dx[2]/box[2][2]-0.5);


	for(int i=0;i<3;i++)
	{
		RCOM[i] = (mass[natom1]*R1[i]+ mass[natom2]*R2[i])/(mass[natom1]+mass[natom2]);
	}

}


// dot product of two vectors from position matrix

double DotVecPos(double **PosIons, int atom1, int atom2)
{

	double prod = 0.00;
	for(int i=0; i<3;i++)
	{
		prod = prod + PosIons[atom1][i]*PosIons[atom2][i];
	
	}

	return prod;

}



int chkMol(double **PosIons, int natoms, int n_molecules,  int *molecule_atom1,int *molecule_atom2, int  fixatom1, int fixatom2, float Rchkmol)
{

	int chmol=1;
	for(int i=0;i<n_molecules;i++)
	{
		if(PosIons[molecule_atom1[i]][0]>Rchkmol && molecule_atom1[i]!=fixatom1)
		{
			chmol=0;
			cout<<"Atom "<<molecule_atom1[i]<<"has x-coordinate("<<PosIons[molecule_atom1[i]][0] <<") greater than "<< Rchkmol<<endl;

		}
			
	
	}


	return chmol;
}



