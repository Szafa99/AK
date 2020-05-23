#include <stdio.h>
#include <math.h>
#include "eval_time.h"

#define SIZE 1024			//rozmiar glownych macierzy

//#define nb 32				//rozmiar bloku alg. blokowego

static double a[SIZE*SIZE];
static double tr_matrix[SIZE*SIZE];
static double b[SIZE*SIZE];
static double c[SIZE*SIZE];
static double d[SIZE*SIZE];
static double e[SIZE*SIZE];

void dgemm_naive_trans(int,double*,double*,double*);
void dgemm_blocked_trans(int,double*,double*,double*);
void dgemm_unroll4_trans(int,double*,double*,double*);
void block_trans(int,int,int,int,double*,double*,double*);
void transpose(int,double*,double*);
void showresult();

unsigned int blocksize;

double suma_b[6]={0,0,0,0,0,0};
double suma_n=0;
double suma_un=0;


int main(void)
{
	unsigned int i,j,n,f;
	double time_tabl[3];


for(int z=0;z<5;z++)
{

	for( int size = 256; size <= 256; size*=2 ) {
		n = size;
		f = 2*n*n*n;		//liczba operacji zmiennoprzecinkowych

		for(i=0;i<n;++i)	//wypelnij a i b jakimis wartosciami poczatkowymi
		    for (j=0;j<n;++j)
		    {
			a[j+i*n]=(double)(i+j*n);
			b[j+i*n]=(double)(j+i*n);
		    }
		//c i d zostaw wyzerowane
/////////////////////////////////////////////NAIVETRANS



		init_time();
		transpose(n,a,tr_matrix); 		
		dgemm_naive_trans(n,tr_matrix,b,c);
		read_time(time_tabl);
		suma_n+=(double)f/time_tabl[1];
		printf("NAIVE:\tMX_SIZE=%d time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );


// //////////////////////////////////// UNROLL4TRANS

		init_time();
		transpose(n,a,tr_matrix); 		
		dgemm_unroll4_trans(n,tr_matrix,b,d);
		read_time(time_tabl);
		suma_un+=(double)f/time_tabl[1];
		printf("UNROLL:\tMX_SIZE=%d time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );
		//sprawdzenie czy oba algorytmy daly ten sam wynik

		for (i=0;i<n*n;++i) 
		    if (fabs(c[i]-d[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}

//////////////////////////////////// BLOCKED TRANS
int count=0;
		for( int block = 4; block <= 1024; block*=2 )
			if( block <= 128 ) {
				blocksize = block;

				init_time();
				transpose(n,a,tr_matrix); 		
				dgemm_blocked_trans(n,tr_matrix,b,e);
				read_time(time_tabl);
                suma_b[count++]+=(double)f/time_tabl[1];
				printf("BLOCKED\tMX_SIZE=%d  BL_SIZE%d time = %2.3lf GFLOPS = %.3lf\n", n, block, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );

		//sprawdzenie czy oba algorytmy daly ten sam wynik
					for (i=0;i<n*n;++i) 
					if (fabs(d[i]-e[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}
					else e[i] = 0.0;
		}
		

        for (i=0;i<n*n;++i) { c[i] = 0.0; d[i] = 0.0; } // zeruj macierz c i d

	}
}

	showresult();
	rtrn:
	return(0);
}


void showresult(){

printf("Srednia GFPS T_naive = %.3lf\n", suma_n/5.0e9 );
printf("Srednia GFPS T_UNROLL = %.3lf\n", suma_un/5.0e9 );

		register int blocksize=2;
		for(int i=0;i<6;i++)
				printf("Srednia GFPS T_BLOCKED dla bloku %d = %.3lf\n",blocksize*=2, suma_b[i]/5.0e9 );
    
}



void dgemm_unroll4_trans(int n, double* A, double* B, double* C)
{
register double reg0,reg1,reg2,reg3;

for(register int i=0;i<n;i+=2)
    for(register int j=0;j<n;j+=2)
    {
	reg0=reg1=reg2=reg3=0.0;
	
	for(register int k=0;k<n;++k)
	{
	    reg0+=A[k+i*n]*B[k+j*n];
	    reg1+=A[k+(i+1)*n]*B[k+j*n];
	    reg2+=A[k+i*n]*B[k+(j+1)*n];
	    reg3+=A[k+(i+1)*n]*B[k+(j+1)*n];
	}
	
	C[i+j*n]+=reg0;
	C[i+1+j*n]+=reg1;
	C[i+(j+1)*n]+=reg2;
	C[i+1+(j+1)*n]+=reg3;
    }
}






void dgemm_naive_trans(int n, double* A, double* B, double* C)
{
register double cij;

for(register int i=0;i<n;++i)
    for(register int j=0;j<n;++j)
    {
	cij=C[i+j*n]; 			
	for(register int k=0;k<n;++k)
	    cij+=A[k+i*n]*B[k+j*n]; 	 
	C[i+j*n]=cij; 			
    }
}




void dgemm_blocked_trans(int n, double* A, double* B, double* C)
{

for( register int i=0;i<n;i+=blocksize)
    for(register int j=0;j<n;j+=blocksize)
	for(register int z=0; z<n; z+=blocksize)
	    block_trans(n,i,j,z,A,B,C);
}


void block_trans(int n, int bi, int bj, int bz, double *A, double *B, double *C)
{

register double cij;

for(register int i=bi;i<bi+blocksize;++i)
    for(register int j=bj;j<bj+blocksize;++j)
    {
	cij=C[i+j*n];
	for(register int z=bz;z<bz+blocksize;++z)
	    cij+=A[z+i*n]*B[z+j*n];
	C[i+j*n]=cij;
    }
}


void transpose(int n, double*X, double*Y )
{
 register int i, j;

 for(i=0;i<n;i++)
   for(j=0;j<n;j++)
     Y[i*n+j] = X[i+j*n];
}
