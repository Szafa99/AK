
#include <stdio.h>
#include <math.h>

#include "eval_time.h"

#define SIZE 2048			//rozmiar glownych macierzy


static double a[SIZE*SIZE];
static double b[SIZE*SIZE];
static double c[SIZE*SIZE];
static double d[SIZE*SIZE];

static double e[SIZE*SIZE];

void dgemm_naive(int,double*,double*,double*);
void dgemm_blocked(int,double*,double*,double*);
void dgemm_unroll4(int,double*,double*,double*);
void block(int,int,int,int,double*,double*,double*);
void showresult();

unsigned int blocksize;

double suma_b[6]={0,0,0,0,0,0};
double suma_n=0;
double suma_un=0;


int main(void)
{
	unsigned int i,j,n,f,size;
	double time_tabl[3];


for(int z=0;z<5;z++){
	for( int size = 512; size <= 512; size*=2 ) {
		n = size;
		f = 2*n*n*n;		//liczba operacji zmiennoprzecinkowych

		for(i=0;i<n;++i)	//wypelnij a i b jakimis wartosciami poczatkowymi
		    for (j=0;j<n;++j)
		    {
			a[j+i*n]=(double)(i+j*n);
			b[j+i*n]=(double)(j+i*n);
		    }

		//c i d zostaw wyzerowane
/////////////////////////////////////////////////// NAIVE
		init_time();
		dgemm_naive(n,a,b,c);
		read_time(time_tabl);
		suma_n+=(double)f/time_tabl[1];
		printf("UNROLL:\tMX_SIZE=%d time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[2]/1.0e9 );


//////////////////////////// UNROLL
		init_time();
		dgemm_unroll4(n,a,b,d);
		read_time(time_tabl);
		suma_un+=(double)f/time_tabl[1];
		printf("UNROLL:\tMX_SIZE=%d time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[2]/1.0e9 );

		//sprawdzenie czy oba algorytmy daly ten sam wynik

		for (i=0;i<n*n;++i) 
		    if (fabs(c[i]-d[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}

// //////////////////////////// BLOCK
int count = 0;
		for( blocksize = 4; blocksize <= 128; blocksize*=2 )
			if( blocksize < n ) {

				init_time();
				dgemm_blocked(n,a,b,e);
				read_time(time_tabl);
  				suma_b[count++]+=(double)f/time_tabl[1];
				printf("BLOCKED\tMX_SIZE=%d  BL_SIZE%d time = %2.3lf GFLOPS = %.3lf\n", n, blocksize, time_tabl[1],(double)f/time_tabl[2]/1.0e9 );

		
				// sprawdzenie czy oba algorytmy daly ten sam wynik
				for (i=0;i<n*n;++i) 
					if (fabs(d[i]-e[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}
					else e[i] = 0.0;

		}
		for (i=0;i<n*n;++i) { c[i] = 0.0; d[i] = 0.0; } 
	}
	}

	showresult();
	rtrn:
	return(0);
}

void showresult(){

printf("Srednia GFPS naive = %.3lf\n", suma_n/5.0e9 );
printf("Srednia GFPS UNROLL = %.3lf\n", suma_un/5.0e9 );

		register int blocksize=2;
		for(int i=0;i<6;i++)
				printf("Srednia GFPS BLOCKED dla bloku %d = %.3lf\n",blocksize*=2, suma_b[i]/5.0e9 );
    
}

void dgemm_unroll4(int n, double* A, double* B, double* C)
{
register int i,j,k;
register double reg0,reg1,reg2,reg3;

for(i=0;i<n;i+=2)
    for(j=0;j<n;j+=2)
    {
	reg0=reg1=reg2=reg3=0.0;
	
	for(k=0;k<n;++k)
	{
	    reg0+=A[i+k*n]*B[k+j*n];
	    reg1+=A[i+1+k*n]*B[k+j*n];
	    reg2+=A[i+k*n]*B[k+(j+1)*n];
	    reg3+=A[i+1+k*n]*B[k+(j+1)*n];
	}
	
	C[i+j*n]+=reg0;
	C[i+1+j*n]+=reg1;
	C[i+(j+1)*n]+=reg2;
	C[i+1+(j+1)*n]+=reg3;
    }
}


void dgemm_naive(int n, double* A, double* B, double* C)
{
register int i,j,k;
register double cij;

for(i=0;i<n;++i)
    for(j=0;j<n;++j)
    {
	cij=C[i+j*n]; 			// cij = C[i][j]
	for(k=0;k<n;++k)
	    cij+=A[i+k*n]*B[k+j*n]; 	// cij += A[i][k]*B[k][j]
	C[i+j*n]=cij; 			// C[i][j] = cij
    }
}

void dgemm_blocked(int n, double* A, double* B, double* C)
{
register int bi,bj,bk;

for(bi=0;bi<n;bi+=blocksize)
    for(bj=0;bj<n;bj+=blocksize)
	for(bk=0;bk<n;bk+=blocksize)
	    block(n,bi,bj,bk,A,B,C);
}

void block(int n, int bi, int bj, int bk, double *A, double *B, double *C)
{

register int i,j,k;
register double cij;

for(i=bi;i<bi+blocksize;++i)
    for(j=bj;j<bj+blocksize;++j)
    {
	cij=C[i+j*n];
	for(k=bk;k<bk+blocksize;++k)
	    cij+=A[i+k*n]*B[k+j*n];
	C[i+j*n]=cij;
    }
}
