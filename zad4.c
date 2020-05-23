

#include <stdio.h>
#include <math.h>
#include "eval_time.h"

#define SIZE 1024			//rozmiar glownych macierzy


static double a[SIZE*SIZE];
static double b[SIZE*SIZE];
static double c[SIZE*SIZE];
static double d[SIZE*SIZE];
static double e[SIZE*SIZE];

void dgemm_naive(int,double*,double*,double*);
void dgemm_blocked(int,double*,double*,double*);
void dgemm_unroll16(int,double*,double*,double*);
void block(int,int,int,int,double*,double*,double*);
void showresult();

static double suma_u16[4];
static double suma_n[4];
static double suma_b4[4];


unsigned int blocksize=4;

int main(void)
{
	unsigned int i,j,n,f;
	double time_tabl[3];
int count=0;

for(int z=0;z<5;z++){
	count=0;	
	for(int size = 128; size <= 1024; size *=2 ) {
		
		n = size;
		f = 2*n*n*n;		//liczba operacji zmiennoprzecinkowych

			for(i=0;i<n;++i)	
				for (j=0;j<n;++j)
				{
				a[j+i*n]=(double)(i+j*n);
				b[j+i*n]=(double)(j+i*n);
				}
			//c i d zostaw wyzerowane

	////////////////////////////////////////// Naive
			init_time();
			dgemm_naive(n,a,b,c);
			read_time(time_tabl);
			suma_n[count]+=(double)f/time_tabl[1];
			printf("Naive:MX_SIZE=%d \ttime = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );


	//////////////////// UNROLL_16
			init_time();
			dgemm_unroll16(n,a,b,d);
			read_time(time_tabl);
			suma_u16[count]+=(double)f/time_tabl[1];
			printf("Unroll_16: MX_SIZE=%d\t time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );


			for (i=0;i<n*n;++i) 
				if (fabs(c[i]-d[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}
						

	///////////////////////////////////// blocked
					init_time();
					dgemm_blocked(n,a,b,e);
					read_time(time_tabl);
					suma_b4[count++]+=(double)f/time_tabl[1];
			printf("Blocked: MX_SIZE=%d \t time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)f/time_tabl[1]/1.0e9 );
					

			
					// sprawdzenie czy oba algorytmy daly ten sam wynik
					for (i=0;i<n*n;++i) 
						if (fabs(d[i]-e[i])>1.0e-9) {printf("Error u16 naive !\n"); goto rtrn;}



			for (i=0;i<n*n;++i) { c[i] = 0.0; d[i] = 0.0;e[i]=0.0; } 
			}
}


	rtrn:
showresult();
	return(0);
}



void showresult(){
int matrixsize=64;
for(int i=0;i<4;i++){
	printf("Matrixsize: %dx%d",matrixsize,matrixsize*=2);
printf("Srednia GFPS naive = %.3lf\n", suma_n[i]/5.0e9);
printf("Srednia GFPS UNROLL16 = %.3lf\n", suma_u16[i]/5.0e9 );
printf("Srednia GFPS BLOCKED blok4 = %.3lf\n", suma_b4[i]/5.0e9 );
}
}



void dgemm_unroll16(int n, double* A, double* B, double* C)
{
register int i,j,k;
register double reg0, reg1, reg2,  reg3,  reg4,  reg5,  reg6,  reg7, reg8, reg9, reg10, reg11, reg12, reg13, reg14, reg15;

for(i=0;i<n;i+=4)
    for(j=0;j<n;j+=4)
    {
	reg0=reg1=reg2=reg3=reg4=reg5=reg6=reg7=reg8=reg9=reg10=reg11=reg12=reg13=reg14=reg15=0.0;
	
	for(k=0;k<n;++k)
	{
	    reg0 +=A[i +k*n]*B[k+j*n];
	    reg1 +=A[i+1+k*n]*B[k+j*n];

	    reg2 +=A[i+2+k*n]*B[k+j*n];
	    reg3 +=A[i+3+k*n]*B[k+j*n];
		
	    reg4 +=A[i  +k*n]*B[k+(j+1)*n];
	    reg5 +=A[i+1+k*n]*B[k+(j+1)*n];

	    reg6 +=A[i+2+k*n]*B[k+(j+1)*n];
	    reg7 +=A[i+3+k*n]*B[k+(j+1)*n];

	    reg8 +=A[i  +k*n]*B[k+(j+2)*n];
	    reg9 +=A[i+1+k*n]*B[k+(j+2)*n];
	    reg10+=A[i+2+k*n]*B[k+(j+2)*n];
	    reg11+=A[i+3+k*n]*B[k+(j+2)*n];

	    reg12+=A[i  +k*n]*B[k+(j+3)*n];
	    reg13+=A[i+1+k*n]*B[k+(j+3)*n];
	    reg14+=A[i+2+k*n]*B[k+(j+3)*n];
	    reg15+=A[i+3+k*n]*B[k+(j+3)*n];

	}
	
	C[(i+0)+(j+0)*n] +=reg0;
	C[(i+1)+(j+0)*n] +=reg1;
	C[(i+2)+(j+0)*n] +=reg2;
	C[(i+3)+(j+0)*n] +=reg3;

	C[(i+0)+(j+1)*n] +=reg4;
	C[(i+1)+(j+1)*n] +=reg5;
	C[(i+2)+(j+1)*n] +=reg6;
	C[(i+3)+(j+1)*n] +=reg7;

	C[(i+0)+(j+2)*n] +=reg8;
	C[(i+1)+(j+2)*n] +=reg9;
	C[(i+2)+(j+2)*n] +=reg10;
	C[(i+3)+(j+2)*n] +=reg11;

	C[(i+0)+(j+3)*n] +=reg12;
	C[(i+1)+(j+3)*n] +=reg13;
	C[(i+2)+(j+3)*n] +=reg14;
	C[(i+3)+(j+3)*n] +=reg15;
	
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
