#include <stdio.h>
#include <math.h>
#include "eval_time.h"
#include <x86intrin.h>

#define SIZE 2048			


static double a[SIZE*SIZE];
static double b[SIZE*SIZE];
static double c[SIZE*SIZE];
static double d[SIZE*SIZE];
static double e[SIZE*SIZE];
static double f[SIZE*SIZE];
static double g[SIZE*SIZE];
static double h[SIZE*SIZE];

void naive_avx(int n, double* A, double* B, double* C);
void unroll_avx(int,double*,double*,double*);
void unroll4_avx(int,double*,double*,double*);
void unroll8_avx(int n, double* A, double* B, double* C);

void dgemm_blocked_sse(int n,int blocksize, double* A, double* B, double* C,void(*blockptr)(int n,int blocksize, int bi, int bj, int bk, double *A, double *B, double *C));
void block_sse(int,int,int,int,int,double*,double*,double*);
void unroll_block_sse(int,int,int,int,int,double*,double*,double*);


void showresult();

static double suma_n[4];
static double suma_u[4];
static double suma_u4[4];
static double suma_u8[4];
static double suma_bl[4][5];
static double suma_ubl[4][5];

  int blocksize;


int main(void)
{
	long int i,j,n,ff;
	double time_tabl[3];
	int count=0;

for(int z=0;z<5;z++){
	count=0;	
	for(int size = 512; size <= 2048; size *=2 ) {
		
		n = size;
		ff = 2*n*n*n;		

			for(i=0;i<n;++i)	
				for (j=0;j<n;++j)
				{
				a[j+i*n]=(double)(i+j*n);
				b[j+i*n]=(double)(j+i*n);
				}
			//c i d zostaw wyzerowane

	////////////////////////////////////////// Naive
			init_time();
			naive_avx(n,a,b,c);
			read_time(time_tabl);
			suma_n[count]+=(double)ff/time_tabl[1];
			printf("Naive:MX_SIZE=%d \ttime = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );


	//////////////////// UNROLL
			init_time();
			unroll_avx(n,a,b,d);
			read_time(time_tabl);
			suma_u[count]+=(double)ff/time_tabl[1];
			printf("Unroll: MX_SIZE=%d\t time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );


			for (i=0;i<n*n;++i) 
				if (fabs(c[i]-d[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}

	//////////////////// UNROLL_4
			init_time();
			unroll4_avx(n,a,b,e);
			read_time(time_tabl);
			suma_u4[count]+=(double)ff/time_tabl[1];
			printf("Unroll4: MX_SIZE=%d\t time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );


			for (i=0;i<n*n;++i) 
				if (fabs(d[i]-e[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}

							//////////////////// UNROLL_8
			init_time();
			unroll8_avx(n,a,b,f);
			read_time(time_tabl);
			suma_u8[count]+=(double)ff/time_tabl[1];
			printf("Unroll8: MX_SIZE=%d\t time = %2.3lf GFLOPS = %.3lf\n", n, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );


			for (i=0;i<n*n;++i) 
				if (fabs(e[i]-f[i])>1.0e-9) {printf("Error!\n"); goto rtrn;}
						

	///////////////////////////////////// blocked
	int countbl=0;
	for( int blocksize=4;blocksize<=128;  blocksize*=2){
					init_time();
					dgemm_blocked_sse(n,blocksize, a,b,g,block_sse);
					read_time(time_tabl);
					suma_bl[count][countbl++]+=(double)ff/time_tabl[1];
			printf("Blocked: MX_SIZE:%d \t BLOCK_SIZE:%d time = %2.3lf GFLOPS = %.3lf\n", n, blocksize, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );
					// sprawdzenie czy oba algorytmy daly ten sam wynik
					for (i=0;i<n*n;++i) 
						if (fabs(f[i]-g[i])>1.0e-9) {printf("Error !\n"); goto rtrn;}
					
			for (i=0;i<n*n;++i)g[i] = 0.0; 
	
	}
	
	//////////////////////////////////////// unroll_blocked
			 countbl=0;
				for( int blocksize=4;blocksize<=128;  blocksize*=2){
					init_time();
					dgemm_blocked_sse(n,blocksize, a,b,h,unroll_block_sse);
					read_time(time_tabl);
					suma_ubl[count][countbl++]+=(double)ff/time_tabl[1]/1.0e9;
			printf("Blocked: MX_SIZE:%d \t BLOCK_SIZE:%d time = %2.3lf GFLOPS = %.3lf\n", n, blocksize, time_tabl[1],(double)ff/time_tabl[1]/1.0e9 );
					//sprawdzenie czy oba algorytmy daly ten sam wynik
					for (i=0;i<n*n;++i) 
						if (fabs(f[i]-g[i])>1.0e-9) {printf("Error !\n"); goto rtrn;}
					
			if(blocksize==128)count++;
			for (i=0;i<n*n;++i)h[i] = 0.0; 
	
	}




			for (i=0;i<n*n;++i) { c[i] = 0.0; d[i] = 0.0;e[i]=0.0;f[i]=0.0; } 
			}
}


showresult();
	rtrn:
	return(0);
}




void naive_avx(int n, double* A, double* B, double* C)
{
__m256d cv;

for(register int i=0;i<n;i+=4)
    for(register int j=0;j<n;++j)
    {

	cv=_mm256_load_pd(C+i+j*n); 	
	for(register int k=0;k<n;++k){
	    cv=_mm256_add_pd(cv, _mm256_mul_pd(_mm256_loadu_pd(A+i+k*n),_mm256_broadcast_sd(B+k+j*n))); 	
	
	}
	_mm256_store_pd(C+i+j*n,cv); 			
    }
}

void unroll_avx(int n, double* A, double* B, double* C)
{
register int i,j,k;
__m256d reg0,reg1,reg2;

for(i=0;i<n;i+=8)
    for(j=0;j<n;++j)
    {
	reg0 = _mm256_load_pd(C+i+0+j*n);
    reg1 = _mm256_load_pd(C+i+4+j*n);
	for(k=0;k<n;++k)
	{
        reg2 = _mm256_broadcast_sd(B+k+j*n);
        reg0 = _mm256_add_pd(reg0,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+0+i), reg2));
        reg1 = _mm256_add_pd(reg1,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+4+i), reg2));
	}
	_mm256_store_pd(C+i+0+j*n, reg0);
    _mm256_store_pd(C+i+4+j*n, reg1);
    }
}



void unroll4_avx(int n, double* A, double* B, double* C)
{
register int i,j,k;
__m256d reg0,reg1,reg2,reg3,reg4;

for(i=0;i<n;i+=16)
    for(j=0;j<n;++j)
    {
	reg0 = _mm256_load_pd(C+i+0+j*n);
    reg1 = _mm256_load_pd(C+i+4+j*n);
    reg2 = _mm256_load_pd(C+i+8+j*n);
    reg3 = _mm256_load_pd(C+i+12+j*n);
	for(k=0;k<n;++k)
	{
        reg4 = _mm256_broadcast_sd(B+k+j*n);

        reg0 = _mm256_add_pd(reg0,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+0+i), reg4));
        reg1 = _mm256_add_pd(reg1,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+4+i), reg4));
        reg2 = _mm256_add_pd(reg2,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+8+i), reg4));
        reg3 = _mm256_add_pd(reg3,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+12+i), reg4));
	}
	_mm256_store_pd(C+i+0+j*n, reg0);
    _mm256_store_pd(C+i+4+j*n, reg1);
    _mm256_store_pd(C+i+8+j*n, reg2);
    _mm256_store_pd(C+i+12+j*n, reg3);
    }
}

void unroll8_avx(int n, double* A, double* B, double* C)
{
register int i,j,k;
__m256d reg0,reg1,reg2,reg3,reg4,reg5,reg6,reg7,reg8;

for(i=0;i<n;i+=32)
    for(j=0;j<n;++j)
    {
	reg0 = _mm256_load_pd(C+i+0+j*n);
    reg1 = _mm256_load_pd(C+i+4+j*n);
    reg2 = _mm256_load_pd(C+i+8+j*n);
    reg3 = _mm256_load_pd(C+i+12+j*n);
	reg4 = _mm256_load_pd(C+i+16+j*n);
    reg5 = _mm256_load_pd(C+i+20+j*n);
    reg6 = _mm256_load_pd(C+i+24+j*n);
    reg7 = _mm256_load_pd(C+i+28+j*n);
	for(k=0;k<n;++k)
	{
        reg8 = _mm256_broadcast_sd(B+k+j*n);

        reg0 = _mm256_add_pd(reg0,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+0+i), reg8));
        reg1 = _mm256_add_pd(reg1,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+4+i), reg8));
        reg2 = _mm256_add_pd(reg2,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+8+i), reg8));
        reg3 = _mm256_add_pd(reg3,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+12+i), reg8));
		reg4 = _mm256_add_pd(reg4,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+16+i), reg8));
        reg5 = _mm256_add_pd(reg5,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+20+i), reg8));
        reg6 = _mm256_add_pd(reg6,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+24+i), reg8));
        reg7 = _mm256_add_pd(reg7,_mm256_mul_pd(_mm256_loadu_pd(A+n*k+28+i), reg8));
	}
	_mm256_store_pd(C+i+0+j*n, reg0);
    _mm256_store_pd(C+i+4+j*n, reg1);
    _mm256_store_pd(C+i+8+j*n, reg2);
    _mm256_store_pd(C+i+12+j*n, reg3);
	_mm256_store_pd(C+i+16+j*n, reg4);
    _mm256_store_pd(C+i+20+j*n, reg5);
    _mm256_store_pd(C+i+24+j*n, reg6);
    _mm256_store_pd(C+i+28+j*n, reg7);
	}
}





void dgemm_blocked_sse(int n,int blocksize, double* A, double* B, double* C,void(*blockptr)(int ,int , int , int , int , double* , double* , double*) )
{
register int bi,bj,bk;

for(bi=0;bi<n;bi+=blocksize)
    for(bj=0;bj<n;bj+=blocksize)
	for(bk=0;bk<n;bk+=blocksize)
	    blockptr(n ,blocksize, bi,bj,bk,A,B,C);
}

void block_sse(int n,int blocksize, int bi, int bj, int bk, double *A, double *B, double *C)
{

__m128d cij;
register int i,j,k;


for(i=bi;i<bi+blocksize;i+=2)
    for(j=bj;j<bj+blocksize;++j)
    {
	cij= _mm_load_pd(C+i+j*n);
	
	for(k=bk;k<bk+blocksize;++k)
	    cij=_mm_add_pd(cij, _mm_mul_pd(_mm_loadu_pd(A+i+k*n),_mm_load1_pd(B+k+j*n))); 	

	_mm_store_pd(C+i+j*n,cij);
    }
}



void unroll_block_sse(int n,int blocksize, int bi, int bj, int bk, double *A, double *B, double *C)
{

__m128d cv,reg0,reg1;
register int i,j,k;


for(i=bi;i<bi+blocksize;i+=4)
    for(j=bj;j<bj+blocksize;++j)
    {
	reg0=_mm_load_pd(C+i+0+j*n);
	reg1=_mm_load_pd(C+i+2+j*n);

	for(k=bk;k<bk+blocksize;++k)
        cv = _mm_load1_pd(B+k+j*n);
        
		reg0 = _mm_add_pd(reg0,_mm_mul_pd(_mm_load_pd(A+n*k+0+i), cv));
        reg1 = _mm_add_pd(reg1,_mm_mul_pd(_mm_load_pd(A+n*k+2+i), cv));

	_mm_store_pd(C+i+0+j*n,cv);
	_mm_store_pd(C+i+2+j*n,cv);
    }
}







void showresult(){
int matrixsize=256;
printf("\n\n---------WYNIK KONCOWY------------\n\n");
for(int i=0;i<3;i++){
	printf("\nMatrixsize: %dx%d\n\n",matrixsize,matrixsize*=2);

printf("Srednia GFPS naive = %.3lf\n", suma_n[i]/5.0e9);
printf("Srednia GFPS UNROLL = %.3lf\n", suma_u[i]/5.0e9 );
printf("Srednia GFPS UNROLL4 = %.3lf\n", suma_u4[i]/5.0e9 );
printf("Srednia GFPS UNROLL8 = %.3lf\n", suma_u8[i]/5.0e9 );
	int bl=2;
for(int b=0;b<5;b++){
	printf("Srednia GFPS BLOCKED blok:%d = %.3lf\n",bl, suma_bl[i][b]/5.0e9 );
	printf("Srednia GFPS unroll_BLOCKED blok:%d = %.3lf\n",bl*=2, suma_ubl[i][b]/5 );
	}
}
}
