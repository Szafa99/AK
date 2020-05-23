

#include <stdio.h>
#include "eval_time.h"


void main(void)
{
	 char napis[10];
	double y, time_tabl[3];


	init_time();
	printf( "Wpisz cos  " ); 
    scanf("%s", &napis );
 	read_time( time_tabl );
	printf("\nBlock I/O     T0 = %lf  T1= %lf  T2 = %lf\n", time_tabl[0], time_tabl[1], time_tabl[2] );


		init_time();
		for(int i = 0; i < 10000; i++ )
			printf("Tekst\n");
	 	read_time( time_tabl );
		printf("I/O: T0 = %lf  T1= %lf  T2 = %lf\n", time_tabl[0], time_tabl[1], time_tabl[2] );


	init_time();
	for(unsigned long int i = 0; i < 1000000000; i++ )
		y = i / ( i + 1 );
	read_time( time_tabl );
	printf("Dzielenie:  T0 = %lf  T1= %lf  T2 = %lf\n", time_tabl[0], time_tabl[1], time_tabl[2] );

	return ;
}
