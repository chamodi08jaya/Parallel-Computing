/* mm_threads.c -  Matrix Multiplication using Threads*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>
#include <pthread.h>
#include <assert.h>

#define RANDLIMIT	5	/* Magnitude limit of generated randno.*/
#define N		1000 	/* Matrix Size is configurable -
				   Should be a multiple of THREADS */
#define THREADS		4	/* Define Number of Threads */
#define NUMLIMIT 	70.0

#ifdef BLOCKED
#define NB 50 /*No. of blocks across matrix and also down matrix in one dimension*/
#define NEIB N/NB /*No. of elements in a direction*/
#define BLOCKSIZE NEIB*NEIB /*No. of elements in a block*/
#endif

void *slave (void *myid);

/*Shared Data*/
float a[N][N];
float b[N][N];
float c[N][N];
pthread_mutex_t mutex;

void *slave( void *myid )
{
	long x, low, high;

	/*Calculate which rows to calculate by each Thread*/

	x = N/THREADS;
	low = (long) myid * x;
	high = low + x;

    	int i, j, k;

	#ifdef BLOCKED
		int p,q,r;
	#endif

    	/*Calculation*/
    	/*The locking and unlocking of the mutex varaible is not needed
    	because the each Thread is accessing a different part of the array.
	Therefore the data integrity issue is not involved.*/

		/*Add a blocked to the thread*/
			
			#ifdef BLOCKED
			assert(N%NB==0);
			for (p=0; p<NB; p++){
				for (q=0; q<NB; q++) {
				  for (r=0; r<NB; r++) {
				for (i=p*NEIB;i<p*NEIB+NEIB;i++)
				  for (j=q*NEIB;j<q*NEIB+NEIB;j++)
					for (k=r*NEIB;k<r*NEIB+NEIB;k++) {
					  c[i][j]=c[i][j] + a[i][k]*b[k][j]; /*standard way*/

					  /*c[j][i]=c[j][i] + a[j][k]*b[k][i];*/ /*makes c col by col*/
					}
				 }
				}
			}
		#else /*not blocked*/
			for (i=low; i<high; i++) {
				for (j=0; j<N; j++) {
					c[i][j] = 0.0;
					for (k=0; k<N; k++){
						//pthread_mutex_lock (&mutex);
							c[i][j] = c[i][j] + a[i][k]*b[k][j];
						//pthread_mutex_unlock (&mutex);
				}
				}
			}
		#endif

   	//printf ("\nI am thread : %d, low = %d, high =%d\n",(int) myid, low, high);

}


int main(int argc, char *argv[])
{

    	struct timeval start, stop;
    	long i,j;
    	pthread_t tid[THREADS];

		#ifdef BLOCKED
			int p,q,r;
		#endif

	assert(N%THREADS==0);

    	pthread_mutex_init (&mutex, NULL); /*initialize mutex*/

    	/* generate mxs randomly */
    	for (i=0; i<N; i++)
   		for (j=0; j<N; j++) {
        		a[i][j] = 1+(int) (NUMLIMIT*rand()/(RAND_MAX+1.0));
               		b[i][j] = (double) (rand() % RANDLIMIT);
     	 	}

#ifdef PRINT

    	/* print matrices A and B */
    	printf("\nMatrix A:\n");
    	for (i=0; i<N; i++){
    		for (j=0; j<N; j++)
        		printf("%.3f\t",a[i][j]);
       		printf("\n");
    	}

    	printf("\nMatrix B:\n");
    	for (i=0; i<N; i++){
    		for (j=0; j<N; j++)
        		printf("%.3f\t",b[i][j]);
        	printf("\n");
    	}


#endif

     	/*Start Timing*/
     	gettimeofday(&start, 0);

     	/*Create Threads*/
     	for ( i=0; i<THREADS ; i++)
		if (pthread_create( &tid[i], NULL, slave, (void *) i) != 0)
			perror ("Pthread create fails");

     	/*Join Threads*/
     	for ( i=0; i<THREADS ; i++)
		if (pthread_join( tid[i], NULL) != 0 )
		    	perror ("Pthread join fails");


    	/*End Timing*/
    	gettimeofday(&stop, 0);


#ifdef PRINT


	/*print results*/
   	printf("\nAnswer = \n");
   	for (i=0; i<N; i++){
       		for (j=0; j<N; j++)
          		printf("%.3f\t",c[i][j]);
       		printf("\n");
   	}


#endif

   	/*Print the timing details*/

    	fprintf(stdout,"Time = %.6f\n\n",
		(stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));


    	return(0);
}
