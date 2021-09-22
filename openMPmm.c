/* A matrix multiplication code using OpenMP */

#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <omp.h>
#include <stdio.h>
#include <assert.h>

#define RANDLIMIT	5	/* Magnitude limit of generated randno.*/
#define N		2000
#define CHUNKSIZE       250
#define NUMLIMIT 	70.0

#ifdef BLOCKED
#define NB 100 /*No. of blocks across matrix and also down matrix in one dimension*/
#define NEIB N/NB /*No. of elements in a direction*/
#define BLOCKSIZE NEIB*NEIB /*No. of elements in a block*/
#endif

float a[N][N];
float b[N][N];
float c[N][N];

int main (int argc, char *argv[]) {
  struct timeval start, stop;
  long i,j,k;
  int id, chunk;
  
  #ifdef BLOCKED
	  int p,q,r;
  #endif

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
  
  chunk = CHUNKSIZE;
  
  /*Add a blocked to the openMP*/
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
	  #pragma omp parallel for shared(a,b,c,chunk) private(i,j,k,id) schedule(static,chunk)
		for (i=0; i<N; i++) {
		  for (j=0; j<N; j++) {
		c[i][j] = 0.0;
		for (k=0; k<N; k++)
		  c[i][j] = c[i][j] + a[i][k]*b[k][j];
		  }
		}
  #endif
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
