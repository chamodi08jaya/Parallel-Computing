 /* mm.c - blocked vs standard matrix Multiplication */


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/time.h>
#include <assert.h>

#define RANDLIMIT	5	/* Magnitude limit of generated randno.*/
#define N	2000 	/* Matrix Size */
#define NUMLIMIT 70.0

#ifdef BLOCKED
#define NB 100 /* no.of blocks across matrix and also down matrix */ 
#define NEIB N/NB /* no. of elements in a direction */
#define BLOCKSIZE NEIB*NEIB /* no. of elements in a block */
#endif

float a[N][N];
float b[N][N];
float c[N][N];


int main(int argc, char *argv[]) {
  struct timeval start, stop;	
  int i,j,k;

#ifdef BLOCKED
  int p,q,r;
#endif
  

  /* generate mxs  */

  for (i=0; i<N; i++)
    for (j=0; j<N; j++) {
      a[i][j] = 1+(int) (NUMLIMIT*rand()/(RAND_MAX+1.0));
      /* a[i][j] = 1.0; 
	 b[i][j] = 2.0;*/
      b[i][j] = (double) (rand() % RANDLIMIT);
      /*c[i][j] = 0.0;*/
    }

#ifdef PRINT
  /* print matrices */
  printf("Matrix A:\n");           
  for (i=0; i<N; i++){      
    for (j=0; j<N; j++)  
      printf("%.3f\t",a[i][j]); 
    printf("\n");  
  }        

  printf("Matrix B:\n");       
  for (i=0; i<N; i++){  
    for (j=0; j<N; j++)
      printf("%.3f\t",b[i][j]);
    printf("\n");
  }   

  printf("Matrix C:\n");       
  for (i=0; i<N; i++){  
    for (j=0; j<N; j++)
      printf("%.3f\t",c[i][j]);
    printf("\n");
  }   
#endif

  gettimeofday(&start, 0);

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
#else /* not blocked */
  for (i=0; i<N; i++) {
    for (j=0; j<N; j++) {
      for (k=0; k<N; k++) {
	c[i][j] = c[i][j] + a[i][k]*b[k][j]; /* Working;standard way */
      /*c[j][i] = c[j][i] + a[j][k]*b[k][i];*/ /* Working; Makes C column by col */
      }
    } /* end j loop */
  }
#endif /* not blocked end */

  gettimeofday(&stop, 0);

#ifdef PRINT
  /* print results*/
  printf("Answer c:\n");
  for (i=0; i<N; i++){
    for (j=0; j<N; j++) 
      printf("%.3f\t",c[i][j]);
    printf("\n");
  }
#endif

  fprintf(stdout,"Time = %.6f\n\n",
	  (stop.tv_sec+stop.tv_usec*1e-6)-(start.tv_sec+start.tv_usec*1e-6));
  return(0);
}
