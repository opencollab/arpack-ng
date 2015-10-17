#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * We're using the non symmetric routines dnaupd and dneupd.
 * This is not efficient since the problem is
 * symmetric but is done to exhibit the bug.
 */

extern void snaupd_(int *, char *, int *, char *, int *,
		   float *, float *, int *, float *,
		   int *, int *, int *, float *,
		   float *, int *, int *);

extern void sneupd_( int*, char*, int *, float *, float *, float *, int*, float *,
		float *, float *, char *, int *, char *, int *, float *, float *, int *,
		float *, int *, int *, int *, float *, float *, int *, int * );


void matVec(float * x, float * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((float) (i+1))*x[i];
};

int main() {
  int ido = 0;
  char bmat[] = "I";
  int N = 1000;
  char which[] = "LM";
  int nev = 9;
  float tol = 0;
  float resid[N];
  int ncv = 2*nev+1;
  float V[ncv*N];
  int ldv = N;
  int iparam[11];
  int ipntr[14];
  float workd[3*N];
  int rvec = 1;
  char howmny[] = "A";
  float* dr = (float*) malloc((nev+1)*sizeof(float));
  float* di = (float*) malloc((nev+1)*sizeof(float));
  int select[3*ncv];
  float z[(N+1)*(nev+1)];
  int ldz = N+1;
  float sigmar=0;
  float sigmai=0;
  float workev[3*ncv];
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0;
  float workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0;
  int lworkl = 3*(ncv*ncv) + 6*ncv;
  int info = 0;

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[6] = 1;

  snaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	 workd, workl, &lworkl, &info);

  while(ido == 1) {

    matVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));

    snaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	 workd, workl, &lworkl, &info);
  }

  sneupd_( &rvec, howmny, select, dr,di, z, &ldz, &sigmar, &sigmai,workev,
	   bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	   workd, workl, &lworkl, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f\n", dr[i]);
    if(fabs(dr[i] - (float)(1000-i))>1e-2){
      free(dr);
      free(di);
      exit(EXIT_FAILURE);
    }
  }
  free(dr);
  free(di);
  return 0;
}
