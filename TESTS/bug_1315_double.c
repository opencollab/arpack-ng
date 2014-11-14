#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/* test program to solve for the 9 largest eigenvalues of
 * A*x = lambda*x where A is the diagonal matrix
 * with entries 1000, 999, ... , 2, 1 on the diagonal.
 * We're using the non symmetric routines dnaupd and dneupd.
 * This is not efficient since the problem is
 * symmetric but is done to exhibit the bug.
 * */

extern void dnaupd_(int *, char *, int *, char *, int *,
		   double *, double *, int *, double *,
		   int *, int *, int *, double *,
		   double *, int *, int *);

extern void dneupd_( int*, char*, int *, double *, double *, double *, int*, double *,
		double *, double *, char *, int *, char *, int *, double *, double *, int *,
		double *, int *, int *, int *, double *, double *, int *, int * );

void matVec(double * x, double * y) {
  int i;
  for ( i = 0; i < 1000; ++i)
    y[i] = ((double) (i+1))*x[i];
};

int main() {
  int ido = 0;
  char bmat[] = "I";
  int N = 1000;
  char which[] = "LM";
  int nev = 9;
  double tol = 0;
  double resid[N];
  int ncv = 2*nev+1;
  double V[ncv*N];
  int ldv = N;
  int iparam[11];
  int ipntr[14];
  double workd[3*N];
  int rvec = 1;
  char howmny[] = "A";
  double* dr = (double*) malloc((nev+1)*sizeof(double));
  double* di = (double*) malloc((nev+1)*sizeof(double));
  int select[3*ncv];
  double z[(N+1)*(nev+1)];
  int ldz = N+1;
  double sigmar=0;
  double sigmai=0;
  double workev[3*ncv];
  int k;
  for (k=0; k < 3*N; ++k )
    workd[k] = 0;
  double workl[3*(ncv*ncv) + 6*ncv];
  for (k=0; k < 3*(ncv*ncv) + 6*ncv; ++k )
    workl[k] = 0;
  int lworkl = 3*(ncv*ncv) + 6*ncv;
  int info = 0;

  iparam[0] = 1;
  iparam[2] = 10*N;
  iparam[3] = 1;
  iparam[6] = 1;

  dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	 workd, workl, &lworkl, &info);

  while(ido == 1) {

    matVec(&(workd[ipntr[0]-1]), &(workd[ipntr[1]-1]));

    dnaupd_(&ido, bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	 workd, workl, &lworkl, &info);
  }

  dneupd_( &rvec, howmny, select, dr,di, z, &ldz, &sigmar, &sigmai,workev,
	   bmat, &N, which, &nev, &tol, resid, &ncv, V, &ldv, iparam, ipntr,
	   workd, workl, &lworkl, &info);
  int i;
  for (i = 0; i < nev; ++i) {
    printf("%f\n", dr[i]);
    if(fabs(dr[i] - (double)(1000-i))>1e-6){
      exit(EXIT_FAILURE);
    }
  }
  return 0;
}
