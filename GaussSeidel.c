//use with BECS_3D.c to solve using Gauss Seidel method
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"

void GaussSeidel(double **A, double *x, double *b, int n);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void mprint(double **matrix, int m, char *label);
void vprint(double *vector, int m, char *label);
/*
//comment out main when used with BECS_3D
int main(int argc, char** argv)
{
  //driver to check GS 3x3 matrix
  double **A, *b, *x;
  int n = 3;
  A = dmatrix(1,n,1,n);
  b = dvector(1,n);
  x = dvector(1,n);

  A[1][1] = 9; A[1][2]= 2; A[1][3] = 2;
  A[2][1] = 5; A[2][2]= 4; A[2][3] = -3;
  A[3][1] = 3; A[3][2]= 3; A[3][3] = 9;

  b[1] =9; b[2] = 1; b[3] = 0;

  mprint(A,n,"A original");
  vprint(b,n,"b original");
  vprint(x,n,"x old");

  GaussSeidel(A,x,b,n);

  vprint(x,n,"x new");
}
*/

void GaussSeidel(double **A, double *x, double *b, int n)
{
  //Solve equation Ax = b using GS method
  //A is of size n by n. x and b are vectors

  int MAX_ITER = 1000; //maximum iterations
  double MAX_ERROR = 0.0001; //permissible error size
  int i,j,k;

  double *xnew;
  double *bguess; //A*xnew, should approximate b
  xnew = dvector(1,n);
  bguess = dvector(1,n);

  //initialize x and xnew
  for(i = 1;i<=n;i++)
    {
      x[i] = 1;
      xnew[i] = 1;
    }

  int iteration;
  double sumerror;
  double sumax; //summation A[i][j]*x[j], j>i
  double sumaxnew; //summation A[i][j]*xnew[j], j<i
  for(iteration = 1;iteration < MAX_ITER; iteration++)
    {
      for(i = 1; i<=n;i++)
	{
	  sumax = 0;
	  sumaxnew = 0;
	  for(j = 1; j<i; j++)
	    {
	      sumaxnew+=A[i][j]*xnew[j];
	    }
	  for(j = i+1; j<=n;j++)
	    {
	      sumax+=A[i][j]*x[j];
	    }

	  xnew[i] = (b[i] - sumax-sumaxnew)/A[i][i];
	}
      VectorMultiplication(A,xnew,bguess,n);
      sumerror = 0;
      //stop if error reaches threshold
      for(j = 1; j<=n; j++)
	{
	  sumerror += fabs(bguess[j]-b[j]);
	}

      if(sumerror < MAX_ERROR)
	{
	  return;
	}
      //set x = xnew
      for(k = 1; k<=n;k++)
	{
	  x[k] = xnew[k];
	}
    }
  free_dvector(xnew,1,n);
  free_dvector(bguess,1,n);
}
/*
void mprint(double **matrix, int m, char *label){
  int i,j;
  printf("%s:\n", label);
  for (i=1;i<=m;++i){
    for(j=1;j<=m;++j){
      printf("%10.5f ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n------------\n");
}

void vprint(double *vector, int m, char *label){
  int i, j;
  printf("%s:\n",label);

  for (i = 1; i <= m; ++i){
    printf("%10.5f ", vector[i]);
  }
  printf("\n------------------------\n");
}
*/
