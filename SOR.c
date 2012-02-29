//use with BECS_3D.c to solve using SOR method
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"

void SOR(double **A, double *x, double *b, int n);
void mprint(double **matrix, int m, char *label);
void vprint(double *vector, int m, char *label);
/*
//comment out main and mprint vprint when used with BECS_3D
int main(int argc, char** argv)
{
  //driver to check SOR 3x3 matrix
  double **A, *b, *x;
  int n = 3;
  A = dmatrix(1,n,1,n);
  b = dvector(1,n);
  x = dvector(1,n);

  A[1][1] = 1; A[1][2]= 1; A[1][3] = 2;
  A[2][1] = 2; A[2][2]= 4; A[2][3] = -3;
  A[3][1] = 3; A[3][2]= 6; A[3][3] = -5;

  b[1] =9; b[2] = 1; b[3] = 0;

  mprint(A,n,"A original");
  vprint(b,n,"b original");

  SOR(A,x,b,n);

  vprint(x,n,"x new");
}
*/

void SOR(double **A, double *x, double *b, int n)
{
  //Solve equation Ax = b using SOR method
  //A is of size n by n. x and b are vectors

  int MAX_ITER = 1000; //maximum iterations
  double MAX_ERROR = 0.0001; //permissible error size
  double w = 1; //relaxation factor
  int i,j,k;

  double *xnew;
  double *Ax; //A*x
  double *Axnew; //A*xnew

  xnew = dvector(1,n);
  Ax = dvector(1,n);
  Axnew = dvector(1,n);

  //initialize x = 1
  for(i = 1;i<=n;i++)
    {
      x[i] = 1;
      xnew[i] = 1;
    }

  int iteration;
  double sumerror;
  double sumax, sumaxnew;
  for(iteration = 1;iteration < MAX_ITER; iteration++)
    {
      for(i = 1; i<=n;i++)
	{
	  sumax = 0;
	  sumaxnew = 0;
	  for(j = 1; j<=n;j++)
	    {
	      if(j>i) sumax += A[i][j]*x[j];
	      else if(j<i) sumaxnew += A[i][j]*xnew[j];
	    }
	  xnew[i] = (1-w)*x[i]+(w/A[i][i])*(b[i]-sumax-sumaxnew);
	}

      VectorMultiplication(A,xnew,Ax, n);
      sumerror = 0;
      //stop if error reaches threshold
      for(j = 1; j<=n; j++)
	{
	  sumerror += fabs(Ax[j]-b[j]);
	}

      if(sumerror < MAX_ERROR) return;

      //set x = xnew
      for(k = 1; k<=n;k++)
	{
	  x[k] = xnew[k];
	}
    }
  //  printf("didn't converge in %d iterations\n", MAX_ITER);

  free_dvector(xnew,1,n);
  free_dvector(Ax,1,n);
  free_dvector(Axnew,1,n);
}
/*
void mprint(double **matrix, int m, char *label){
  int i,j;
  printf("%s:\n", label);
  for (i=1;i<=m;++i){
    for(j=1;j<=m;++j){
      printf("%10.2f ", matrix[i][j]);
    }
    printf("\n");
  }
  printf("\n------------\n");
}

void vprint(double *vector, int m, char *label){
  int i, j;
  printf("%s:\n",label);

  for (i = 1; i <= m; ++i){
    printf("%10.2f ", vector[i]);
  }

  printf("\n------------------------\n");
}
*/
