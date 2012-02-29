//use with BECS_3D.c to solve using Jacobi method
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"

void Jacobi(double **A, double *x, double *b, int n);
void mprint(double **matrix, int m, char *label);
void vprint(double *vector, int m, char *label);

/*
  //comment out main when used with BECS_3D
  int main(int argc, char** argv)
  {
    //driver to check Jacobi 3x3 matrix
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
    Jacobi(A,x,b,n);
    vprint(x,n,"x new");
  }
*/
void Jacobi(double **A, double *x, double *b, int n)
{
  //Solve equation Ax = b using Jacobi method
  //A is of size n by n. x and b are vectors

  int MAX_ITER = 1000; //maximum iterations
  double MAX_ERROR = 0.0001; //permissible error size
  int i,j,k;

  double *Tx; //product T and xold, n by 1 vector
  double *C; //Dinv*b
  double **Dinv; //Dinv*R
  double **negDinv; //negative Dinv
  double **T; //-Dinv*R
  double **R; //non-diagonal part of A
  double *xnew;
  double *Ax; //A*x

  Tx = dvector(1,n);
  C = dvector(1,n);
  Dinv = dmatrix(1,n,1,n);
  negDinv = dmatrix(1,n,1,n);
  T = dmatrix(1,n,1,n);
  R = dmatrix(1,n,1,n);
  xnew = dvector(1,n);
  Ax = dvector(1,n);

  //find Dinv
  for(i = 1;i<=n;i++)
    {
      Dinv[i][i] = 1/A[i][i];
      negDinv[i][i] = -Dinv[i][i];
      for(j = 1;j <=n;j++)
	{
	  if(i!=j)
	    {
	      R[i][j] = A[i][j];
	    }
	}
    }

  VectorMultiplication(Dinv, b, C, n); //find C
  MatrixMultiplication(negDinv,R, T, n); //find T

  //only need T, Tx and C and Ax
  free_dmatrix(Dinv, 1,n,1,n);
  free_dmatrix(negDinv, 1,n,1,n);
  free_dmatrix(R, 1,n,1,n);

  //initialize x = 1
  for(i = 1;i<=n;i++)
    {
      x[i] = 1;
    }

  int iteration;
  double sumerror;
  for(iteration = 1;iteration < MAX_ITER; iteration++)
    {
      VectorMultiplication(T,x,Tx,n);
      for(i = 1; i<=n;i++)
	{
	  xnew[i] = Tx[i] + C[i];
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

  free_dvector(Tx,1,n);
  free_dvector(C,1,n);
  free_dmatrix(T,1,n,1,n);
  free_dvector(Ax,1,n);
  free_dvector(xnew,1,n);
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
