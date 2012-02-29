#include <stdio.h>
#include <math.h>

void MatrixMultiplication(double **A, double **B, double **C, int n); //multiply 2 n*n matrices;
void VectorMultiplication(double **A, double *b, double *c, int n); //multiply n*n matrix with n*1 vector

void MatrixMultiplication(double **A, double **B, double **C, int n)
{
  double tempsum;
  int i,j,k;

  //A*B = C
  for(i = 1; i<=n;i++)
    {
      for(j = 1; j<=n;j++)
	{
	  tempsum = 0;
	  for(k = 1; k<=n;k++)
	    {
	      tempsum+= A[i][k]*B[k][j];
	    }
	  C[i][j] = tempsum;
	}
    }
}

void VectorMultiplication(double **A, double *b, double *c, int n)
{
  double tempsum;
  int i,j;

  //A*b = c
  for(i = 1;i<=n;i++)
    {
      tempsum = 0;
      for(j = 1;j<=n;j++)
	{
	  tempsum+= A[i][j]*b[j];
	}
      c[i] = tempsum;
    }
}
