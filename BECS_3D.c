//Solve 3D BECS
//Output z slice of the Temp n+1 vector
//command line argument tells you whether it's Jacobi, GaussSeidel, SOR method

#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include "nrutil.h"
#include <time.h>

void mprint(double **matrix, int m, char *label);
void vprint(double *vector, int m, char *label);
void free_dmatrix(double **m, long nrl, long nrh, long ncl, long nch);
void free_dvector(double *v, long nl, long nh);
void BECS(int nx, int ny, int nz, double Lx, double Ly, double Lz, int tsteps, double dt, double alpha, char* method);
void Jacobi(double **A, double *x, double *b, int n);
void GaussSeidel(double **A, double *x, double *b, int n);
void SOR(double **A, double *x, double *b, int n);

int main(int argc, char** argv)
{
  int nx, ny, nz;
  double Lx, Ly, Lz; //size of x,y,z dimensions
  Lx = 1;
  Ly = 1;
  Lz = 1;

  int tsteps = 100; //number of time steps
  double dt = 0.0005; //size of each time step
  double alpha = 0.001; //diffusivity constant
  time_t start, end;

  for(nx = 1; nx <= 25; nx++)
    {
      time(&start);

      //temporarily set nx = ny = nz for timing purposes
      ny = nx;
      nz = nx;

      //check valid command line argument
      if((strcmp(argv[1],"SOR")!=0)&&(strcmp(argv[1],"Jacobi")!=0)&&(strcmp(argv[1],"GaussSeidel")!=0))
	{
	  printf("Please input solver method: Jacobi, GaussSeidel, SOR\n");
	  return;
	}

      //run BECS3D on size nx x elements, cubic domain
      BECS(nx, ny, nz, Lx, Ly, Lz, tsteps, dt, alpha, argv[1]);
      time(&end);
      printf("%d: %f\n", nx, difftime(end, start));
    }
}

void BECS(int nx, int ny, int nz, double Lx, double Ly, double Lz, int tsteps, double dt, double alpha, char* method)
{
  double dx, dy, dz; //width of each slice in x,y,z dim
  dx = Lx/nx;
  dy = Ly/ny;
  dz = Lz/nz;

  double Cx, Cy, Cz; //constants, stable if Cx + Cy + Cz <= 0.5
  Cx = alpha*dt/dx/dx;
  Cy = alpha*dt/dy/dy;
  Cz = alpha*dt/dz/dz;

  //Bu = b. Given B and b, solve for u.
  double **B;
  double **C; //C so that B doesn't get changed during Gauss Elim
  double *u;
  double *b;

  B = dmatrix(1,nx*ny*nz+1,1,nx*ny*nz+1);
  C = dmatrix(1,nx*ny*nz+1,1,nx*ny*nz+1);
  b = dvector(1,nx*ny*nz+1);
  u = dvector(1,nx*ny*nz+1);

  //Populate coefficient matrix B with 7 diagonals (3 centered, 4 on sides)
  double diagonal = 1+2*Cx+2*Cy+2*Cz;
  long int x, y, z;
  for(x = 1; x <= nx*ny*nz; x++)
    {
      for(y = 1; y <= nx*ny*nz; y++)
	{
	  if(x == y)
	    {
	      B[x][y] = diagonal;
	    }
	  else
	    {
	      if(abs(x-y) == 1) //adjacent to diagonal diagonal
		{
		  B[x][y] = -Cx;
		}
	      else if(abs(x-y) == nx) //nx away from diagonal on left or right
		{
		  B[x][y] = -Cy;
		}
	      else if(abs(x-y) == nx*ny) //nx*ny away from diagonal on left or right
		{
		  B[x][y] = -Cz;
		}
	      else B[x][y] = 0;
	    }
	}
    }

  long int i,j,k;
  double dirichlet = 0;
  double tmpx, tmpy, tmpz;

  //Initialize b, T at time 0, Gaussian plus noise
  for(i = 1; i <= nx*ny*nz; i++)
    {
      //calculate x, y, z coordinates
      x = (i % (nx*ny)) % nx;
      y = (i % (nx*ny)) / nx + 1;
      z = i / (nx*ny) + 1;

      //begin dirichlet boundary
      if((x == 0) || (x == nx))
	{
	  b[i] = dirichlet;
	}
      else if((y == 0) || (y == ny))
	{
	  b[i] = dirichlet;
	}
      else if((z == 0) || (z == nz))
	{
	  b[i] = dirichlet;
	}
      else
	{
	  tmpx = 5*x*dx-2.5;
	  tmpy = 5*y*dy-2.5;
	  tmpz = 5*z*dz-2.5;
	  b[i] = exp(-(tmpx*tmpx))*exp(-(tmpz*tmpz))*exp(-(tmpz*tmpz));
	}
      //end dirichlet boundary
    }

  //Find u at each time step from T = 0 to T = dt*tsteps
  int t, counter;;
  double rowsum;

  for(t = 1; t <= tsteps; t++)
    {
      //time elapsed t*dt

      //create copy of B in C for iterative solver
      for(i = 1; i <= nx*ny*nz; i++)
	{
	  for(j = 1; j <= nx*ny*nz; j++)
	    {
	      C[i][j] = B[i][j];
	    }
	}

      //solver in command line
      if(strcmp(method, "Jacobi")==0)
	{
	  Jacobi(C,u,b,nx*ny*nz);
	}
      else if(strcmp(method, "GaussSeidel")==0)
	{
	  GaussSeidel(C,u,b,nx*ny*nz);
	}
      else if(strcmp(method, "SOR")==0)
	{
	  SOR(C,u,b,nx*ny*nz);
	}

      //set b =u;
      for(k = 1; k <= nx*ny*nz; k++)
	{
	  b[k] = u[k];
	}
    }
  free_dmatrix(B, 1, nx*ny*nz+1, 1, nx*ny*nz+1);
  free_dmatrix(C, 1, nx*ny*nz+1, 1, nx*ny*nz+1);
  free_dvector(u, 1, nx*ny*nz+1);
  free_dvector(b, 1, nx*ny*nz+1);
}

/* print a matrix in not too horrible a way */
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

/* print a vector in not too horrible a way */
void vprint(double *vector, int m, char *label){
  int i, j;
  printf("%s:\n",label);

  for (i = 1; i <= m; ++i){
    printf("%10.2f ", vector[i]);
  }

  printf("\n------------------------\n");
}
