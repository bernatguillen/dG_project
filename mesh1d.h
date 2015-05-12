#ifndef MESH1D_H_
#define MESH1D_H_

#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>

int JacobiP(int N, double a, double b, double *x,double *fx, int npoints); //Jacobi polynomials degree N
int GradJacobiP(int N, double a, double b, double *x, double *fx, int npoints); //Jacobi polynomials derivative
int JacobiGQ(int N, double a, double b, double *x, double *w); //Jacobi Gauss Quadratures
int JacobiGL(int N, double a, double b, double *x); //Gauss Lobatto Quadratures 
int Vandermonde1D(int N, double *x, double **V, int npoints); //Generalized Vandermonde matrix
int GradVandermonde1D(int N, double *x, double **V, int npoints); //Derivative of Generalized Vandermonde matrix
int Dmatrix1D(int N, double *x, double **V, double **D, int npoints); //Derivative Matrix (u_h' = Du_h)



int MeshGen1D(double xmin, double xmax, int K, int *Nv, double *Vx, int **EToV); //Uniform mesh generator 1D


class dG1D_Framework{
public:
	dG1D_Framework(int K_in, int Nv_in, int N_in, double *Vx_in, int **EToV_in);
	~dG1D_Framework();
	int K,Nv,N,Np,Nfp,Nfaces;
	double NODETOL;
	double *Vx; //Vertex nodes
	int **EToV;
	double *r; //LGL points
	double **V;
	double **Vinv;
	double **x; //Physical points
	double **Dr; //Dmatrix
	double **LIFT; //LIFT Matrix
	double **rx;
	double **J; //Geometric Factors
	double **nx; //outward normals
	double **Fx; //masked array of points (only vertices of each node)
	double **Fscale; //masked inverse jacobian J
	int **Fmask;
	int DFmask;
	int **EToF;
	int **EToE;
	int *vmapM;
	int *vmapP;
	int dimmapB;
	int *mapB;
	int *vmapB;
	int mapI, mapO, vmapI, vmapO;
	double rk4a[5];
  	double rk4b[5];
  	double rk4c[5];
};

#endif //MESH1D_H_
