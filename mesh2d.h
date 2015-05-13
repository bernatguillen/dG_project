#ifndef MESH2D_H_
#define MESH2D_H_

#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include "mesh1d.h"

int Simplex2DP(int i, int j, double **a, double **fx, int Npoints);


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
