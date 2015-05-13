#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_cblas.h>
#include<limits>
#include<iostream>
#include<fstream>

#include "mesh1d.h"
#include "d1discont.h"

int AdvecRHS1D(double **u, double t, double a, dG1D_Framework *mesh, double **rhsu){
  double alpha = 0.5;
  double **du = new double*[mesh->Nfp*mesh->Nfaces];
  for(int i = 0; i<mesh->Nfp*mesh->Nfaces; ++i){
    du[i] = new double[mesh->K];
  }
  int iM, jM, iP, jP,i,j;
  for(int m = 0; m<mesh->Nfp*mesh->Nfaces*mesh->K; ++m){
    iM = mesh->vmapM[m]%(mesh->Np);
    jM = mesh->vmapM[m]/mesh->Np;
    iP = mesh->vmapP[m]%mesh->Np;
    jP = mesh->vmapP[m]/mesh->Np;
    i = m%(mesh->Nfp*mesh->Nfaces);
    j = m/(mesh->Nfp*mesh->Nfaces);
 //   printf("m=%d,i=%d,j=%d,iM=%d,jM=%d,iP=%d,jP=%d\n",m,i,j,iM,jM,iP,jP);
    du[i][j] = (u[iM][jM]-u[iP][jP])*(a*mesh->nx[i][j]-(1-alpha)*fabs(a*mesh->nx[i][j]))/2.0;
  }
  i = mesh->vmapO%(mesh->Np);
  j = mesh->vmapO/(mesh->Np);
  double uin = u[i][j];
  i = mesh->mapI%(mesh->Nfp*mesh->Nfaces);
  j = mesh->mapI/(mesh->Nfp*mesh->Nfaces);
  iM = mesh->vmapI%mesh->Np;
  jM = mesh->vmapI/mesh->Np;
// printf("%d,%d,%d\n",i,j,mesh->mapI);
  du[i][j] = (u[iM][jM]-uin)*(a*mesh->nx[i][j]-(1-alpha)*fabs(a*mesh->nx[i][j]))/2.0;
  i = mesh->mapO%(mesh->Nfp*mesh->Nfaces);
  j = mesh->mapO/(mesh->Nfp*mesh->Nfaces);
//  printf("%d,%d,%d\n",i,j,mesh->mapO);
  du[i][j] = 0.0;
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      rhsu[i][j] = 0;
      for(int l = 0; l<mesh->Nfp*mesh->Nfaces; ++l){
        rhsu[i][j] += mesh->LIFT[i][l]*mesh->Fscale[l][j]*du[l][j];
      }
      for(int l = 0; l<mesh->Np; ++l){
        rhsu[i][j] -= a*mesh->Dr[i][l]*u[l][j]*mesh->rx[i][j];
      }
    }
  }
  for(int i = 0; i<mesh->Nfp*mesh->Nfaces; ++i){
    delete [] du[i];
  }
  delete [] du;
  return 1;
}

int Advec1D(double **u, double a, double FinalTime, dG1D_Framework *mesh){
  double **rhsu = new double*[mesh->Np];
  double **resu = new double*[mesh->Np];
  for(int i = 0; i<mesh->Np; ++i){ 
    rhsu[i] = new double[mesh->K];
    resu[i] = new double[mesh->K];
    for(int j = 0; j<mesh->K; ++j){
      resu[i][j] = 0.0;
    }
  }

  double t = 0;
  double timelocal;
  double xmin = std::numeric_limits<double>::max();
  for(int i = 1; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      if(fabs(mesh->x[i][j]-mesh->x[i-1][j])<xmin) xmin = fabs(mesh->x[i][j]-mesh->x[i-1][j]);
    }
  }

  double CFL = 0.75;
  double dt = CFL/(2*a)*xmin;
  
  int Nsteps = int(FinalTime/dt) + 2;
  dt = FinalTime/(Nsteps-1);
  for(int tstep = 0; tstep<Nsteps; ++tstep){
    t = tstep*dt;
    for(int INTRK = 0; INTRK<5; INTRK++){
      timelocal = t + mesh->rk4c[INTRK]*dt;
      AdvecRHS1D(u,timelocal,a,mesh,rhsu);
      for(int i = 0; i<mesh->Np; ++i){
        for(int j = 0; j<mesh->K; ++j){
          resu[i][j] = mesh->rk4a[INTRK]*resu[i][j] + dt*rhsu[i][j];
          u[i][j] += mesh->rk4b[INTRK]*resu[i][j];
        }
      }
    }
  }
  for(int i = 0; i < mesh->Np; ++i){
    delete [] rhsu[i];
    delete [] resu[i];
  }
  delete [] rhsu;
  delete [] resu;
  return 1;
}