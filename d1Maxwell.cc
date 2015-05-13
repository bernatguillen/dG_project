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
#include "d1Maxwell.h"

int printmatrix(double **u, int n, int m){
  for(int i = 0; i<n; ++i){
    for(int j = 0; j<m; ++j){
      printf("%f\t",u[i][j]);
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  return 1;
}
int printvector(double *u, int n){
  for(int i = 0; i<n; ++i){
      printf("%f\t",u[i]);
  }
  std::cout<<std::endl<<std::endl;
  return 1;
}
int printvector2(int *u, int n){
  for(int i = 0; i<n; ++i){
      printf("%d\t",u[i]);
  }
  std::cout<<std::endl<<std::endl;
  return 1;
}
int printmatrixtofile(double **u, int n, int m, const char *filename){
  std::ofstream file;
  file.open(filename);
  for(int i = 0; i<n; ++i){
    for(int j = 0; j<m; ++j){
      file << u[i][j] << "\t";
    }
    file<<std::endl;
  }
  file.close();
  return 1;
}


int printmatrix2(int **u, int n, int m){
  for(int i = 0; i<n; ++i){
    for(int j = 0; j<m; ++j){
      printf("%d\t",u[i][j]);
    }
    std::cout<<std::endl;
  }
  std::cout<<std::endl;
  return 1;
}

int MaxwellRHS1D(double **E, double **H, double **eps, double **mu, dG1D_Framework *mesh, double **rhsE, double **rhsH){
  double **dE = new double*[mesh->Nfp*mesh->Nfaces];
  double **dH = new double*[mesh->Nfp*mesh->Nfaces];
  double **Zimpm = new double*[mesh->Nfp*mesh->Nfaces];
  double **Zimpp = new double*[mesh->Nfp*mesh->Nfaces];
  double **Yimpm = new double*[mesh->Nfp*mesh->Nfaces];
  double **Yimpp = new double*[mesh->Nfp*mesh->Nfaces];
  double **fluxE = new double*[mesh->Nfp*mesh->Nfaces];
  double **fluxH = new double*[mesh->Nfp*mesh->Nfaces];

  //1/Z = Y
for(int i = 0; i<mesh->Nfp*mesh->Nfaces; ++i){
    dE[i] = new double[mesh->K];
    dH[i] = new double[mesh->K];
    fluxE[i] = new double[mesh->K];
    fluxH[i] = new double[mesh->K];
    Zimpm[i] = new double[mesh->K];
    Zimpp[i] = new double[mesh->K];
    Yimpm[i] = new double[mesh->K];
    Yimpp[i] = new double[mesh->K];
  }

  int iM, jM, iP, jP,i,j;
  for(int m = 0; m<mesh->Nfp*mesh->Nfaces*mesh->K; ++m){
    iM = mesh->vmapM[m]%(mesh->Np);
    jM = mesh->vmapM[m]/mesh->Np;
    iP = mesh->vmapP[m]%mesh->Np;
    jP = mesh->vmapP[m]/mesh->Np;
    i = m%(mesh->Nfp*mesh->Nfaces);
    j = m/(mesh->Nfp*mesh->Nfaces);
//    printf("m=%d,i=%d,j=%d,iM=%d,jM=%d,iP=%d,jP=%d\n",m,i,j,iM,jM,iP,jP);
    dE[i][j] = (E[iM][jM]-E[iP][jP]);
    dH[i][j] = (H[iM][jM]-H[iP][jP]);
    Zimpm[i][j] = sqrt(mu[iM][jM]/eps[iM][jM]);
    Zimpp[i][j] = sqrt(mu[iP][jP]/eps[iP][jP]);
    Yimpm[i][j] = 1.0/Zimpm[i][j];
    Yimpp[i][j] = 1.0/Zimpp[i][j];
  }

  double Ebc[mesh->dimmapB];
  double Hbc[mesh->dimmapB];
  int iB,jB;
  for(int m = 0; m<mesh->dimmapB;++m){
    iB = mesh->vmapB[m]%(mesh->Np);
    jB = mesh->vmapB[m]/mesh->Np;
    i = mesh->mapB[m]%(mesh->Nfp*mesh->Nfaces);
    j = mesh->mapB[m]/(mesh->Nfp*mesh->Nfaces);
    Ebc[m] = -E[iB][jB];
    Hbc[m] = H[iB][jB];
    dE[i][j] = E[iB][jB] -Ebc[m];
    dH[i][j] = H[iB][jB] - Hbc[m];
  }
  for(int i = 0; i<mesh->Nfp*mesh->Nfaces; ++i){
    for(int j = 0; j<mesh->K;++j){
      fluxE[i][j] = 1.0/(Zimpm[i][j]+Zimpp[i][j])*(mesh->nx[i][j]*Zimpp[i][j]*dH[i][j]-dE[i][j]);
      fluxH[i][j] = 1.0/(Yimpm[i][j]+Yimpp[i][j])*(mesh->nx[i][j]*Yimpp[i][j]*dE[i][j]-dH[i][j]);
    }
  }
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      rhsE[i][j] = 0;
      rhsH[i][j] = 0;
      for(int l = 0; l<mesh->Nfp*mesh->Nfaces; ++l){
        rhsE[i][j] += mesh->LIFT[i][l]*mesh->Fscale[l][j]*fluxE[l][j];
        rhsH[i][j] += mesh->LIFT[i][l]*mesh->Fscale[l][j]*fluxH[l][j];
      }
      for(int l = 0; l<mesh->Np; ++l){
        rhsE[i][j] -= mesh->Dr[i][l]*H[l][j]*mesh->rx[i][j];
        rhsH[i][j] -= mesh->Dr[i][l]*E[l][j]*mesh->rx[i][j];
      }
      rhsE[i][j] /= eps[i][j];
      rhsH[i][j] /= mu[i][j];
    }
  }
  for(int i = 0; i<mesh->Nfp*mesh->Nfaces; ++i){
    delete [] dH[i];
    delete [] dE[i];
    delete [] Zimpp[i];
    delete [] Zimpm[i];
    delete [] Yimpm[i];
    delete [] Yimpp[i];
    delete [] fluxH[i];
    delete [] fluxE[i];
  }
  delete [] dH;
  delete [] dE;
  delete [] Zimpm;
  delete [] Zimpp;
  delete [] Yimpm;
  delete [] Yimpp;
  delete [] fluxE;
  delete [] fluxH;
  return 1;
}

int Maxwell1D(double **E, double **H, double **eps, double **mu, double FinalTime, dG1D_Framework *mesh){
  double **rhsE = new double*[mesh->Np];
  double **resE = new double*[mesh->Np];
  double **rhsH = new double*[mesh->Np];
  double **resH = new double*[mesh->Np];
  for(int i = 0; i<mesh->Np; ++i){ 
    rhsE[i] = new double[mesh->K];
    resE[i] = new double[mesh->K];
    rhsH[i] = new double[mesh->K];
    resH[i] = new double[mesh->K];
    for(int j = 0; j<mesh->K; ++j){
      resE[i][j] = 0.0;
      resH[i][j] = 0.0;
    }
  }

  double xmin = std::numeric_limits<double>::max();
  for(int i = 1; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      if(fabs(mesh->x[i][j]-mesh->x[i-1][j])<xmin) xmin = fabs(mesh->x[i][j]-mesh->x[i-1][j]);
    }
  }
  double CFL = 0.75;
  double dt = CFL*xmin;
  int Nsteps = int(FinalTime/dt) + 2;
  dt = FinalTime/(Nsteps-1);
  for(int tstep = 0; tstep<Nsteps; ++tstep){
    for(int INTRK = 0; INTRK<5; INTRK++){
      MaxwellRHS1D(E,H,eps,mu,mesh,rhsE,rhsH);
      for(int i = 0; i<mesh->Np; ++i){
        for(int j = 0; j<mesh->K; ++j){
          resE[i][j] = mesh->rk4a[INTRK]*resE[i][j] + dt*rhsE[i][j];
          resH[i][j] = mesh->rk4a[INTRK]*resH[i][j] + dt*rhsH[i][j];
          E[i][j] += mesh->rk4b[INTRK]*resE[i][j];
          H[i][j] += mesh->rk4b[INTRK]*resH[i][j];
        }
      }
    }
  }
  for(int i = 0; i<mesh->Np; ++i){ 
    delete [] rhsE[i];
    delete [] resE[i];
    delete [] rhsH[i];
    delete [] resH[i];
  }  
  delete [] rhsE;
  delete [] resE;
  delete [] rhsH;
  delete [] resH;
  return 1;
}