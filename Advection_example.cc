#include<stdlib.h>
#include<stdio.h>
#include<cmath>
#include<iostream>
#include<fstream>
#include "mesh1d.h"
#include "d1Advection.h"


int main(int argc, char **argv){
  int N = 8;
  double xmin = 0.0;
  double xmax = 2.0;
  const double pi = 3.141592653589793238463;
  int K = 10;
  int Nv;
  double Vx[5];
  int **EToV = new int*[K];
  for(int k =0 ; k<K; ++k){
    EToV[k] = new int[2];
  }
  MeshGen1D(xmin, xmax, K, &Nv, Vx, EToV);
  dG1D_Framework *mesh = new dG1D_Framework(K,Nv,N,Vx,EToV);
  
  double **u = new double*[mesh->Np];
  for(int i =0 ; i<mesh->Np; ++i){
    u[i] = new double[mesh->K];
    for(int j = 0; j< mesh->K; ++j){
      u[i][j] = sin(mesh->x[i][j]);
    }
  }
  double  FinalTime = 10;
  Advec1D(u,2*pi,FinalTime,mesh);
  std::ofstream file;
  file.open("Prova.dat");
  std::cout << "Writing to file..." <<std::endl;
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << u[i][j] << "\t";
    }
    file << std::endl;
  }
  std::cout << "File wrote!" << std::endl;

  for(int i = 0; i<mesh->Np; ++i){
    delete [] u[i];
  }
  delete [] u;
  for(int k = 0 ; k<K; ++k){
    delete [] EToV[k];
  }
  delete []EToV;
  return(1);
}
