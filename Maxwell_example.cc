#include<stdlib.h>
#include<stdio.h>
#include<cmath>
#include<iostream>
#include<fstream>
#include "mesh1d.h"
#include "d1Maxwell.h"


int main(int argc, char **argv){
  int N = 6;
  double xmin = -2.0;
  const double pi = 3.141592653589793238463;
  double xmax = 2.0;
  int K = 80;
  int Nv;
  double Vx[K+1];
  int **EToV = new int*[K];
  for(int k =0 ; k<K; ++k){
    EToV[k] = new int[2];
  }
  MeshGen1D(xmin, xmax, K, &Nv, Vx, EToV);
  dG1D_Framework *mesh = new dG1D_Framework(K,Nv,N,Vx,EToV);
  double **eps = new double*[mesh->Np];
  double **mu = new double*[mesh->Np];
  double **E = new double*[mesh->Np];
  double **H = new double*[mesh->Np];
  for(int i =0 ; i<mesh->Np; ++i){
    eps[i] = new double[mesh->K];
    mu[i] = new double[mesh->K];
    E[i] = new double[mesh->K];
    H[i] = new double[mesh->K];
    for(int j = 0; j< mesh->K; ++j){
      eps[i][j] = 1 + ((j+1) > mesh->K/2);
      mu[i][j] = 1;
      E[i][j] = sin(pi*mesh->x[i][j])*(mesh->x[i][j] < 0);
      H[i][j] = 0.0;
    }
  }
  std::ofstream file;
  file.open("x.dat");
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << mesh->x[i][j] << "\t";
    }
    file << std::endl;
  }
  file.close();  
  
  file.open("InitE.dat");
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << E[i][j] << "\t";
    }
    file << std::endl;
  }
  file.close();  
  file.open("InitH.dat");
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << H[i][j] << "\t";
    }
    file << std::endl;
  }
  file.close();  
  double  FinalTime = 10;
  Maxwell1D(E,H,eps,mu,FinalTime,mesh);

  file.open("EndE.dat");
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << E[i][j] << "\t";
    }
    file << std::endl;
  }
  file.close(); 
  file.open("EndH.dat");
  for(int i = 0; i<mesh->Np; ++i){
    for(int j = 0; j<mesh->K; ++j){
      file << H[i][j] << "\t";
    }
    file << std::endl;
  }
  file.close();  
 
  std::cout << "File wrote!" << std::endl;

  for(int i = 0; i<mesh->Np; ++i){
    delete [] E[i];
    delete [] H[i];
    delete [] eps[i];
    delete [] mu[i];
  }
  delete [] E;
  delete [] H;
  delete [] eps;
  delete [] mu;
  for(int k = 0 ; k<K; ++k){
    delete [] EToV[k];
  }
  delete []EToV;
  return(1);
}
