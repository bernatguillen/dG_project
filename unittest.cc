#include<stdlib.h>
#include<stdio.h>
#include<iostream>
#include "mesh1d.h"


/*Unit testing for the different functions developed for the dG method*/
int TestJacobiP(){
  const int N = 100;
  double x[N];
  for(int i = 0; i<N; ++i){
    x[i] = -1.0 + double(i)/(N-1)*2;
  }
  double **fx;
  fx = new double*[3];
  for(int n = 0; n<3; ++n){
    fx[n] = new double[N];
    JacobiP(n,1,1,x,fx[n],N);
  }
  double aux = 0;
  int val;
  
  for(int n1 = 0; n1<3; ++n1){
    for(int n2 = 0; n2<3; ++n2){
      val = (n1==n2);
      aux = 0;
      for(int i = 0; i<N-1; ++i){
	aux+=((1-x[i])*(1+x[i])*fx[n1][i]*fx[n2][i]+(1-x[i+1])*(1+x[i+1])*fx[n1][i+1]*fx[n2][i+1])*(x[i+1]-x[i])/2;
      }
      std::cout<<n1<<" times "<<n2<<" equals "<<aux<<std::endl;
      if(abs(aux-val)>1E-1) return 1;
    }
  }
  return 0;
}

int TestJacobiGQ(){
  const int N=4;
  double x[N+1];
  double w[N+1];
  double xref[] = {-0.830223896278566929872,-0.4688487934707142138038,0,0.4688487934707142138038,0.830223896278566929872};
  double wref[] = {0.086017682122807453242,0.336839460734335403901,0.487619047619047619048,0.336839460734335403901,0.086017682122807453242};

  JacobiGQ(N,1,1,x,w);
  for (int i = 0; i < N+ 1; ++i){
    std::cout<<"(xi,wi): " << x[i] << " " << w[i] <<std::endl;
    if(abs(xref[i]-x[i])>1e-8 || abs(wref[i]-w[i])>1e-8) return 1;
  }

  return 0;
}

int TestJacobiGL(){
  const int N = 4;
  double x[N+1];
  double xref[] = {-1,-0.6546536707079771437983,0,0.6546536707079771437983,1};
  JacobiGL(N,0,0,x);
  for (int i = 0; i < N+1; ++i){
    std::cout<<"xi :"<<x[i]<<std::endl;
    if(abs(xref[i]-x[i])>1e-8) return 1;
  }
  return 0;
}

int TestDmatrix1D(){
  const int N = 2;
  double x[N+1];
  double **V;
  double **D;
  double Dref[9]={-1.5,2,-0.5,-0.5,0,0.5,0.5,-2,1.5};
  V = new double*[N+1];
  D = new double*[N+1];
  for(int i = 0; i<N+1; ++i){
    V[i] = new double[N+1];
    D[i] = new double[N+1];
  }
  JacobiGL(N,0,0,x);
  Vandermonde1D(N,x,V,N+1);
  Dmatrix1D(N, x, V,D, N+1);
  for (int i = 0; i < N+1; ++i){
    for(int j = 0; j< N+1; ++j){
      std::cout << D[i][j] << "\t";
      if(abs(D[i][j]-Dref[3*i+j])>1e-8) return 1;
    }
    std::cout << std::endl;
  }
  for(int i = 0; i<N+1; ++i){
    delete [] V[i];
    delete [] D[i];
  }
  delete [] V;
  delete [] D;
  return 0;
}

int main(int argc, char **argv){
  if(!TestJacobiP()){
    std::cout<<"TestJacobiP passed!"<<std::endl;
  }else{
    std::cout<<"TestJacobiP failed!"<<std::endl;
  }
  if(!TestJacobiGQ()){
    std::cout<<"TestJacobiGQ passed!"<<std::endl;
  }else{
    std::cout<<"TestJacobiGQ failed!"<<std::endl;
  }
  if(!TestJacobiGL()){
    std::cout<<"TestJacobiGL passed!"<<std::endl;
  }else{
    std::cout<<"TestJacobiGL failed!"<<std::endl;
  }
  if(!TestDmatrix1D()){
    std::cout<<"TestDmatrix1D passed!"<<std::endl;
  }else{
    std::cout<<"TestDmatrix1D failed!"<<std::endl;
  }

  return(1);
}
