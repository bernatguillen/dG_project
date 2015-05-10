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

int main(int argc, char **argv){
  if(!TestJacobiP()){
    std::cout<<"TestJacobiP passed!"<<std::endl;
  }else{
    std::cout<<"TestJacobiP failed!"<<std::endl;
  }
  return(1);
}
