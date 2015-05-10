#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg>
#include "mesh1d.h"

int JacobiP(int N, double a, double b, double *x,double *fx, int npoints){
  const double Gamma0 = pow(2,a+b+1)/(a+b+1)*tgamma(a+1)*tgamma(b+1)/tgamma(a+b+1);
  double **PL;
  PL = new double*[N+1];
  for(int n = 0; n<N+1; ++n){
    PL[n] = new double[npoints];
  }
  //Initial values N=0 , 1
  for(int i = 0; i<npoints; ++i){
    PL[0][i] = 1.0/sqrt(Gamma0);
  }
  if(N==0){
    for(int i = 0; i<npoints; ++i){
      fx[i] = PL[0][i];
    }
    return(1);
  }  
  double Gamma1 = Gamma0*(a+1)*(b+1)/(a+b+3);
  for(int i = 0; i<npoints; ++i){
    PL[1][i] = ((a+b+2)*x[i]/2+(a-b)/2)/sqrt(Gamma1);
  }
  if(N==1){
    for(int i = 0; i<npoints; ++i){
      fx[i] = PL[1][i];
    }
    return(1);
  }
  double aold = 2/(2+a+b)*sqrt((a+1)*(b+1)/(a+b+3));
  double anew,h1,bnew;
  for(int n = 0; n<N-1; ++n){
    h1 = 2*(n+1)+a+b;
    anew = 2/(h1+2)*sqrt((n+2)*(n+2+a+b)*(n+2+a)*(n+2+b)/(h1+1)/(h1+3));
    bnew = - (a*a - b*b)/h1/(h1+2);
    for(int i = 0; i<npoints; ++i){
      PL[n+2][i] = 1.0/anew*(-aold*PL[n][i]+(x[i]-bnew)*PL[n+1][i]);
    }
  }
  for(int i = 0; i<npoints; ++i){
    fx[i] = PL[N][i];
  }
  for(int n = 0; n<N+1; ++n){
    delete [] PL[n];
  }
  delete [] PL;
  return(1);
}
int JacobiGQ(int N, double a, double b, double *x, double *w){
  //Compute N'th order Gauss quadrature points, x, and weights, w, associated with JacobiP(a,b,N)
  if(N==0){
    x[0] = (a-b)/(a+b+2);
    w[0] = 2;
    return(1);
  }

}
int JacobiGL(int N, double a, double b, double *x, double *w);

JacobiIterator::JacobiIterator(int m)
  : m_(m){
  xaux_ = new double[m_];
  h_ = 1.0/float(m_-1);
}

JacobiIterator::~JacobiIterator() {
  delete [] xaux_;
}
int JacobiIterator::Step(double *x,double *b, int m){
  double h = 1.0/double(m-1);
  for(int i = 1; i<m-1;++i){
    xaux_[i] = 0.5*(x[i-1]+x[i+1]-h*h*b[i]);
  }
  for(int i = 1; i<m-1; ++i){
    x[i] = xaux_[i];
  }
  return 1;
}
int JacobiIterator::Solve(double *x0, double *b,int m, int nsteps){
  for(int k = 0; k<nsteps; ++k){
    Step(x0,b,m);
  }
  return 1;
}
