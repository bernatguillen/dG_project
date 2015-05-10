#ifndef MESH1D_H_
#define MESH1D_H_

#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>

int JacobiP(int N, double a, double b, double *x,double *fx, int npoints);
int JacobiGQ(int N, double a, double b, double *x, double *w);
int JacobiGL(int N, double a, double b, double *x, double *w);

class JacobiIterator{
 public:
  JacobiIterator(int m);
  ~JacobiIterator();
  int Step(double *x, double *b, int m);
  int Solve(double *x0, double *b, int nsteps, int stride);
 private:
  const int m_;
  double h_;
  double *xaux_;
};

#endif //MESH1D_H_
