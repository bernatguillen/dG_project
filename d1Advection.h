#ifndef D1ADV_H_
#define D1ADV_H_

#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include "mesh1d.h"

int AdvecRHS1D(double **u, double time, double a, dG1D_Framework *mesh, double **rhsu);
int Advec1D(double **u, double a, double FinalTime, dG1D_Framework *mesh);

#endif //D1ADV_H_
