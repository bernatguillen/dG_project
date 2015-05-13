#ifndef D1ADV_H_
#define D1ADV_H_

#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include "mesh1d.h"

int MaxwellRHS1D(double **E, double **H, double **eps, double **mu, dG1D_Framework *mesh, double **rhsE, double **rhsH);
int Maxwell1D(double **E, double **H, double **eps, double **mu, double FinalTime, dG1D_Framework *mesh);

#endif //D1ADV_H_
