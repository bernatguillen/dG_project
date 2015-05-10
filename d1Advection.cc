#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <iomanip>

int main(int argc, char *argv[]) {
  /*Declaring variables m,nsteps,bd,solver,u,f*/
  if (argc != 3) {
    printf("USAGE: %s <meshsize> <T>\n", argv[0]); //p1 case, nsteps = 20
    exit(1);
  }
  const int m = atoi(argv[1]);
  const double T = atof(argv[2]);
  const int nsteps;

  double bd[2];
  bd[0] = 1.0;
  bd[1] = 3.0;

  JacobiIterator *solver;
  solver = new JacobiIterator(m);

  double u[m];
  double f[m];
  for(int i = 0; i<m; ++i){
    u[i] = float(i)/float(m-1);
  }
  Myrhs(u,f,m);
  /*Solving the equation*/
  for(int i = 0; i<m; ++i){
    u[i] = bd[0] + (bd[1]-bd[0])*float(i)/float(m-1);
  }
  solver->Solve(u,f,m,nsteps);
  std::string str1 = "elliptic";
  std::string str2 = argv[2];
  std::string str3 = ".dat";
  std::ofstream myfile((str1+str2+str3).c_str());
  if(myfile.is_open()){
    myfile << std::setprecision (15);
    for(int i = 0; i<m; ++i){
      myfile << u[i] << std::endl;
    }
  }
  delete solver;
  return 0;
}
