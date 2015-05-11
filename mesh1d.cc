#include<stdlib.h>
#include<cmath>
#include<gsl/gsl_linalg.h>
#include<gsl/gsl_eigen.h>
#include<gsl/gsl_sort.h>
#include<gsl/gsl_matrix.h>
#include<gsl/gsl_cblas.h>

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

int GradJacobiP(int N, double a, double b, double *x, double *fx, int npoints){
  //Derivative of JacobiP
  if(N==0){
    for(int i = 0; i<npoints; ++i){
      fx[i] = 0.0;
    }
  }else{
    JacobiP(N-1,a+1,b+1,x,fx,npoints);
    for(int i = 0; i<npoints; ++i){
      fx[i] *= sqrt(N*(N+a+b+1));
    }
  }
  return 1;
}


int JacobiGQ(int N, double a, double b, double *x, double *w){
  //Compute N'th order Gauss quadrature points, x, and weights, w, associated with JacobiP(a,b,N)
  if(N==0){
    x[0] = (a-b)/(a+b+2);
    w[0] = 2;
    return(1);
  }
  gsl_matrix *T = gsl_matrix_calloc(N+1,N+1);
  for(int i = 0; i<N+1; ++i){
    gsl_matrix_set(T,i,i,-(a*a-b*b)/((2*i+a+b)*(2*i+a+b+2)));
  }
  if(a+b <= 1e-15) gsl_matrix_set(T,0,0,0);
  for (int i = 1; i <N+1; ++i)
  {
    gsl_matrix_set(T,i,i-1,2.0/(2*i+a+b)*sqrt(i*(i+a+b)*(i+a)*(i+b)/((2*i+a+b-1)*(2*i+a+b+1))));
    gsl_matrix_set(T,i-1,i,gsl_matrix_get(T,i,i-1));
  }
  gsl_vector *eigval = gsl_vector_alloc(N+1);
  gsl_matrix *eigvec = gsl_matrix_alloc(N+1,N+1);
  gsl_eigen_symmv_workspace *wrk = gsl_eigen_symmv_alloc(N+1);
  gsl_eigen_symmv(T,eigval,eigvec,wrk);
  gsl_eigen_symmv_free(wrk);
  for (int i = 0; i < N+1; ++i){
    x[i] = gsl_vector_get(eigval,i);
    w[i] = pow(gsl_matrix_get(eigvec,0,i),2)*pow(2,a+b+1)/(a+b+1)*tgamma(a+1)*tgamma(b+1)/tgamma(a+b+1);
  }
  gsl_sort2(x,1,w,1,N+1);
  gsl_vector_free(eigval);
  gsl_matrix_free(eigvec);
  gsl_matrix_free(T);
  return(1);
}

int JacobiGL(int N, double a, double b, double *x){
  double w[N+1];
  if(N==1){
    x[0] = -1.0;
    x[1] = 1.0;
    return 1;
  }
  JacobiGQ(N-2,a+1,b+1,x,w);
  for (int i = N-1; i > 0; --i){
    x[i] = x[i-1];
  }
  x[0] = -1.0;
  x[N] = 1.0;
  return 1;
}

int InvMatrix(int N, double **M, double **M_inv){
  gsl_matrix *M_gsl = gsl_matrix_alloc(N,N);
  gsl_matrix *M_gsl_inv = gsl_matrix_alloc(N,N);
  gsl_permutation *perm = gsl_permutation_alloc(N);
  for(int i = 0; i<N; ++i){
    for(int j = 0; j<N; ++j){
      gsl_matrix_set(M_gsl, i,j, M[i][j]);
    }
  }
  int s;
  gsl_linalg_LU_decomp(M_gsl, perm, &s);
  gsl_linalg_LU_invert(M_gsl, perm, M_gsl_inv);
  for(int i =0; i<N; ++i){
    for(int j =0; j<N; ++j){
      M_inv[i][j] = gsl_matrix_get(M_gsl_inv, i,j);
    }
  }
  gsl_matrix_free(M_gsl_inv);
  gsl_matrix_free(M_gsl);
  gsl_permutation_free(perm);
  return 1;
}

int Vandermonde1D(int N, double *x, double **V, int npoints){
  double fx[N+1];
  for (int i = 0; i < npoints; ++i){
    for (int j = 0; j < N+1; ++j){
      V[i][j] = 0;
    }
  }
  for(int j = 0; j<N+1; ++j){
    JacobiP(j,0,0,x,fx,npoints);
    for(int i = 0; i<npoints; ++i){
      V[i][j] = fx[i];
    }
  }
  return 1;
}

int GradVandermonde1D(int N, double *x, double **V, int npoints){
  double fx[N+1];
  for (int i = 0; i < npoints; ++i){
    for (int j = 0; j < N+1; ++j){
      V[i][j] = 0;
    }
  }
  for(int j = 0; j<N+1; ++j){
    GradJacobiP(j,0,0,x,fx,npoints);
    for(int i = 0; i<npoints; ++i){
      V[i][j] = fx[i];
    }
  }
  return 1;  
}

int Dmatrix1D(int N, double *x, double **V, double **D, int npoints){
  double **Vr;
  Vr = new double*[N+1];
  for(int i = 0; i<N+1; ++i){
    Vr[i] = new double[N+1];
  }
  if(npoints!=N+1) return -1;
  GradVandermonde1D(N, x, Vr, npoints);
  gsl_matrix *V_gsl = gsl_matrix_alloc(N+1,N+1);
  gsl_matrix *V_inv = gsl_matrix_alloc(N+1,N+1);
  gsl_permutation *perm = gsl_permutation_alloc(N+1);
  for(int i = 0; i<N+1; ++i){
    for(int j = 0; j<N+1; ++j){
      gsl_matrix_set(V_gsl, i,j, V[i][j]);
    }
  }
  int s;
  gsl_linalg_LU_decomp(V_gsl, perm, &s);
  gsl_linalg_LU_invert(V_gsl, perm, V_inv);
  for(int i = 0; i<N+1; ++i){
    for(int j = 0; j<N+1; ++j){
      D[i][j] = 0;
      for(int k = 0; k<N+1; ++k){
        D[i][j] += Vr[i][k]*gsl_matrix_get(V_inv,k,j);
      }
    }
  }
  gsl_matrix_free(V_inv);
  gsl_matrix_free(V_gsl);
  gsl_permutation_free(perm);
  for(int i = 0; i<N+1; ++i){
    delete [] Vr[i];
  }
  delete [] Vr;
  return 1;
}

int MeshGen1D(double xmin, double xmax, int K, int *Nv, double *Vx, int **EToV){
  *Nv = K+1;
  for(int i = 0; i < *Nv;++i){
    Vx[i] = (xmax-xmin)*double(i)/K + xmin;
  }
  for(int k = 0; k < K; ++k){
    EToV[k][0] = k;
    EToV[k][1] = k+1;
  }
  return 1;
}

dG1D_Framework::dG1D_Framework(int K_in, int Nv_in, int N_in, double *Vx_in, int **EToV_in)
  : K(K_in),
  Nv(Nv_in),
  N(N_in){
    //Npoints per node
    Np = N+1;
    //Faces at point p, faces at node K
    Nfp = 1;
    Nfaces = 2;
    //Node tolerance
    NODETOL = 1e-10;
    Vx = new double[Nv];
    for(int i = 0; i<Nv; ++i) Vx[i] = Vx_in[i];
    EToV = new int* [K];
    for(int k = 0; k<K; ++k){
      EToV[k] = new int[2];
      EToV[k][0] = EToV_in[k][0];
      EToV[k][1] = EToV_in[k][1];
    }
    //r are Gauss Lobatto quadrature points
    r = new double[Np];
    JacobiGL(N,0,0,r);
    //Vandermonde and inverse of Vandermonde Matrix
    V = new double*[Np];
    Vinv = new double*[Np];
    for(int i = 0; i<Np; ++i){
      V[i]=new double[Np];
      Vinv[i]=new double[Np];
    }
    Vandermonde1D(N, r, V, Np);
    InvMatrix(N,V,Vinv);
    //Dr derivative matrix
    Dr = new double*[Np];
    for(int i = 0; i<Np; ++i){
      Dr[i] = new double[Np];
    }
    Dmatrix1D(N, r, V, Dr, Np);
    //LIFT Matrix (Surface (boundary) terms) V*Vt*Emat
    //Emat = zeros(Np,Nfaces*Nfp) , [1 0 0 0 ... 0],[0 0 ... 0 1]-> In this case it is not necessary to calculate
    LIFT = new double*[Np];
    for(int i = 0;i <Np; ++i){
      LIFT[i] = new double[2];
      LIFT[i][0] = 0;
      LIFT[i][1] = 0;
    }
    for(int i = 0; i<Np; ++i){
      for(int j = 0; j<Np; ++j){
        LIFT[i][0] += V[0][j]*V[i][j];
        LIFT[i][1] += V[Np-1][j]*V[i][j];
      }
    }
    //create x, matrix of physical points (affine traslation of r for each node), each column is a node
    x = new double*[Np];
    for(int i = 0; i<Np;++i){
      x[i] = new double[K];
    }
    for(int k = 0; k<K; ++k){
      for(int i = 0; i<Np; ++i){
        x[i][k] = Vx[EToV[k][0]] + 0.5*(1+r[i])*(Vx[EToV[k][1]]-Vx[EToV[k][0]]);
      }
    }
    //create geometric factors rx and J
    rx = new double*[Np];
    J = new double*[Np];
    for(int i = 0; i<Np;++i){
      rx[i] = new double[K];
      J[i] = new double[K];
      for(int k = 0; k<K; ++k){
        J[i][k] = 0;
        for(int j = 0; j<Np; ++j){
          J[i][k]+=Dr[i][j]*x[j][k];
        }
        rx[i][k] = 1/J[i][k];
      }
    }
    //create mask and masked vectors/matrices
    int res1 = 0; //(finding how many times r = -1)
    int res2 = 0; //(finding how many times r = 1)
    for(int i = 0; i<Np; ++i){
      res1+= (abs(r[i]+1)<NODETOL);
      res2+= (abs(r[i]-1)<NODETOL);
    }
    Fmask = new int[res1+res2];
    DFmask = res1+res2;
    int i1 = 0;
    int i2 = res1;
    for(int i = 0; i<Np; ++i){
      if(abs(r[i]+1)<NODETOL){
        Fmask[i1] = i;
        ++i1;
      } 
      if(abs(r[i]-1)<NODETOL){
        Fmask[i2] = i;
        ++i2;
      }
    }
    Fx = new double*[DFmask];
    Fscale = new double*[DFmask];
    for(int d = 0; d<DFmask; ++d){
      Fx[d] = new double[K];
      Fscale[d] = new double[K];
      for(int k = 0; k<K; ++k){
        Fx[d][k] = x[Fmask[d]][k];
        Fscale[d][k] = rx[Fmask[d]][k];
      }
    }
    //Surface normals
    nx = new double*[Nfp*Nfaces];
    for(int n=0; n<Nfp*Nfaces; ++n){
      nx[n] = new double[K];
      for(int k = 0; k<K; ++k){
        nx[0][k] = -1.0;
        nx[1][k] = 1.0;
      }
    }
    //Connectivity arrays
    int **SpFToV = new int*[2*K*Nfaces];
    for(int i = 0; i<2*K*Nfaces; ++i){
      SpFToV[i] = new int[2];
    }
    int sk = 0;
    for(int k = 0; k<K; ++k){
      for(int fac = 0; fac<Nfaces; ++fac){
        SpFToV[sk][0] = sk;
        SpFToV[sk][1] = EToV[k][fac];
        ++sk; 
      }
    }
    int **SpFToF = new int *[2*K*Nfaces];
    for(int i = 0; i<2*Nfaces*K; ++i){
      SpFToF[i] = new int[2];
    }
    sk = 0;
    for(int skv = 0; skv<K*Nfaces; ++skv){
      for(int skvaux = 0; skvaux<K*Nfaces; ++skvaux){
        if(skvaux != skv && SpFToV[skv][1] == SpFToV[skvaux][1]){
          SpFToF[sk][0] = skv;
          SpFToF[sk][1] = skvaux;
          ++sk;
        }
      }
    }
    for(int i =0 ; i<2*K*Nfaces; ++i){
      delete [] SpFToV[i];
    }
    delete [] SpFToV;
    
  }



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
