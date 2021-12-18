/***************************************************************************                               
                    Molecular Dynamics Potentials (MDP)
                           CESMIX-MIT Project  
 
 Contributing authors: Ngoc-Cuong Nguyen (cuongng@mit.edu, exapde@gmail.com)
 ***************************************************************************/

// clang++ -std=c++11 -Wall -Wextra -pedantic -c -fPIC cpuSnap.cpp -o cpusnap.o
// clang++ -shared cpuSnap.o -o cpuSnap.dylib

// snaplib = Libdl.dlopen("cpusnap.dylib")
// idxbcount = Libdl.dlsym(snaplib, :idxbcount)
// count = ccall(idxbcount, Cint, (Cint,), 8)
// idxucount = Libdl.dlsym(snaplib, :idxucount)
// count = ccall(idxucount, Cint, (Cint,), 8)
// initsna = Libdl.dlsym(snaplib, :cpuInitSna)
// Libdl.dlclose(snaplib)

#ifndef CPUSNAP
#define CPUSNAP

#include "cpuSnap.h"

#include <stdio.h>
#include <math.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))
#define MAX(a,b) ((a) > (b) ? (a) : (b))

int idxbcount(int twojmax)
{
    int count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
                if (j >= j1) count++;

    return count;
};

int idxcgcount(int twojmax)
{
    int count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
                for (int m1 = 0; m1 <= j1; m1++)
                    for (int m2 = 0; m2 <= j2; m2++)
                        count++;
            }
    return count;
};

int idxzcount(int twojmax)
{
    int count = 0;
    for(int j1 = 0; j1 <= twojmax; j1++)
        for(int j2 = 0; j2 <= j1; j2++)
            for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
                for (int mb = 0; 2*mb <= j; mb++)
                    for (int ma = 0; ma <= j; ma++)
                        count++;

    return count;
};

int idxucount(int twojmax)
{
    int count = 0;
    for(int j = 0; j <= twojmax; j++) 
        for(int mb = 0; mb <= j; mb++)
            for(int ma = 0; ma <= j; ma++)
                count++;
  
    return count;
};

void cpuBuildIndexList(int *idx_max, int *idxz, int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, 
        int *idxcg_block, int twojmax)
{
  // index list for cglist

  int jdim = twojmax + 1;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxcg_block[j + j2*jdim + j1*jdim*jdim] = idxcg_count;  
        for (int m1 = 0; m1 <= j1; m1++)
          for (int m2 = 0; m2 <= j2; m2++)
            idxcg_count++;
      }
  idx_max[0] = idxcg_count;
          
  int idxu_count = 0;

  for(int j = 0; j <= twojmax; j++) {
    idxu_block[j] = idxu_count;
    for(int mb = 0; mb <= j; mb++)
      for(int ma = 0; ma <= j; ma++)
        idxu_count++;
  }
  //idxu_max = idxu_count;
  idx_max[1] = idxu_count;
  
  // index list for beta and B

  int idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) idxb_count++;

  int idxb_max = idxb_count;
  idx_max[2] = idxb_max;

  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        if (j >= j1) {
          idxb[idxb_count*3 + 0] = j1;
          idxb[idxb_count*3 + 1] = j2;
          idxb[idxb_count*3 + 2] = j;  
          idxb_count++;
        }
  
  idxb_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        if (j >= j1) {
          idxb_block[j + j2*jdim + j1*jdim*jdim] = idxb_count;    
          idxb_count++;
        }
      }

  // index list for zlist

  int idxz_count = 0;

  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2)
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++)
            idxz_count++;

  int idxz_max = idxz_count;
  idx_max[3] = idxz_max;

  idxz_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        idxz_block[j + j2*jdim + j1*jdim*jdim] = idxz_count;    

        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {

            idxz[idxz_count*10 + 0] = j1;
            idxz[idxz_count*10 + 1] = j2;
            idxz[idxz_count*10 + 2] = j;
            idxz[idxz_count*10 + 3] = MAX(0, (2 * ma - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 4] = (2 * ma - j - (2 * idxz[idxz_count*10 + 3] - j1) + j2) / 2;
            idxz[idxz_count*10 + 5] = MIN(j1, (2 * ma - j + j2 + j1) / 2) - idxz[idxz_count*10 + 3] + 1;
            idxz[idxz_count*10 + 6] = MAX(0, (2 * mb - j - j2 + j1) / 2);
            idxz[idxz_count*10 + 7] = (2 * mb - j - (2 * idxz[idxz_count*10 + 6] - j1) + j2) / 2;
            idxz[idxz_count*10 + 8] = MIN(j1, (2 * mb - j + j2 + j1) / 2) - idxz[idxz_count*10 + 6] + 1;
            // apply to z(j1,j2,j,ma,mb) to unique element of y(j)

            const int jju = idxu_block[j] + (j+1)*mb + ma;
            idxz[idxz_count*10 + 9] = jju;
              
            idxz_count++;
          }
      }
};

void cpuInitRootpqArray(double *rootpqarray, int twojmax)
{
  int jdim = twojmax + 1;
  for (int p = 1; p <= twojmax; p++)
    for (int q = 1; q <= twojmax; q++)
      rootpqarray[p*jdim + q] = sqrt(static_cast<double>(p)/q);
};

double cpuDeltacg(double *factorial, int j1, int j2, int j)
{
  double sfaccg = factorial[(j1 + j2 + j) / 2 + 1];
  return sqrt(factorial[(j1 + j2 - j) / 2] *
              factorial[(j1 - j2 + j) / 2] *
              factorial[(-j1 + j2 + j) / 2] / sfaccg);
};

void cpuInitClebschGordan(double *cglist, double *factorial, int twojmax)
{
  double sum,dcg,sfaccg;
  int m, aa2, bb2, cc2;
  int ifac;

  int idxcg_count = 0;
  for(int j1 = 0; j1 <= twojmax; j1++)
    for(int j2 = 0; j2 <= j1; j2++)
      for(int j = j1 - j2; j <= MIN(twojmax, j1 + j2); j += 2) {
        for (int m1 = 0; m1 <= j1; m1++) {
          aa2 = 2 * m1 - j1;

          for (int m2 = 0; m2 <= j2; m2++) {

            bb2 = 2 * m2 - j2;
            m = (aa2 + bb2 + j) / 2;

            if(m < 0 || m > j) {
              cglist[idxcg_count] = 0.0;
              idxcg_count++;
              continue;
            }

            sum = 0.0;

            for (int z = MAX(0, MAX(-(j - j2 + aa2)
                                    / 2, -(j - j1 - bb2) / 2));
                 z <= MIN((j1 + j2 - j) / 2,
                          MIN((j1 - aa2) / 2, (j2 + bb2) / 2));
                 z++) {
              ifac = z % 2 ? -1 : 1;
              sum += ifac /
                (factorial[z] *
                 factorial[(j1 + j2 - j) / 2 - z] *
                 factorial[(j1 - aa2) / 2 - z] *
                 factorial[(j2 + bb2) / 2 - z] *
                 factorial[(j - j2 + aa2) / 2 + z] *
                 factorial[(j - j1 - bb2) / 2 + z]);
            }

            cc2 = 2 * m - j;
            dcg = cpuDeltacg(factorial, j1, j2, j);
            sfaccg = sqrt(factorial[(j1 + aa2) / 2] *
                          factorial[(j1 - aa2) / 2] *
                          factorial[(j2 + bb2) / 2] *
                          factorial[(j2 - bb2) / 2] *
                          factorial[(j  + cc2) / 2] *
                          factorial[(j  - cc2) / 2] *
                          (j + 1));

            cglist[idxcg_count] = sum * dcg * sfaccg;
            idxcg_count++;
          }
        }
      }
}

void cpuInitSna(double *rootpqarray, double *cglist, double *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax)
{
    cpuBuildIndexList(idx_max, idxz, idxz_block, idxb, 
            idxb_block, idxu_block, idxcg_block, twojmax);
    
    cpuInitRootpqArray(rootpqarray, twojmax);
    cpuInitClebschGordan(cglist, factorial, twojmax);        
}

void cpuNeighPairList(int *pairnum, int *pairlist, double *x, double *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int *alist, int inum,  int dim, int ntypes)
{    
    int sumcount = 0;
    for (int ii=0; ii<inum; ii++) {
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        //int m = neighnum[i];     // number of neighbors around i      
        int n1 = neighnum[i];    
        int m = neighnum[i+1] - n1;
        int count = 0;              
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i {
            //int j = neighlist[l + jnum*i];         
            int j = neighlist[n1 + l];         
            int jtype  = atomtype[alist[j]];        
            // distance between atom i and atom j                                    
            double xij0 = x[j*dim] - x[i*dim];  // xj - xi
            double xij1 = x[j*dim+1] - x[i*dim+1]; // xj - xi               
            double xij2 = x[j*dim+2] - x[i*dim+2]; // xj - xi               
            double dij = xij0*xij0 + xij1*xij1 + xij2*xij2;                        
            if (dij < rcutsq[jtype + itype*(ntypes+1)] && dij>1e-20) {
                pairlist[count + sumcount] = j;  // atom j     
                count += 1;
            }
        }        
        pairnum[ii] = count;       
        sumcount += count; 
    }    
};


void cpuNeighPairs(double *xij, double *x, int *aii, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *ilist, int *atomtype, int *alist, int inum, int dim)
{        
    for (int ii=0; ii<inum; ii++) {  // for each atom i in the simulation box     
        int i = ilist[ii];       // atom i
        int itype = atomtype[i];
        //int m = pairnum[ii];        // number of neighbors around i             
        int start = pairnumsum[ii];   
        int m = pairnumsum[ii+1] - start;
        for (int l=0; l<m ; l++) {   // loop over each atom around atom i
            int j = pairlist[l + start];  // atom j              
            int k = start + l;                                     
            aii[k]       = ii;
            ai[k]        = i;
            aj[k]        = alist[j];          
            ti[k]        = itype;       
            tj[k]        = atomtype[alist[j]];        
            for (int d=0; d<dim; d++) 
                xij[k*dim+d]   = x[j*dim+d] -  x[i*dim+d];  // xj - xi            
        }
    }    
};

void cpuSnapComputeUi(double *Utotr, double *Utoti, double *rootpqarray, double *rij, double *wjelem, double *radelem, 
        double rmin0, double rfac0, double rcutfac, int *map, int *aii, int *ti, int *tj, 
        int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag)                
{    
  for(int ij=0; ij<ijnum; ij++) {              
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];    
    double r = sqrt(x * x + y * y + z * z);

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    //double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    //double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
            
    double sfac = 0.0;
    if (switchflag == 0) 
        sfac = 1.0;    
    else if (switchflag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
        }
        else if(r > rcutij) {
            sfac = 1.0;
        }
        else {
            double rcutfac = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
        }
    } 
    sfac *= wjelem[tj[ij]];
    
    double a_r, a_i, b_r, b_i;
    double rootpq;

    //double r0inv;    
    rcutij = 1.0 / sqrt(r * r + z0 * z0);
    a_r = rcutij * z0;
    a_i = -rcutij * z;
    b_r = rcutij * y;
    b_i = -rcutij * x;

    // 2Jmax = 10
    double Pr[11], Pi[11], Qr[9], Qi[9];
    Pr[0] = 1.0;
    Pi[0] = 0.0;    
    
    int jdim = twojmax + 1;
    int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
    int i = aii[ij] + njelem;                        
    Utotr[i] += sfac; // atomic add   
    
    //printf("%i %i %i\n", ij, i, twojmax);  

    int mb = 0;    
    for (int j = 1; j <= twojmax; j++) {        
        // fill in left side of matrix layer from previous layer
        int ma = 0;
        // x y z z0 
        // double p_r, p_i; // -> x, y
        // double u_r = Pr[ma]; // -> z
        // double u_i = Pi[ma]; // -> z0
        z = Pr[ma];
        z0 = Pi[ma];
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        Pr[ma] = rootpq * (a_r * z + a_i * z0);
        Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
        for (ma = 1; ma < j; ma++) {
            x = Pr[ma];
            y = Pi[ma];
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
            Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
            z = x;
            z0 = y;
        }
        ma = j;
        rcutij = rootpqarray[ma*jdim + (j - mb)];
        Pr[ma] = -rcutij * (b_r * z + b_i * z0);
        Pi[ma] = -rcutij * (b_r * z0 - b_i * z);                        
                                
        if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
            int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
            for (int ma = 0; ma <= j; ma++) {     
                if (mapar == 1) {                    
                    Qr[j-ma] = Pr[ma];
                    Qi[j-ma] = -Pi[ma];
                } else {
                    Qr[j-ma] = -Pr[ma];
                    Qi[j-ma] =  Pi[ma];
                }
                mapar = -mapar;
            }                                                
        }
        
        int k =  1 + (j+1)*mb;
        for (int ma = 2; ma <= j; ma++)
            k += ma*ma;                    
        for (int ma = 0; ma <= j; ma++) {
            int in = i + inum*k;                
            Utotr[in] += sfac*Pr[ma]; // atomic add   
            Utoti[in] += sfac*Pi[ma]; // atomic add                       
            k += 1;
        }                   
    }
    
    for (mb = 1; 2*mb <= twojmax; mb++) {     
        for (int ma = 0; ma < 2*mb; ma++) {                      
            Pr[ma] = Qr[ma];
            Pi[ma] = Qi[ma];
        }                
        for (int j = 2*mb; j <= twojmax; j++) { 
            int ma = 0;
            //double p_r, p_i;
            //double u_r = Pr[ma];
            //double u_i = Pi[ma];
            z = Pr[ma];
            z0 = Pi[ma];
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * z + a_i * z0);
            Pi[ma] = rootpq * (a_r * z0 - a_i * z);            
            for (ma = 1; ma < j; ma++) {
                x = Pr[ma];
                y = Pi[ma];
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                rcutij = rootpqarray[ma*jdim + (j - mb)];
                Pr[ma] = rootpq * (a_r * x + a_i * y) -rcutij * (b_r * z + b_i * z0);
                Pi[ma] = rootpq * (a_r * y - a_i * x) -rcutij * (b_r * z0 - b_i * z);
                z = x;
                z0 = y;
            }
            ma = j;
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = -rcutij * (b_r * z + b_i * z0);
            Pi[ma] = -rcutij * (b_r * z0 - b_i * z);       
            
            if (j==(2*mb)) {
                int mapar = 1;
                for (int ma = 0; ma <= j/2; ma++) {
                    if (mapar == 1) {                    
                        Pr[j/2+ma] = Pr[j/2-ma];
                        Pi[j/2+ma] = -Pi[j/2-ma];
                    } else {
                        Pr[j/2+ma] = -Pr[j/2-ma];
                        Pi[j/2+ma] = Pi[j/2-ma];
                    }
                    mapar = -mapar;        
                }                                                        
            }
            
            // store Qr, Qi, for the next mb level
            if (j==(2*mb+1)) {
                int mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
                for (int ma = 0; ma <= j; ma++) {     
                    if (mapar == 1) {                    
                        Qr[j-ma] = Pr[ma];
                        Qi[j-ma] = -Pi[ma];
                    } else {
                        Qr[j-ma] = -Pr[ma];
                        Qi[j-ma] =  Pi[ma];
                    }
                    mapar = -mapar;
                }                                                
            }
            
            int k =  1 + (j+1)*mb;
            for (int ma = 2; ma <= j; ma++)
                k += ma*ma;            
            for (int ma = 0; ma <= j; ma++) {
                int in = i + inum*k;                
                Utotr[in] += sfac*Pr[ma]; // atomic add   
                Utoti[in] += sfac*Pi[ma]; // atomic add                       
                k += 1; 
            }                                                           
        }
    }        
  }
};

void cpuAddWself2Ui(double *Utotr, double *Utoti, double wself, int *idxu_block, int *type, int *map, 
        int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum)
{
    int N1 = inum;
    int N2 = N1*(twojmax+1);
    int N3 = N2*nelements;                                
    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;  // inum*(twojmax+1)
        int ii = l%N1;    // inum
        int j = (l-ii)/N1; // (twojmax+1)
        int jelem = (idx-l)/N2; // nelements   
        int ielem = (chemflag) ? map[type[ii]]: 0;                
        int nmax = ii + inum*idxu_max*jelem;
        
        int jju = idxu_block[j];                
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {                
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)                        
                        Utotr[inum*jju + nmax] += wself;                                     
                jju++;                                
            }
        }
        
        // copy left side to right side with inversion symmetry VMK 4.4(2)
        // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
        
        jju = idxu_block[j];        
        int jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb < j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju =  inum*jju + nmax;
                int njjup = inum*jjup + nmax;
                if (mapar == 1) {
                    Utotr[njjup] = Utotr[njju];
                    Utoti[njjup] = -Utoti[njju];
                } else {
                    Utotr[njjup] = -Utotr[njju];
                    Utoti[njjup] =  Utoti[njju];
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }                    
};

void cpuSnapComputeEi(double *eatom, double *Utotr, double *Utoti, double *cglist, double *bzero, 
        double *coeffelem, int *ilist, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, 
        int twojmax, int idxb_max, int idxu_max, int nelements, int ncoeffall, int bnorm_flag, 
        int bzero_flag, int wselfall_flag, int inum)
{    
    int nelemsq = nelements*nelements;    
    int nu_max = idxu_max*inum;
    //int nb_max = idxb_max*inum;    
    int jdim = twojmax+1;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        //int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  
        const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];

        int itype = type[ii]; // element type of atom i
        int ielem = map[itype];  // index of that element type        
        int icoeff = jjb + idxb_max*(jelem);
        int itriple = 1+icoeff+ielem*ncoeffall;
        
        int jju = idxu_block[j];
        int ii1 = ii + inum*idxu_max*elem1;
        int ii2 = ii + inum*idxu_max*elem2;        
        int idu, jju1, jju2, icga, icgb;        
        int ia, ib, ma, mb, ma1min, ma2max, na, mb1min, mb2max, nb, ma1, ma2;
        double ztmp_r, ztmp_i, sumzu = 0.0;
        for (mb = 0; 2 * mb < j; mb++)
            for (ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                // calculate Z_{j1,j2,j}^{elem1,elem2}(ma,mb)   
                ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];                                              
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            mb = j / 2;
            for (ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;                
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];                                
                jju++;
            }
            ma = mb;
            idu = ii + inum*jju + nu_max*elem3;        
            ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
            mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
            jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            icgb = mb1min * (j2 + 1) + mb2max;
            ztmp_r = 0.0;
            ztmp_i = 0.0;
            for (ib = 0; ib < nb; ib++) {
                ma1 = ma1min;
                ma2 = ma2max;
                icga = ma1min * (j2 + 1) + ma2max;
                for (ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia
                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib
            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
            sumzu += 0.5*(Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz]);            
        } // end if jeven
        
        sumzu *= 2.0;                
        
        if (bzero_flag) {
          if (!wselfall_flag) {
            if (elem1 == elem2 && elem1 == elem3) {
              sumzu -= bzero[j];
            }
          } 
          else {
            sumzu -= bzero[j];
          }
        }
        
        //blist[idx] = sumzu;                                
        
        if (icoeff==0)
            eatom[ilist[ii]] += coeffelem[ielem*ncoeffall];
        eatom[ilist[ii]] += coeffelem[itriple]*sumzu;                           
    }
}

void cpuSnapComputeBi(double *blist, double *Utotr, double *Utoti, double *cglist, double *bzero, 
        int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int bzero_flag, int wselfall_flag, int inum)
{    
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    //int nb_max = idxb_max*inum;    
    int jdim = twojmax+1;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  
        const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int ii1 = ii + inum*idxu_max*elem1;
        int ii2 = ii + inum*idxu_max*elem2;        
        int idu, idz, jju1, jju2, icga, icgb;        
        int ia, ib, ma, mb, ma1min, ma2max, na, mb1min, mb2max, nb, ma1, ma2;
        double ztmp_r, ztmp_i, sumzu = 0.0;
        for (mb = 0; 2 * mb < j; mb++)
            for (ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                // calculate Z_{j1,j2,j}^{elem1,elem2}(ma,mb)   
                ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            mb = j / 2;
            for (ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
                ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
                na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
                mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
                mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
                nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
                jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
                jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
                icgb = mb1min * (j2 + 1) + mb2max;
                ztmp_r = 0.0;
                ztmp_i = 0.0;
                for (ib = 0; ib < nb; ib++) {
                    ma1 = ma1min;
                    ma2 = ma2max;
                    icga = ma1min * (j2 + 1) + ma2max;
                    for (ia = 0; ia < na; ia++) {
                        ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                        ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                      ma1++;
                      ma2--;
                      icga += j2;
                    } // end loop over ia
                    jju1 += j1 + 1;
                    jju2 -= j2 + 1;
                    icgb += j2;
                } // end loop over ib
                if (bnorm_flag) {
                  ztmp_i /= j+1;
                  ztmp_r /= j+1;
                }            
                sumzu += Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i;                
                //sumzu += Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            ma = mb;
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            ma1min = MAX(0, (2 * ma - j - j2 + j1) / 2);
            ma2max = (2 * ma - j - (2 * ma1min - j1) + j2) / 2;
            na = MIN(j1, (2 * ma - j + j2 + j1) / 2) - ma1min + 1;
            mb1min = MAX(0, (2 * mb - j - j2 + j1) / 2);
            mb2max = (2 * mb - j - (2 * mb1min - j1) + j2) / 2;
            nb = MIN(j1, (2 * mb - j + j2 + j1) / 2) - mb1min + 1;                
            jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            icgb = mb1min * (j2 + 1) + mb2max;
            ztmp_r = 0.0;
            ztmp_i = 0.0;
            for (ib = 0; ib < nb; ib++) {
                ma1 = ma1min;
                ma2 = ma2max;
                icga = ma1min * (j2 + 1) + ma2max;
                for (ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia
                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib
            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
            sumzu += 0.5*(Utotr[idu] * ztmp_r + Utoti[idu] * ztmp_i);            
            //sumzu += 0.5 * (Utotr[idu] * zlist_r[idz] + Utoti[idu] * zlist_i[idz]);
        } // end if jeven
        
        sumzu *= 2.0;                
        
        if (bzero_flag) {
          if (!wselfall_flag) {
            if (elem1 == elem2 && elem1 == elem3) {
              sumzu -= bzero[j];
            }
          } 
          else {
            sumzu -= bzero[j];
          }
        }
        
        blist[idx] = sumzu;                                        
    }
};

void cpuSnapComputeYi(double *ylist_r, double *ylist_i, double *Utotr, double *Utoti, double *cglist, double* beta, 
        int *map, int *type, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, int inum)
{    
    
    int N1 = idxu_max*nelements*inum;
    for (int i = 0; i< N1; i++) {
        ylist_r[i] = 0.0;
        ylist_i[i] = 0.0;
    }
    //cpuArraySetValue(ylist_r, (double) 0.0, N1);
    //cpuArraySetValue(ylist_i, (double) 0.0, N1);
    //printf("%i %i %i %i %i %i\n", inum, idxu_max, nelements, N1, idxz_max, idxb_max);
    
    int jdim = twojmax + 1;         
    int N2 = idxz_max*inum;                          
    for (int idx=0; idx < N2; idx++) {
      int ii = idx%inum;              
      int jjz = (idx-ii)/inum;         
      int jjz10 = jjz*10;
      const int j1 = idxz[jjz10+0];
      const int j2 = idxz[jjz10+1];
      const int j = idxz[jjz10+2];
      const int ma1min = idxz[jjz10+3];
      const int ma2max = idxz[jjz10+4];
      const int na = idxz[jjz10+5];
      const int mb1min = idxz[jjz10+6];
      const int mb2max = idxz[jjz10+7];
      const int nb = idxz[jjz10+8];          
      const int jju = idxz[jjz10+9];
      const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
      
      int itype = type[ii]; // element type of atom i
      int ielem = map[itype];  // index of that element type                        

      for(int elem1 = 0; elem1 < nelements; elem1++)
        for (int elem2 = 0; elem2 < nelements; elem2++) {        
          
            double ztmp_r = 0.0;
            double ztmp_i = 0.0;
            
            int ii1 = ii + inum*idxu_max*elem1;
            int ii2 = ii + inum*idxu_max*elem2;
            int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
            int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
            int icgb = mb1min * (j2 + 1) + mb2max;
            for (int ib = 0; ib < nb; ib++) {
                int ma1 = ma1min;
                int ma2 = ma2max;
                int icga = ma1min * (j2 + 1) + ma2max;
                for (int ia = 0; ia < na; ia++) {
                    ztmp_r += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)] - Utoti[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)]);
                    ztmp_i += cgblock[icgb]*cgblock[icga] * (Utotr[ii1+inum*(jju1+ma1)] * Utoti[ii2+inum*(jju2+ma2)] + Utoti[ii1+inum*(jju1+ma1)] * Utotr[ii2+inum*(jju2+ma2)]);                    
                  ma1++;
                  ma2--;
                  icga += j2;
                } // end loop over ia

                jju1 += j1 + 1;
                jju2 -= j2 + 1;
                icgb += j2;
            } // end loop over ib

            if (bnorm_flag) {
              ztmp_i /= j+1;
              ztmp_r /= j+1;
            }            
                               
            for(int elem3 = 0; elem3 < nelements; elem3++) {
              int itriple;  
              double betaj;
              if (j >= j1) {
                const int jjb = idxb_block[j + j2*jdim + j1*jdim*jdim]; 
                itriple = ((elem1 * nelements + elem2) * nelements + elem3) * idxb_max + jjb + 1 + ielem*ncoeffall;
                if (j1 == j) {
                  if (j2 == j) betaj = 3*beta[itriple];
                  else betaj = 2*beta[itriple];
                } else betaj = beta[itriple];          
              } else if (j >= j2) {
                const int jjb = idxb_block[j1 + j2*jdim + j*jdim*jdim];
                itriple = ((elem3 * nelements + elem2) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
                if (j2 == j) betaj = 2*beta[itriple];
                else betaj = beta[itriple];
              } else {
                const int jjb = idxb_block[j1 + j*jdim + j2*jdim*jdim];
                itriple = ((elem2 * nelements + elem3) * nelements + elem1) * idxb_max + jjb + 1 + ielem*ncoeffall;
                betaj = beta[itriple];
              }
              
              if (!bnorm_flag && j1 > j)
                betaj *= (j1 + 1) / (j + 1.0);
                         
              //printf("%i %i %i %i %i %i\n", inum, ii, jju, elem1, elem2, elem3);
              ylist_r[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_r;
              ylist_i[ii + inum*jju + inum*idxu_max*elem3] += betaj * ztmp_i;        
           }
        }         
    }  
};

void cpuSnapComputeFi(double *fatom, double *vatom, double *ylist_r, double *ylist_i, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem,  double rmin0, double rfac0, double rcutfac, int *map, int *aii, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int anum, int ijnum, int switchflag, int chemflag) 
{                 
  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];    
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);
    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double z0 = r / tan((r - rmin0) * rfac0 * M_PI / (rcutij - rmin0));                
    double dz0dr = z0 / r - (r*rfac0 * M_PI / (rcutij - rmin0)) * (rsq + z0 * z0) / rsq;

    double sfac = 0.0, dsfac = 0.0;        
    if (switchflag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switchflag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else if(r > rcutij) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else {
            double rcutfac = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
            dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
        }
    } 
    sfac *= wjelem[tj[ij]];
    dsfac *= wjelem[tj[ij]];
    
    double a_r, a_i, b_r, b_i, rootpq;        
    rcutij = 1.0 / sqrt(r * r + z0 * z0);
    a_r = rcutij * z0;
    a_i = -rcutij * z;
    b_r = rcutij * y;
    b_i = -rcutij * x;

    double u_r, u_i, ux_r, ux_i, uy_r, uy_i, uz_r, uz_i;
    double w_r, w_i, wx_r, wx_i, wy_r, wy_i, wz_r, wz_i;
    u_r = -pow(rcutij, 3.0) * (r + z0 * dz0dr);
    wx_r = u_r * ux;
    wy_r = u_r * uy;
    wz_r = u_r * uz;
    ux_r = dz0dr * ux;
    uy_r = dz0dr * uy;
    uz_r = dz0dr * uz;

    double dardx, daidx, dardy, daidy, dardz, daidz;
    dardx = ux_r * rcutij + z0 * wx_r;
    daidx = -z * wx_r;
    dardy = uy_r * rcutij + z0 * wy_r;
    daidy = -z * wy_r;
    dardz = uz_r * rcutij + z0 * wz_r;
    daidz = -z * wz_r;    
    daidz += -rcutij;

    double dbrdx, dbidx, dbrdy, dbidy, dbrdz, dbidz;
    dbrdx = y * wx_r;
    dbidx = -x * wx_r;    
    dbrdy = y * wy_r;
    dbidy = -x * wy_r;    
    dbrdz = y * wz_r;
    dbidz = -x * wz_r;        
    dbidx += -rcutij;
    dbrdy += rcutij;
    
    // 2Jmax = 10    
    double Pr[11], Pi[11], Qr[9], Qi[9];
    double Prx[11], Pix[11], Qrx[9], Qix[9];
    double Pry[11], Piy[11], Qry[9], Qiy[9];        
    double Prz[11], Piz[11], Qrz[9], Qiz[9];
    Pr[0] = 1.0; Pi[0] = 0.0;    
    Prx[0] = 0.0; Pix[0] = 0.0;    
    Pry[0] = 0.0; Piy[0] = 0.0;    
    Prz[0] = 0.0; Piz[0] = 0.0;        
    
    int jdim = twojmax + 1;
    int njelem = (chemflag==0) ? 0 : (inum*idxu_max*map[tj[ij]]);    
    int i = aii[ij] + njelem;                        
                    
    double dedx, dedy, dedz;       
    u_r = 0.5*ylist_r[i];  
    //e    =  sfac*u_r;
    dedx = (dsfac * ux) * u_r;
    dedy = (dsfac * uy) * u_r;
    dedz = (dsfac * uz) * u_r; 
            
    int j, k, ma, mb, mapar;    
    mb = 0;
    for (j = 1; j <= twojmax; j++) {        
        // fill in left side of matrix layer from previous layer
        ma = 0;
        u_r = Pr[ma];
        u_i = Pi[ma];
        ux_r = Prx[ma];
        ux_i = Pix[ma];            
        uy_r = Pry[ma];
        uy_i = Piy[ma];            
        uz_r = Prz[ma];
        uz_i = Piz[ma];                    
        rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
        Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
        Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);        
        Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
        Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
        Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
        Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
        Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
        Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                    
        for (ma = 1; ma < j; ma++) {
            w_r = Pr[ma];
            w_i = Pi[ma];
            wx_r = Prx[ma];
            wx_i = Pix[ma];            
            wy_r = Pry[ma];
            wy_i = Piy[ma];            
            wz_r = Prz[ma];
            wz_i = Piz[ma];                                        
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
            Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
            Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
            Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
            Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
            Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
            Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
            Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);            
            u_r = w_r;
            u_i = w_i;
            ux_r = wx_r;
            ux_i = wx_i;
            uy_r = wy_r;
            uy_i = wy_i;
            uz_r = wz_r;
            uz_i = wz_i;            
        }
        ma = j;
        rcutij = rootpqarray[ma*jdim + (j - mb)];
        Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
        Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);                        
        Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
        Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
        Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
        Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
        Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
        Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
                                
        if (j==(2*mb+1)) { // store Qr, Qi, for the next mb level
            mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
            for (ma = 0; ma <= j; ma++) {     
                if (mapar == 1) {                    
                    Qr[j-ma] = Pr[ma];
                    Qi[j-ma] = -Pi[ma];
                    Qrx[j-ma] =  Prx[ma];
                    Qix[j-ma] = -Pix[ma];
                    Qry[j-ma] =  Pry[ma];
                    Qiy[j-ma] = -Piy[ma];
                    Qrz[j-ma] =  Prz[ma];
                    Qiz[j-ma] = -Piz[ma];                    
                } else {
                    Qr[j-ma] = -Pr[ma];
                    Qi[j-ma] =  Pi[ma];
                    Qrx[j-ma] = -Prx[ma];
                    Qix[j-ma] =  Pix[ma];
                    Qry[j-ma] = -Pry[ma];
                    Qiy[j-ma] =  Piy[ma];
                    Qrz[j-ma] = -Prz[ma];
                    Qiz[j-ma] =  Piz[ma];                    
                }
                mapar = -mapar;
            }                              
        }
        
        k =  1 + (j+1)*mb;
        for (ma = 2; ma <= j; ma++)
            k += ma*ma;                    
        for (ma = 0; ma <= j; ma++) {                            
            rsq = ylist_r[i + inum*k]; 
            rinv = ylist_i[i + inum*k];
            //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
            dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
            dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
            dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;            
            k += 1;
        }                   
    }
    
    for (mb = 1; 2*mb <= twojmax; mb++) {     
        for (ma = 0; ma < 2*mb; ma++) {                      
            Pr[ma] = Qr[ma];
            Pi[ma] = Qi[ma];
            Prx[ma] = Qrx[ma];
            Pix[ma] = Qix[ma];
            Pry[ma] = Qry[ma];
            Piy[ma] = Qiy[ma];
            Prz[ma] = Qrz[ma];
            Piz[ma] = Qiz[ma];            
        }                
        for (j = 2*mb; j <= twojmax; j++) { 
            ma = 0;
            u_r = Pr[ma];
            u_i = Pi[ma];
            ux_r = Prx[ma];
            ux_i = Pix[ma];            
            uy_r = Pry[ma];
            uy_i = Piy[ma];            
            uz_r = Prz[ma];
            uz_i = Piz[ma];                                
            rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
            Pr[ma] = rootpq * (a_r * u_r + a_i * u_i);
            Pi[ma] = rootpq * (a_r * u_i - a_i * u_r);            
            Prx[ma] = rootpq * (dardx * u_r + daidx * u_i + a_r * ux_r + a_i * ux_i);
            Pix[ma] = rootpq * (dardx * u_i - daidx * u_r + a_r * ux_i - a_i * ux_r);
            Pry[ma] = rootpq * (dardy * u_r + daidy * u_i + a_r * uy_r + a_i * uy_i);
            Piy[ma] = rootpq * (dardy * u_i - daidy * u_r + a_r * uy_i - a_i * uy_r);
            Prz[ma] = rootpq * (dardz * u_r + daidz * u_i + a_r * uz_r + a_i * uz_i);
            Piz[ma] = rootpq * (dardz * u_i - daidz * u_r + a_r * uz_i - a_i * uz_r);                                
            for (ma = 1; ma < j; ma++) {
                w_r = Pr[ma];
                w_i = Pi[ma];
                wx_r = Prx[ma];
                wx_i = Pix[ma];            
                wy_r = Pry[ma];
                wy_i = Piy[ma];            
                wz_r = Prz[ma];
                wz_i = Piz[ma];                                            
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                rcutij = rootpqarray[ma*jdim + (j - mb)];
                Pr[ma] = rootpq * (a_r * w_r + a_i * w_i) -rcutij * (b_r * u_r + b_i * u_i);
                Pi[ma] = rootpq * (a_r * w_i - a_i * w_r) -rcutij * (b_r * u_i - b_i * u_r);
                Prx[ma] = rootpq * (dardx * w_r + daidx * w_i + a_r * wx_r + a_i * wx_i) -rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
                Pix[ma] = rootpq * (dardx * w_i - daidx * w_r + a_r * wx_i - a_i * wx_r) -rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
                Pry[ma] = rootpq * (dardy * w_r + daidy * w_i + a_r * wy_r + a_i * wy_i) -rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
                Piy[ma] = rootpq * (dardy * w_i - daidy * w_r + a_r * wy_i - a_i * wy_r) -rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
                Prz[ma] = rootpq * (dardz * w_r + daidz * w_i + a_r * wz_r + a_i * wz_i) -rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
                Piz[ma] = rootpq * (dardz * w_i - daidz * w_r + a_r * wz_i - a_i * wz_r) -rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);                            
                u_r = w_r;
                u_i = w_i;
                ux_r = wx_r;
                ux_i = wx_i;
                uy_r = wy_r;
                uy_i = wy_i;
                uz_r = wz_r;
                uz_i = wz_i;                
            }
            ma = j;
            rcutij = rootpqarray[ma*jdim + (j - mb)];
            Pr[ma] = -rcutij * (b_r * u_r + b_i * u_i);
            Pi[ma] = -rcutij * (b_r * u_i - b_i * u_r);       
            Prx[ma] =-rcutij * (dbrdx * u_r + dbidx * u_i + b_r * ux_r + b_i * ux_i);
            Pix[ma] =-rcutij * (dbrdx * u_i - dbidx * u_r + b_r * ux_i - b_i * ux_r);
            Pry[ma] =-rcutij * (dbrdy * u_r + dbidy * u_i + b_r * uy_r + b_i * uy_i);
            Piy[ma] =-rcutij * (dbrdy * u_i - dbidy * u_r + b_r * uy_i - b_i * uy_r);
            Prz[ma] =-rcutij * (dbrdz * u_r + dbidz * u_i + b_r * uz_r + b_i * uz_i);
            Piz[ma] =-rcutij * (dbrdz * u_i - dbidz * u_r + b_r * uz_i - b_i * uz_r);
            
            if (j==(2*mb)) {
                mapar = 1;
                for (ma = 0; ma <= j/2; ma++) {
                    if (mapar == 1) {                    
                        Pr[j/2+ma] = Pr[j/2-ma];
                        Pi[j/2+ma] = -Pi[j/2-ma];
                        Prx[j/2+ma] = Prx[j/2-ma];
                        Pix[j/2+ma] = -Pix[j/2-ma];
                        Pry[j/2+ma] = Pry[j/2-ma];
                        Piy[j/2+ma] = -Piy[j/2-ma];
                        Prz[j/2+ma] = Prz[j/2-ma];
                        Piz[j/2+ma] = -Piz[j/2-ma];                        
                    } else {
                        Pr[j/2+ma] = -Pr[j/2-ma];
                        Pi[j/2+ma] = Pi[j/2-ma];
                        Prx[j/2+ma] = -Prx[j/2-ma];
                        Pix[j/2+ma] =  Pix[j/2-ma];
                        Pry[j/2+ma] = -Pry[j/2-ma];
                        Piy[j/2+ma] =  Piy[j/2-ma];
                        Prz[j/2+ma] = -Prz[j/2-ma];
                        Piz[j/2+ma] =  Piz[j/2-ma];                        
                    }
                    mapar = -mapar;        
                }                                                        
            }
            
            if (j==(2*mb)) {
                mapar = 1;
                for (ma = 0; ma <= j; ma++) {
                    if (mapar == 1) {                    
                        Prx[j/2+ma] = Prx[j/2-ma];
                        Pix[j/2+ma] = -Pix[j/2-ma];
                        Pry[j/2+ma] = Pry[j/2-ma];
                        Piy[j/2+ma] = -Piy[j/2-ma];
                        Prz[j/2+ma] = Prz[j/2-ma];
                        Piz[j/2+ma] = -Piz[j/2-ma];                        
                    } else {
                        Prx[j/2+ma] = -Prx[j/2-ma];
                        Pix[j/2+ma] =  Pix[j/2-ma];
                        Pry[j/2+ma] = -Pry[j/2-ma];
                        Piy[j/2+ma] =  Piy[j/2-ma];
                        Prz[j/2+ma] = -Prz[j/2-ma];
                        Piz[j/2+ma] =  Piz[j/2-ma];                        
                    }
                    mapar = -mapar;        
                }                                                        
            }
                        
            // store Qr, Qi, for the next mb level
            if (j==(2*mb+1)) {
                mapar = (((j+1)/2)%2 == 0) ? -1 : 1;
                for (ma = 0; ma <= j; ma++) {     
                    if (mapar == 1) {                    
                        Qr[j-ma] = Pr[ma];
                        Qi[j-ma] = -Pi[ma];  
                        Qrx[j-ma] =  Prx[ma];
                        Qix[j-ma] = -Pix[ma];
                        Qry[j-ma] =  Pry[ma];
                        Qiy[j-ma] = -Piy[ma];
                        Qrz[j-ma] =  Prz[ma];
                        Qiz[j-ma] = -Piz[ma];                                            
                    } else {
                        Qr[j-ma] = -Pr[ma];
                        Qi[j-ma] =  Pi[ma];
                        Qrx[j-ma] = -Prx[ma];
                        Qix[j-ma] =  Pix[ma];
                        Qry[j-ma] = -Pry[ma];
                        Qiy[j-ma] =  Piy[ma];
                        Qrz[j-ma] = -Prz[ma];
                        Qiz[j-ma] =  Piz[ma];                                            
                    }
                    mapar = -mapar;
                }                                                
            }
            
            k =  1 + (j+1)*mb;
            for (ma = 2; ma <= j; ma++)
                k += ma*ma;                            
            if (j==(2*mb)) {
                for (ma = 0; ma < mb; ma++) {
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];         
                    //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }     
                ma = mb;
                rsq = 0.5*ylist_r[i + inum*k]; 
                rinv = 0.5*ylist_i[i + inum*k];        
                //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
            }
            else {
                for (ma = 0; ma <= j; ma++) {
                    rsq = ylist_r[i + inum*k]; 
                    rinv = ylist_i[i + inum*k];             
                    //e += sfac*(Pr[ma] * rsq + Pi[ma] * rinv);
                    dedx += (dsfac * Pr[ma] * ux + sfac * Prx[ma]) * rsq + (dsfac * Pi[ma] * ux + sfac * Pix[ma]) * rinv;
                    dedy += (dsfac * Pr[ma] * uy + sfac * Pry[ma]) * rsq + (dsfac * Pi[ma] * uy + sfac * Piy[ma]) * rinv;
                    dedz += (dsfac * Pr[ma] * uz + sfac * Prz[ma]) * rsq + (dsfac * Pi[ma] * uz + sfac * Piz[ma]) * rinv;
                    k += 1; 
                }                
            }            
        }
    }      

    ma = ai[ij];        
    mb = aj[ij];                    
    //eatom[ma] += 2.0*e;
    
    dedx = 2.0*dedx;      
    dedy = 2.0*dedy;      
    dedz = 2.0*dedz;             
    fatom[0+3*ma] += dedx;
    fatom[1+3*ma] += dedy;
    fatom[2+3*ma] += dedz;
    fatom[0+3*mb] -= dedx;
    fatom[1+3*mb] -= dedy;
    fatom[2+3*mb] -= dedz;    
    
    rootpq = -0.5;
    dbrdx = rootpq*x*dedx;
    dbidx = rootpq*y*dedy;
    dbrdy = rootpq*z*dedz;
    dbidy = rootpq*x*dedy;
    dbrdz = rootpq*x*dedz;
    dbidz = rootpq*y*dedz;        
    vatom[0*anum+ma] += dbrdx;
    vatom[1*anum+ma] += dbidx;
    vatom[2*anum+ma] += dbrdy;
    vatom[3*anum+ma] += dbidy;
    vatom[4*anum+ma] += dbrdz;
    vatom[5*anum+ma] += dbidz;        
    vatom[0*anum+mb] += dbrdx;
    vatom[1*anum+mb] += dbidx;
    vatom[2*anum+mb] += dbrdy;
    vatom[3*anum+mb] += dbidy;
    vatom[4*anum+mb] += dbrdz;
    vatom[5*anum+mb] += dbidz;            
  }   
};

void cpuComputeSij(double *Sr, double *Si, double *Srx, double *Six, double *Sry, double *Siy, double *Srz, double *Siz, 
        double *rootpqarray, double *rij, double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag)                
{    
  for(int ij=0; ij<ijnum; ij++) {        
    double x = rij[ij*3+0];
    double y = rij[ij*3+1];
    double z = rij[ij*3+2];
    double rsq = x * x + y * y + z * z;
    double r = sqrt(rsq);
    double rinv = 1.0 / r;
    double ux = x * rinv;
    double uy = y * rinv;
    double uz = z * rinv;

    double rcutij = (radelem[ti[ij]]+radelem[tj[ij]])*rcutfac; //(radelem[type[ii]]+radelem[type[jj]])*rcutfac;
    double rscale0 = rfac0 * M_PI / (rcutij - rmin0);
    double theta0 = (r - rmin0) * rscale0;
    double z0 = r / tan(theta0);                
    double dz0dr = z0 / r - (r*rscale0) * (rsq + z0 * z0) / rsq;
            
    double sfac = 0.0, dsfac = 0.0;        
    if (switch_flag == 0) {
        sfac = 1.0;
        dsfac = 0.0;
    }
    else if (switch_flag == 1) {
        if (r <= rmin0) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else if(r > rcutij) {
            sfac = 1.0;
            dsfac = 0.0;
        }
        else {
            double rcutfac = M_PI / (rcutij - rmin0);
            sfac =  0.5 * (cos((r - rmin0) * rcutfac) + 1.0);   
            dsfac = -0.5 * sin((r - rmin0) * rcutfac) * rcutfac;
        }
    } 
    sfac *= wjelem[tj[ij]];
    dsfac *= wjelem[tj[ij]];

    //sfac = 1.0; 
    //dsfac = 0.0;
    
    double r0inv, dr0invdr;
    double a_r, a_i, b_r, b_i;
    double da_r[3], da_i[3], db_r[3], db_i[3];
    double dz0[3], dr0inv[3];
    double rootpq;
    int jdim = twojmax + 1;
  
    r0inv = 1.0 / sqrt(r * r + z0 * z0);
    a_r = r0inv * z0;
    a_i = -r0inv * z;
    b_r = r0inv * y;
    b_i = -r0inv * x;

    dr0invdr = -pow(r0inv, 3.0) * (r + z0 * dz0dr);

    dr0inv[0] = dr0invdr * ux;
    dr0inv[1] = dr0invdr * uy;
    dr0inv[2] = dr0invdr * uz;

    dz0[0] = dz0dr * ux;
    dz0[1] = dz0dr * uy;
    dz0[2] = dz0dr * uz;

    for (int k = 0; k < 3; k++) {
        da_r[k] = dz0[k] * r0inv + z0 * dr0inv[k];
        da_i[k] = -z * dr0inv[k];
    }
    da_i[2] += -r0inv;

    for (int k = 0; k < 3; k++) {
        db_r[k] = y * dr0inv[k];
        db_i[k] = -x * dr0inv[k];
    }
    db_i[0] += -r0inv;
    db_r[1] += r0inv;
    
    Sr[ij+0*ijnum] = 1.0;
    Si[ij+0*ijnum] = 0.0;
    Srx[ij+0*ijnum] = 0.0;
    Six[ij+0*ijnum] = 0.0;
    Sry[ij+0*ijnum] = 0.0;
    Siy[ij+0*ijnum] = 0.0;
    Srz[ij+0*ijnum] = 0.0;
    Siz[ij+0*ijnum] = 0.0;
    for (int j = 1; j <= twojmax; j++) {
        int jju = idxu_block[j];
        int jjup = idxu_block[j-1];
        
        // fill in left side of matrix layer from previous layer
        for (int mb = 0; 2*mb <= j; mb++) {
            Sr[ij+jju*ijnum] = 0.0;
            Si[ij+jju*ijnum] = 0.0;
            Srx[ij+jju*ijnum] = 0.0;
            Six[ij+jju*ijnum] = 0.0;
            Sry[ij+jju*ijnum] = 0.0;
            Siy[ij+jju*ijnum] = 0.0;
            Srz[ij+jju*ijnum] = 0.0;
            Siz[ij+jju*ijnum] = 0.0;
            for (int ma = 0; ma < j; ma++) {
                rootpq = rootpqarray[(j - ma)*jdim + (j - mb)];
                int njju = ij+jju*ijnum;
                int njju1 = ij+(jju+1)*ijnum;
                int njjup = ij+jjup*ijnum;
                double u_r = Sr[njjup];
                double u_i = Si[njjup];
                double ux_r = Srx[njjup];
                double ux_i = Six[njjup];
                double uy_r = Sry[njjup];
                double uy_i = Siy[njjup];
                double uz_r = Srz[njjup];
                double uz_i = Siz[njjup];

                Sr[njju] += rootpq * (a_r * u_r + a_i * u_i);
                Si[njju] += rootpq * (a_r * u_i - a_i * u_r);
                Srx[njju] += rootpq * (da_r[0] * u_r + da_i[0] * u_i + a_r * ux_r + a_i * ux_i);
                Six[njju] += rootpq * (da_r[0] * u_i - da_i[0] * u_r + a_r * ux_i - a_i * ux_r);
                Sry[njju] += rootpq * (da_r[1] * u_r + da_i[1] * u_i + a_r * uy_r + a_i * uy_i);
                Siy[njju] += rootpq * (da_r[1] * u_i - da_i[1] * u_r + a_r * uy_i - a_i * uy_r);
                Srz[njju] += rootpq * (da_r[2] * u_r + da_i[2] * u_i + a_r * uz_r + a_i * uz_i);
                Siz[njju] += rootpq * (da_r[2] * u_i - da_i[2] * u_r + a_r * uz_i - a_i * uz_r);

                rootpq = rootpqarray[(ma + 1)*jdim + (j - mb)];
                Sr[njju1] = -rootpq * (b_r * u_r + b_i * u_i);
                Si[njju1] = -rootpq * (b_r * u_i - b_i * u_r);
                Srx[njju1] = -rootpq * (db_r[0] * u_r + db_i[0] * u_i + b_r * ux_r + b_i * ux_i);
                Six[njju1] = -rootpq * (db_r[0] * u_i - db_i[0] * u_r + b_r * ux_i - b_i * ux_r);
                Sry[njju1] = -rootpq * (db_r[1] * u_r + db_i[1] * u_i + b_r * uy_r + b_i * uy_i);
                Siy[njju1] = -rootpq * (db_r[1] * u_i - db_i[1] * u_r + b_r * uy_i - b_i * uy_r);
                Srz[njju1] = -rootpq * (db_r[2] * u_r + db_i[2] * u_i + b_r * uz_r + b_i * uz_i);
                Siz[njju1] = -rootpq * (db_r[2] * u_i - db_i[2] * u_r + b_r * uz_i - b_i * uz_r);
                jju++;
                jjup++;
            }
            jju++;
        }

    // copy left side to right side with inversion symmetry VMK 4.4(2)
    // u[ma-j][mb-j] = (-1)^(ma-mb)*Conj([u[ma][mb])
                   
        jju = idxu_block[j];
        jjup = jju+(j+1)*(j+1)-1;
        int mbpar = 1;
        for (int mb = 0; 2*mb <= j; mb++) {
            int mapar = mbpar;
            for (int ma = 0; ma <= j; ma++) {
                int njju = ij+jju*ijnum;
                int njjup = ij+jjup*ijnum;
                if (mapar == 1) {
                    Sr[njjup] = Sr[njju];
                    Si[njjup] = -Si[njju];
                    if (j%2==1 && mb==(j/2)) {                    
                    Srx[njjup] =  Srx[njju];
                    Six[njjup] = -Six[njju];
                    Sry[njjup] =  Sry[njju];
                    Siy[njjup] = -Siy[njju];
                    Srz[njjup] =  Srz[njju];
                    Siz[njjup] = -Siz[njju];
                    }
                } else {
                    Sr[njjup] = -Sr[njju];
                    Si[njjup] =  Si[njju];
                    if (j%2==1 && mb==(j/2)) {
                    Srx[njjup] = -Srx[njju];
                    Six[njjup] =  Six[njju];
                    Sry[njjup] = -Sry[njju];
                    Siy[njjup] =  Siy[njju];
                    Srz[njjup] = -Srz[njju];
                    Siz[njjup] =  Siz[njju];                    
                    }
                }
                mapar = -mapar;
                jju++;
                jjup--;
            }
            mbpar = -mbpar;
        }        
    }        

    for (int j = 0; j <= twojmax; j++) {
        int jju = idxu_block[j];
        for (int mb = 0; 2*mb <= j; mb++)
          for (int ma = 0; ma <= j; ma++) {
            int ijk = ij+jju*ijnum;               
            Srx[ijk] = dsfac * Sr[ijk] * ux + sfac * Srx[ijk]; 
            Six[ijk] = dsfac * Si[ijk] * ux + sfac * Six[ijk]; 
            Sry[ijk] = dsfac * Sr[ijk] * uy + sfac * Sry[ijk]; 
            Siy[ijk] = dsfac * Si[ijk] * uy + sfac * Siy[ijk]; 
            Srz[ijk] = dsfac * Sr[ijk] * uz + sfac * Srz[ijk]; 
            Siz[ijk] = dsfac * Si[ijk] * uz + sfac * Siz[ijk];                  
            jju++;
          }
    }
    
    for (int k=0; k<idxu_max; k++) {
        int ijk = ij + ijnum*k;
        Sr[ijk] = sfac*Sr[ijk];
        Si[ijk] = sfac*Si[ijk];
    }            
  }
};

void cpuZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum)
{
    int N1 = inum;
    int N2 = N1*(twojmax+1);
    int N3 = N2*nelements;                                
    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;  // inum*(twojmax+1)
        int ii = l%N1;    // inum
        int j = (l-ii)/N1; // (twojmax+1)
        int jelem = (idx-l)/N2; // nelements   
        int ielem = (chemflag) ? map[type[ai[ii]]]: 0;                
        int jju = idxu_block[j];        
        for (int mb = 0; mb <= j; mb++) {
            for (int ma = 0; ma <= j; ma++) {
                int n = ii + inum*jju + inum*idxu_max*jelem;        
                Stotr[n] = 0.0;
                Stoti[n] = 0.0;
                if (jelem == ielem || wselfall_flag)
                    if (ma==mb)
                        Stotr[n] = wself; ///// double check this
                jju++;
            }
        }
        
    }                    
};

 void cpuKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *map, int *ai, int *tj, int inum, int ijnum, int N1, int N2, int chemflag)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max   
        int jelem = (chemflag) ? map[tj[ij]] : 0;     
        int i = ai[ij] + inum*jju + N1*jelem;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
 void cpuKernelAddUarraytot(double *Stotr, double *Stoti, double *Sr, double *Si, 
        int *ai, int inum, int ijnum, int N2)
{    
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;  // ijnum
        int jju = (idx-ij)/ijnum;    // idxu_max        
        int i = ai[ij] + inum*jju;                
        Stotr[i] += Sr[idx];
        Stoti[i] += Si[idx];                
    }
};
void cpuAddUarraytot(double *Stotr, double *Stoti, double *Sr, 
        double *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag)
{   
    int N1 = inum*idxu_max;    
    int N2 = ijnum*idxu_max;    
    if (chemflag==0) {
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, ai, inum, ijnum, N2);          
    } else
        cpuKernelAddUarraytot(Stotr, Stoti, Sr, Si, map, ai, tj, inum, ijnum, N1, N2, chemflag);  
};

void cpuComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum)
{
    int jdim = twojmax + 1;    
    int N1 = inum;    
    int N2 = N1*idxz_max;
    int N3 = N2*nelements*nelements;                                
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;   //  inum*idxz_max
        int ii = l%inum;    // inum
        int jjz = (l-ii)/inum; // idxz_max
        int ielem = (idx-l)/N2;  // nelements*nelements  
        int elem2 = ielem%nelements; // nelements
        int elem1 = (ielem-elem2)/nelements; // nelements
              
        const int j1 = idxz[jjz*10+0];
        const int j2 = idxz[jjz*10+1];
        const int j = idxz[jjz*10+2];
        const int ma1min = idxz[jjz*10+3];
        const int ma2max = idxz[jjz*10+4];
        const int na = idxz[jjz*10+5];
        const int mb1min = idxz[jjz*10+6];
        const int mb2max = idxz[jjz*10+7];
        const int nb = idxz[jjz*10+8];
        int jju1 = idxu_block[j1] + (j1 + 1) * mb1min;
        int jju2 = idxu_block[j2] + (j2 + 1) * mb2max;
        int icgb = mb1min * (j2 + 1) + mb2max;

        const double *cgblock = cglist + idxcg_block[j + j2*jdim + j1*jdim*jdim];
        double qr = 0.0;
        double qi = 0.0;          
        for (int ib = 0; ib < nb; ib++) {
            double suma1_r = 0.0;
            double suma1_i = 0.0;

            // Stotr: inum*idxu_max*nelements  
            const double *u1_r = &Stotr[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u1_i = &Stoti[ii + inum*jju1 + inum*idxu_max*elem1];
            const double *u2_r = &Stotr[ii + inum*jju2 + inum*idxu_max*elem2];
            const double *u2_i = &Stoti[ii + inum*jju2 + inum*idxu_max*elem2];

            int ma1 = ma1min;
            int ma2 = ma2max;
            int icga = ma1min * (j2 + 1) + ma2max;

            for (int ia = 0; ia < na; ia++) {
                suma1_r += cgblock[icga] * (u1_r[inum*ma1] * u2_r[inum*ma2] - u1_i[inum*ma1] * u2_i[inum*ma2]);
                suma1_i += cgblock[icga] * (u1_r[inum*ma1] * u2_i[inum*ma2] + u1_i[inum*ma1] * u2_r[inum*ma2]);
                ma1++;
                ma2--;
                icga += j2;
            } // end loop over ia

            qr += cgblock[icgb] * suma1_r;
            qi += cgblock[icgb] * suma1_i;

            jju1 += j1 + 1;
            jju2 -= j2 + 1;
            icgb += j2;
        } // end loop over ib
        
        if (bnorm_flag) {
            qr /= (j+1);
            qi /= (j+1);
        }        
        
        zlist_r[idx] = qr;
        zlist_i[idx] = qi;          
    }
};

void cpuKernelComputeBi1(double *blist, double *zlist_r, double *zlist_i, 
        double *Stotr, double *Stoti, int *idxb, int *idxu_block, int *idxz_block, int jdim,         
        int nelements, int nelemsq, int nz_max, int nu_max, int inum, int N2, int N3)
{    
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;        
        int ii = l%inum;
        int jjb = (l-ii)/inum;
        int jelem = (idx-l)/N2;                    
        int k = jelem%nelemsq;      
        int elem3 = k%nelements;
        int elem2 = (k-elem3)/nelements;
        int elem1 = (jelem-k)/nelemsq;    
        //int itriple = elem3 + nelements*elem2 + nelemsq*elem1;    
        int idouble = elem2 + nelements*elem1;  
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

        int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
        int jju = idxu_block[j];
        int idu;
        int idz;
        double sumzu = 0.0;
        for (int mb = 0; 2 * mb < j; mb++)
            for (int ma = 0; ma <= j; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            } // end loop over ma, mb

        // For j even, handle middle column
        if (j % 2 == 0) {
            int mb = j / 2;
            for (int ma = 0; ma < mb; ma++) {
                idu = ii + inum*jju + nu_max*elem3;        
                idz = ii + inum*jjz + nz_max*idouble;        
                sumzu += Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz];
                jjz++;
                jju++;
            }
            idu = ii + inum*jju + nu_max*elem3;        
            idz = ii + inum*jjz + nz_max*idouble;        
            sumzu += 0.5 * (Stotr[idu] * zlist_r[idz] + Stoti[idu] * zlist_i[idz]);
        } // end if jeven

        blist[idx] = 2.0 * sumzu;                
    }
}
void cpuKernelComputeBi2(double *blist, double *bzero,int *ilist, int *type,
       int *map, int *idxb, int nelements, int nb_max, int inum, int N2, int chemflag)
{        
    for (int idx=0; idx < N2; idx++) {        
        int ii = idx%inum;        
        int jjb = (idx-ii)/inum;    
        
        int ielem = (chemflag) ? map[type[ilist[ii]]]: 0;                
        int itriple = (ielem*nelements+ielem)*nelements+ielem;

        const int j = idxb[jjb*3 + 2];  
        blist[ii + inum*jjb + nb_max*itriple] -= bzero[j];                
    }
}
void cpuKernelComputeBi4(double *blist, double *bzero,
       int *idxb, int inum, int N2, int N3)
{        
    for (int idx=0; idx < N3; idx++) {
        int l = idx%N2;
        int ii = l%inum;        
        int jjb = (l-ii)/inum;    
        int j = idxb[jjb*3 + 2];  
        blist[idx] -= bzero[j];        
    }
}
void cpuComputeBi2(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum)
{                
    int nelemsq = nelements*nelements;
    int nz_max = idxz_max*inum;
    int nu_max = idxu_max*inum;
    int nb_max = idxb_max*inum;
    int N2 = inum*idxb_max;
    int N3 = N2*nelements*nelemsq;
    int jdim = twojmax+1;

    cpuKernelComputeBi1(blist, zlist_r, zlist_i, Stotr, Stoti, 
            idxb, idxu_block, idxz_block, jdim, nelements, nelemsq, nz_max, nu_max, inum, N2, N3);

    if (bzero_flag) {
        if (!wselfall_flag) {
            cpuKernelComputeBi2(blist, bzero, ilist, type, map, 
                    idxb, nelements, nb_max, inum, N2, chemflag);
        }
        else {
            cpuKernelComputeBi4(blist, bzero, idxb, inum, N2, N3);            
        }
    }
};

void cpuComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, 
        double *dulist_r, double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum)
{                    
    int nz_max = idxz_max*inum;
    int nb_max = idxb_max*ijnum;    
    int nu_max = idxu_max*ijnum;    
    int N2 = ijnum*idxb_max;
    int jdim = twojmax+1;

    for (int i=0; i<nb_max*3*nelements*nelements*nelements; i++)
        dblist[i] = 0.0;
    //cpuArraySetValue(dblist, (double) 0.0, nb_max*3*nelements*nelements*nelements);
        
    for (int idx=0; idx < N2; idx++) {
        int ij = idx%ijnum;              
        int jjb = (idx-ij)/ijnum;                              
        int elem3 = (chemflag) ? map[tj[ij]] : 0;//(chemflag) ? map[type[alist[aj[ij]]]] : 0;
        int i = ai[ij]; // atom i
        const int j1 = idxb[jjb*3 + 0];
        const int j2 = idxb[jjb*3 + 1];
        const int j = idxb[jjb*3 + 2];  

       // Sum terms Conj(dudr(j,ma,mb))*z(j1,j2,j,ma,mb)
        for(int elem1 = 0; elem1 < nelements; elem1++)
            for(int elem2 = 0; elem2 < nelements; elem2++) {

            int jjz = idxz_block[j + j2*jdim + j1*jdim*jdim];
            int jju = idxu_block[j];
            int idouble = elem1*nelements+elem2;
            int itriple = (elem1*nelements+elem2)*nelements+elem3;
            int nimax = nz_max*idouble;                      

            double *dbdr = &dblist[nb_max*3*itriple];
            double sumzdu_r[3];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j; mb++)
              for (int ma = 0; ma <= j; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j even, handle middle column

            if (j % 2 == 0) {
              int mb = j / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;
                jjz++;
                jju++;
              }

              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if jeven
                                
            for (int k = 0; k < 3; k++)
              dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
            
            // Sum over Conj(dudr(j1,ma1,mb1))*z(j,j2,j1,ma1,mb1)
            double j1fac = (j + 1) / (j1 + 1.0);
            idouble = elem1*nelements+elem2;
            itriple = (elem3*nelements+elem2)*nelements+elem1;            
            //jjz = idxz_block[j][j2][j1];
            jjz = idxz_block[j1 + j2*jdim + j*jdim*jdim];
            jju = idxu_block[j1];

            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j1; mb++)
              for (int ma = 0; ma <= j1; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j1 even, handle middle column

            if (j1 % 2 == 0) {
              int mb = j1 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                       
                jjz++;
                jju++;
              }
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j1even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j1fac;

            // Sum over Conj(dudr(j2,ma2,mb2))*z(j,j1,j2,ma2,mb2)
            double j2fac = (j + 1) / (j2 + 1.0);
            idouble = elem2*nelements+elem1;
            itriple = (elem1*nelements+elem3)*nelements+elem2;
            //jjz = idxz_block[j][j1][j2];
            jjz = idxz_block[j2 + j1*jdim + j*jdim*jdim];        
            jju = idxu_block[j2];
            
            //dbdr = &dblist[nb_max*3*itriple];
            for (int k = 0; k < 3; k++)
              sumzdu_r[k] = 0.0;

            for (int mb = 0; 2 * mb < j2; mb++)
              for (int ma = 0; ma <= j2; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                            
                jjz++;
                jju++;
              } //end loop over ma mb

            // For j2 even, handle middle column

            if (j2 % 2 == 0) {
              int mb = j2 / 2;
              for (int ma = 0; ma < mb; ma++) {
                int n1 = i + inum*jjz + nimax;
                int n2 = ij+ ijnum*jju;
                double z_r = zlist_r[n1];
                double z_i = zlist_i[n1];
                for (int k = 0; k < 3; k++)
                  sumzdu_r[k] += dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i;                                                                 
                jjz++;
                jju++;
              }
              
              int n1 = i + inum*jjz + nimax;
              int n2 = ij+ ijnum*jju;
              double z_r = zlist_r[n1];
              double z_i = zlist_i[n1];
              for (int k = 0; k < 3; k++)
                sumzdu_r[k] += (dulist_r[n2 + nu_max*k] * z_r + dulist_i[n2 + nu_max*k] * z_i) * 0.5;
              // jjz++;
              // jju++;
            } // end if j2even

            for (int k = 0; k < 3; k++)
              if (bnorm_flag)
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k];
              else
                dbdr[ij + ijnum*k + ijnum*3*jjb] += 2.0 * sumzdu_r[k] * j2fac;
          }        
    }
}

void cpuSnapTallyBispectrum(double *bi, double *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{      
    //cpuArraySetValue(bi, (double) 0.0, nperdim*ntype);
    for (int i=0; i<nperdim*ntype; i++)
        bi[i] = 0.0;
    
    int N2 = inum*ncoeff;
    for (int idx=0; idx<N2; idx++) {        
        int ii = idx%inum;
        int icoeff = (idx-ii)/inum;
        int i = ilist[ii]; // index of atom i
        int itype = type[i]; // element type of atom i        
        int n = nperdim*(itype-1);
        bi[icoeff+n] += bispectrum[ii + inum*icoeff];                  
        //printf("%i %i %i %i %i %i %i \n", idx, ii, icoeff, i, itype, n, quadraticflag);
        if (quadraticflag==1) {
            int k = n+ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            double bveci = bispectrum[ii + inum*icoeff];
            bi[k] += 0.5*bveci*bveci;
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                double bvecj = bispectrum[ii + inum*jcoeff];
                bi[k++] += bveci*bvecj;
            }
        }                
    }
}

void cpuSnapTallyBispectrumDeriv(double *db, double *bispectrum, double *dbdr, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{   
    //cpuArraySetValue(db, (double) 0.0, inum*3*nperdim*ntype);
    for (int i=0; i<inum*3*nperdim*ntype; i++)
        db[i] = 0.0;
        
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int ii = aii[ij]; // index of atom i
        int i = ai[ij]; // index of atom i
        int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = nperdim*(itype-1);        
        int nii = inum*3*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        //printf("%i %i %i %i %i %i %i %i %i %i %i \n", idx, ij, icoeff, ii, i, j, itype, n, nii, nij, quadraticflag);

        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];        
        db[0 + 3*i + nii] += bix; 
        db[1 + 3*i + nii] += biy;
        db[2 + 3*i + nii] += biz;
        db[0 + 3*j + nii] -= bix;
        db[1 + 3*j + nii] -= biy;
        db[2 + 3*j + nii] -= biz;
        
        if (quadraticflag) {
            double bi = bispectrum[ii + inum*icoeff];
            double dbxtmp = bi*bix;
            double dbytmp = bi*biy;
            double dbztmp = bi*biz;
            int k = ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            nii = inum*3*(k + n);  
            db[0 + 3*i + nii] += dbxtmp;
            db[1 + 3*i + nii] += dbytmp;
            db[2 + 3*i + nii] += dbztmp;
            db[0 + 3*j + nii] -= dbxtmp;
            db[1 + 3*j + nii] -= dbytmp;
            db[2 + 3*j + nii] -= dbztmp;          
         
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                int nj = ijnum*3*jcoeff;
                double bjx = dbdr[ij + ijnum*0 + nj];
                double bjy = dbdr[ij + ijnum*1 + nj];
                double bjz = dbdr[ij + ijnum*2 + nj];        
                double bj = bispectrum[ii + inum*jcoeff];                
                dbxtmp = bi*bjx + bix*bj;
                dbytmp = bi*bjy + biy*bj;
                dbztmp = bi*bjz + biz*bj;

                k += 1;
                nii = inum*3*(k + n);  
                db[0 + 3*i + nii] += dbxtmp;
                db[1 + 3*i + nii] += dbytmp;
                db[2 + 3*i + nii] += dbztmp;
                db[0 + 3*j + nii] -= dbxtmp;
                db[1 + 3*j + nii] -= dbytmp;
                db[2 + 3*j + nii] -= dbztmp;                          
//                 db[i + inum*0 + nii] += dbxtmp;
//                 db[i + inum*1 + nii] += dbytmp;
//                 db[i + inum*2 + nii] += dbztmp;
//                 db[j + inum*0 + nii] -= dbxtmp;
//                 db[j + inum*1 + nii] -= dbytmp;
//                 db[j + inum*2 + nii] -= dbztmp;                                          
            }            
        }        
    }
}

void cpuSnapTallyBispectrumVirial(double *bv, double *bispectrum, double *dbdr, double *rij, int *aii, 
        int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag)
{      
    //cpuArraySetValue(bv, (double) 0.0, 6*nperdim*ntype);
    for (int i=0; i<6*nperdim*ntype; i++)
        bv[i] = 0.0;
    
    int N2 = ijnum*ncoeff;
    for (int idx=0; idx<N2; idx++) {
        int ij = idx%ijnum;
        int icoeff = (idx-ij)/ijnum;        
        int ii = aii[ij]; // index of atom i
        //int i = ai[ij]; // index of atom i
        //int j = aj[ij]; // index of atom i
        int itype = ti[ij]; // element type of atom i       
        int n = nperdim*(itype-1);        
        int nii = 6*(icoeff + n);  
        int nij = ijnum*3*icoeff;
        
        double factor = 1.0;
        double dx = -rij[0+3*ij];
        double dy = -rij[1+3*ij];
        double dz = -rij[2+3*ij];                    
        double bix = dbdr[ij + ijnum*0 + nij];
        double biy = dbdr[ij + ijnum*1 + nij];
        double biz = dbdr[ij + ijnum*2 + nij];                
        double v0 = factor*dx*bix;
        double v1 = factor*dy*biy;
        double v2 = factor*dz*biz;
        double v3 = factor*dx*biy;
        double v4 = factor*dx*biz;
        double v5 = factor*dy*biz;        
        
        bv[0 + nii] += v0;
        bv[1 + nii] += v1;
        bv[2 + nii] += v2;
        bv[3 + nii] += v3;
        bv[4 + nii] += v4;
        bv[5 + nii] += v5;        
//         bv[0 + nii] += v0;
//         bv[1 + nii] += v1;
//         bv[2 + nii] += v2;
//         bv[3 + nii] += v3;
//         bv[4 + nii] += v4;
//         bv[5 + nii] += v5;                       
        if (quadraticflag) {
            double bi = bispectrum[ii + inum*icoeff];
            double dbxtmp = bi*bix;
            double dbytmp = bi*biy;
            double dbztmp = bi*biz;
            int k = ncoeff + ncoeff*(ncoeff+1)/2 - (ncoeff-icoeff)*(ncoeff-icoeff+1)/2;
            nii = 6*(k + n);  
            double v0 = factor*dx*dbxtmp;
            double v1 = factor*dy*dbytmp;
            double v2 = factor*dz*dbztmp;
            double v3 = factor*dx*dbytmp;
            double v4 = factor*dx*dbztmp;
            double v5 = factor*dy*dbztmp;        
            bv[0 + nii] += v0;
            bv[1 + nii] += v1;
            bv[2 + nii] += v2;
            bv[3 + nii] += v3;
            bv[4 + nii] += v4;
            bv[5 + nii] += v5;        
//             bv[0 + nii] += v0;
//             bv[1 + nii] += v1;
//             bv[2 + nii] += v2;
//             bv[3 + nii] += v3;
//             bv[4 + nii] += v4;
//             bv[5 + nii] += v5;                                   
            for (int jcoeff = icoeff+1; jcoeff < ncoeff; jcoeff++) {
                int nj = ijnum*3*jcoeff;
                double bjx = dbdr[ij + ijnum*0 + nj];
                double bjy = dbdr[ij + ijnum*1 + nj];
                double bjz = dbdr[ij + ijnum*2 + nj];        
                double bj = bispectrum[ii + inum*jcoeff];                
                dbxtmp = bi*bjx + bix*bj;
                dbytmp = bi*bjy + biy*bj;
                dbztmp = bi*bjz + biz*bj;

                k += 1;                
                nii = 6*(k + n);  
                double v0 = factor*dx*dbxtmp;
                double v1 = factor*dy*dbytmp;
                double v2 = factor*dz*dbztmp;
                double v3 = factor*dx*dbytmp;
                double v4 = factor*dx*dbztmp;
                double v5 = factor*dy*dbztmp;        
                bv[0 + nii] += v0;
                bv[1 + nii] += v1;
                bv[2 + nii] += v2;
                bv[3 + nii] += v3;
                bv[4 + nii] += v4;
                bv[5 + nii] += v5;        
//                 bv[0 + nii] += v0;
//                 bv[1 + nii] += v1;
//                 bv[2 + nii] += v2;
//                 bv[3 + nii] += v3;
//                 bv[4 + nii] += v4;
//                 bv[5 + nii] += v5;                                                   
            }            
        }        
    }
}


#endif

