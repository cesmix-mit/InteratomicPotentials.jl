#pragma once

#ifdef __cplusplus
extern "C" {
#endif

int idxcgcount(int twojmax);

int idxucount(int twojmax);

int idxbcount(int twojmax);

int idxzcount(int twojmax);

void cpuInitSna(double *rootpqarray, double *cglist, double *factorial, int *idx_max, int *idxz, 
      int *idxz_block, int *idxb, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax);

void cpuNeighPairList(int *pairnum, int *pairlist, double *x, double *rcutsq, int *ilist, int *neighlist, 
        int *neighnum, int *atomtype, int *alist, int inum,  int dim, int ntypes);

void cpuNeighPairs(double *xij, double *x, int *aii, int *ai, int *aj,  int *ti, int *tj, 
        int *pairlist, int *pairnumsum, int *ilist, int *atomtype, int *alist, int inum, int dim);

void cpuSnapComputeUi(double *Utotr, double *Utoti, double *rootpqarray, double *rij, double *wjelem, double *radelem, 
        double rmin0, double rfac0, double rcutfac, int *map, int *aii, int *ti, int *tj, 
        int twojmax, int idxu_max, int inum, int ijnum, int switchflag, int chemflag);

void cpuAddWself2Ui(double *Utotr, double *Utoti, double wself, int *idxu_block, int *type, int *map, 
        int wselfall_flag, int chemflag, int idxu_max, int nelements, int twojmax, int inum);

void cpuSnapComputeEi(double *eatom, double *Utotr, double *Utoti, double *cglist, double *bzero, 
        double *coeffelem, int *ilist, int *map, int *type, int *idxb, int *idxcg_block, int *idxu_block, 
        int twojmax, int idxb_max, int idxu_max, int nelements, int ncoeffall, int bnorm_flag, 
        int bzero_flag, int wselfall_flag, int inum);

void cpuSnapComputeBi(double *blist, double *Utotr, double *Utoti, double *cglist, double *bzero, 
        int *idxb, int *idxcg_block, int *idxu_block, int *idxz_block, int twojmax, int idxb_max, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int bzero_flag, int wselfall_flag, int inum);

void cpuSnapComputeYi(double *ylist_r, double *ylist_i, double *Utotr, double *Utoti, double *cglist, double* beta, 
        int *map, int *type, int *idxz, int *idxb_block, int *idxu_block, int *idxcg_block, int twojmax, 
        int idxb_max, int idxu_max, int idxz_max, int nelements, int ncoeffall, int bnorm_flag, int inum);

void cpuSnapComputeFi(double *fatom, double *vatom, double *ylist_r, double *ylist_i, double *rootpqarray, double *rij, 
        double *wjelem, double *radelem,  double rmin0, double rfac0, double rcutfac, int *map, int *aii, int *ai, int *aj, 
        int *ti, int *tj, int twojmax, int idxu_max, int inum, int anum, int ijnum, int switchflag, int chemflag); 

void cpuComputeSij(double *Sr, double *Si, double *Srx, double *Six, double *Sry, double *Siy, double *Srz, double *Siz, 
        double *rootpqarray, double *rij, double *wjelem, double *radelem, double rmin0, double rfac0, double rcutfac, int *idxu_block,  
        int *ti, int *tj, int twojmax, int idxu_max, int ijnum, int switch_flag);

void cpuZeroUarraytot2(double *Stotr, double *Stoti, double wself, int *idxu_block, 
        int *type, int *map, int *ai, int wselfall_flag, int chemflag, int idxu_max, int nelements, 
         int twojmax, int inum);

void cpuAddUarraytot(double *Stotr, double *Stoti, double *Sr, 
        double *Si, int *map, int *ai, int *tj, int idxu_max, int inum, int ijnum, int chemflag);

void cpuComputeZi2(double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *cglist, int *idxz, int *idxu_block, int *idxcg_block, int twojmax, int idxu_max, 
        int idxz_max, int nelements, int bnorm_flag, int inum);

void cpuComputeBi2(double *blist, double *zlist_r, double *zlist_i, double *Stotr, double *Stoti, 
        double *bzero, int *ilist, int *type, int *map, int *idxb, int *idxu_block, int *idxz_block, 
        int twojmax, int idxb_max, int idxu_max, int idxz_max, int nelements, int bzero_flag, 
        int wselfall_flag, int chemflag, int inum);

void cpuComputeDbidrj(double *dblist, double *zlist_r, double *zlist_i, 
        double *dulist_r, double *dulist_i, int *idxb, int *idxu_block, int *idxz_block, 
        int *map, int *ai, int *tj, int twojmax, int idxb_max, int idxu_max, int idxz_max, 
        int nelements, int bnorm_flag, int chemflag, int inum, int ijnum);

void cpuSnapTallyBispectrum(double *bi, double *bispectrum, int *ilist, 
        int *type, int inum, int ncoeff, int nperdim, int ntype, int quadraticflag);

void cpuSnapTallyBispectrumDeriv(double *db, double *bispectrum, double *dbdr, int *aii, 
        int *ai, int *aj, int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

void cpuSnapTallyBispectrumVirial(double *bv, double *bispectrum, double *dbdr, double *rij, int *aii, 
        int *ti, int inum, int ijnum, int ncoeff, int nperdim, int ntype, int quadraticflag);

#ifdef __cplusplus
}
#endif


