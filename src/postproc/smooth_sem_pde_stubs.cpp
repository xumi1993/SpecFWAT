#include <stdio.h>
// #ifdef WITH_MPI
#include <mpi.h>
// #endif
const int  NGLL =  5;
const int  NGLLX = 5;
const int  NGLL2 = 25;
const int  NGLL3 = 125;
const int NGLL3PAD = 128;

typedef float realw;
typedef const realw* __restrict__ realw_const_p;
typedef realw* __restrict__ realw_p;

extern "C" void
smooth_pde_cuda_(int *h_is_sph, int *h_nspec, int *h_nglob,int *nprocs,
                int *h_nstep,int *nspec_max,int *phase_ispec_inner_elastic,
                int *h_nspec_outer,int *h_nspec_inner,
                const realw* xix,const realw* xiy,const realw* xiz,
                const realw* etax,const realw* etay,const realw* etaz,
                const realw* gamx,const realw* gamy,const realw* gamz,
                const realw *jaco, const realw* rotate0, 
                const realw* wgllwgll_xy,const realw* wgllwgll_xz,
                const realw* wgllwgll_yz,const realw* hprimeT,
                const realw* hprime_wgll,const int *ibool,
                const int *is_CPML, const realw *h_cv,const realw *h_ch,
                int *h_num_intfs,int *h_max_nibool,const int *my_neighbors,
                const int *nibool_intf,const int* ibool_intf,
                const realw* dat_bak, realw_p dat,const realw *rvol,
                realw_p dat_glob,realw_p ddat_glob){}