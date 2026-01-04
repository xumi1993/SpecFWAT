#include <cuda_runtime.h>
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



__device__ __forceinline__ void 
mxm_optx(const realw *u, const realw *hprime,int I,int J, int K,
        realw_p temp1)
{
    *temp1 = 0.;
    for(int l = 0; l < NGLL; l ++) {
        *temp1 += u[K*NGLL2+J*NGLL+l] * hprime[I*NGLL+l];
    }
}

__device__ __forceinline__ void 
mxm_opty(const realw *u, const realw *hprime,int I,int J, int K,
        realw_p temp1)
{
    *temp1 = 0.;
    for(int l = 0; l < NGLL; l ++) {
        *temp1 += u[K*NGLL2+l*NGLL+I] * hprime[J*NGLL+l];
    }
}


__device__ __forceinline__ void 
mxm_optz(const realw *u, const realw *hprime,int I,int J, int K,
        realw_p temp1)
{
    *temp1 = 0.;
    for(int l = 0; l < NGLL; l ++) {
        *temp1 += u[l*NGLL2+J*NGLL+I] * hprime[K*NGLL+l];
    }
}

__global__ void 
kernel_project(int num_elmts,const int *ibool,const realw* dat,
              const realw *rvol_local,realw_p dat_glob)
{
    int ispec = blockIdx.x + gridDim.x*blockIdx.y;
    int tx = threadIdx.x;
    if(tx >=NGLL3) tx = NGLL3 - 1;
    if(ispec >= num_elmts ) return;
    int offset = ispec * NGLL3PAD + tx;
    int iglob = ibool[offset] - 1;
    if(threadIdx.x < NGLL3) {
        realw temp = dat[offset] * rvol_local[offset];
        atomicAdd(&dat_glob[iglob],temp);
    }
}


__global__ static void 
kernel_prepare_boundary_matrix (realw* array_val,int ndim,int num_interfaces,
                         int max_nibool_interfaces, const int *nibool_interfaces,
                          const int *ibool_interfaces,
                          realw * buffer_send_matrix,int recv_stage)
{
    int id = threadIdx.x + (blockIdx.x + blockIdx.y*gridDim.x)*blockDim.x;

    // shape is like (num_interface,max_nibool,ndim)
    int iinterface = id / (max_nibool_interfaces * ndim);
    id = id % (max_nibool_interfaces * ndim);
    int ipoin = id / ndim;
    int i = id % ndim;

    if(iinterface >= num_interfaces || i >= ndim || ipoin >=nibool_interfaces[iinterface]) return;

    // entry in interface array
    int ientry = ipoin + max_nibool_interfaces*iinterface;
    // global index in wavefield
    int iglob = ibool_interfaces[ientry] - 1;

    if(recv_stage == 0) {
        buffer_send_matrix[ndim*ientry + i] = array_val[ndim*iglob + i];
    }
    else {
        atomicAdd(&array_val[ndim*iglob + i],buffer_send_matrix[ndim*ientry + i]);
    }
}


__global__ void 
kernel_smooth_sph_pde(int is_sph, int num_elmts,int iphase,int num_phase_ispec,
                    const int* phase_ispec_inner_elastic,
                    const realw* xix,const realw* xiy,const realw* xiz,
                    const realw* etax,const realw* etay,const realw* etaz,
                    const realw* gamx,const realw* gamy,const realw* gamz,
                    const realw *jaco, const realw* rotate, 
                    const realw* wgllwgll_xy,const realw* wgllwgll_xz,
                    const realw* wgllwgll_yz,const realw* hprimeT,
                    const realw* hprime_wgll,const int *ibool,
                    const realw cv,const realw ch,
                    const realw* dat_glob,realw_p ddat_glob,
                    const int * irregular_element_number, realw xix_regular)
{
    int bx = blockIdx.x + gridDim.x*blockIdx.y;
    if(bx >= num_elmts) return;
    int ispec = phase_ispec_inner_elastic[bx + (iphase - 1) * num_phase_ispec] - 1;
    int tx = threadIdx.x;
    if(tx >= NGLL3) tx = NGLL3 - 1; 
    int offset = ispec * NGLL3PAD + tx;
    int iglob = ibool[offset] - 1;

    // ijk
    int K,J,I;
    K = (tx/NGLL2);
    J = ((tx-K*NGLL2)/NGLLX);
    I = (tx-K*NGLL2-J*NGLLX);

    // shared memory
    __shared__ realw sh_u[NGLL3];
    __shared__ realw sh_hprimeT[NGLL2], sh_hprimewgll[NGLL2];

    // copy shared memory
    if(threadIdx.x < NGLL3) {
        sh_u[tx] = dat_glob[iglob];
        if(threadIdx.x < NGLL2) {
            sh_hprimeT[tx] = hprimeT[tx];
            sh_hprimewgll[tx] = hprime_wgll[tx];
        }
    }
    __syncthreads();

    // compute gradient
    realw dudx{},dudy{},dudz{};
    realw temp1{},temp2{},temp3{};
    realw rxl,ryl,rzl,xixl,xiyl,xizl,etaxl,etayl,
          etazl,gamxl,gamyl,gamzl,jacobianl;

    // rotate with shape(nspec)
    rxl = rotate[(ispec*3+0)*NGLL3PAD + tx];
    ryl = rotate[(ispec*3+1)*NGLL3PAD + tx];
    rzl = rotate[(ispec*3+2)*NGLL3PAD + tx];

    int ispec_irreg = irregular_element_number[ispec];
    if (ispec_irreg != 0) {
        int offset_irreg = (ispec_irreg - 1) * NGLL3PAD + tx;
        xixl = xix[offset_irreg];
        xiyl = xiy[offset_irreg];
        xizl = xiz[offset_irreg];
        etaxl = etax[offset_irreg];
        etayl = etay[offset_irreg];
        etazl = etaz[offset_irreg];
        gamxl = gamx[offset_irreg];
        gamyl = gamy[offset_irreg];
        gamzl = gamz[offset_irreg];
    } else {
        xixl = xix_regular;
        xiyl = xix_regular;
        xizl = xix_regular;
        etaxl = xix_regular;
        etayl = xix_regular;
        etazl = xix_regular;
        gamxl = xix_regular;
        gamyl = xix_regular;
        gamzl = xix_regular;
    }
    jacobianl = jaco[offset];
    mxm_optx(sh_u,sh_hprimeT,I,J,K,&temp1);
    mxm_opty(sh_u,sh_hprimeT,I,J,K,&temp2);
    mxm_optz(sh_u,sh_hprimeT,I,J,K,&temp3);
    dudx = temp1 * xixl + temp2 * etaxl + temp3 * gamxl;
    dudy = temp1 * xiyl + temp2 * etayl + temp3 * gamyl;
    dudz = temp1 * xizl + temp2 * etazl + temp3 * gamzl;

    // compute new terms
    __syncthreads();
    if(threadIdx.x < NGLL3) {
        if (is_sph){
            sh_u[tx] = ((cv-ch) * (rxl*xixl+ryl*xiyl+rzl*xizl) * 
                        (rxl*dudx+ryl*dudy+rzl*dudz) +
                        ch * (xixl*dudx+xiyl*dudy + 
                        xizl*dudz)) * jacobianl;
        } else {
            sh_u[tx] = ((cv-ch)*xizl*dudz +
                        ch*(xixl*dudx + xiyl*dudy +
                        xizl*dudz)) * jacobianl;
        }
    }
    __syncthreads();
    mxm_optx(sh_u,sh_hprimewgll,I,J,K,&temp1);

    __syncthreads();
    if(threadIdx.x < NGLL3) {
        if (is_sph){
            sh_u[tx] = ((cv-ch) * (rxl*etaxl+ryl*etayl+rzl*etazl) * 
                        (rxl*dudx+ryl*dudy+rzl*dudz) +
                        ch * (etaxl*dudx+etayl*dudy +
                        etazl*dudz)) * jacobianl;
        } else {
            sh_u[tx] = ((cv-ch)*etazl*dudz +
                        ch*(etaxl*dudx + etayl*dudy +
                        etazl*dudz)) * jacobianl;
        }
        
    }
    __syncthreads();
    mxm_opty(sh_u,sh_hprimewgll,I,J,K,&temp2);

    __syncthreads();
    if(threadIdx.x < NGLL3) {
        if (is_sph) {
            sh_u[tx] = ((cv-ch) * (rxl*gamxl+ryl*gamyl+rzl*gamzl) * 
                        (rxl*dudx+ryl*dudy+rzl*dudz) +
                        ch * (gamxl*dudx+gamyl*dudy + 
                        gamzl*dudz)) * jacobianl;
        } else {
            sh_u[tx] = ((cv-ch)*gamzl*dudz +
                        ch*(gamxl*dudx + gamyl*dudy +
                        gamzl*dudz)) * jacobianl;
        }
        
    }
    __syncthreads();
    mxm_optz(sh_u,sh_hprimewgll,I,J,K,&temp3);

    // update ddat_glob
    realw fac1 = wgllwgll_yz[K*NGLL + J];
    realw fac2 = wgllwgll_xz[K*NGLL + I];
    realw fac3 = wgllwgll_xy[J*NGLL + I];
    realw sum_term = fac1*temp1+ fac2*temp2 + fac3*temp3;

    if(threadIdx.x < NGLL3)
        atomicAdd(&ddat_glob[iglob],-sum_term);
}

__global__ void 
apply_mass(int nglob,const realw * rvol, realw_p dat_glob,
            realw_p ddat_glob)
{
    int id = threadIdx.x + blockIdx.x * blockDim.x; 
    if(id < nglob) {
        ddat_glob[id] *= rvol[id];
        dat_glob[id] += ddat_glob[id];
    }
}

__global__ void 
keep_pml(int nglob,int nspec, const int *ibool,
        const int * is_cpml,
        const realw* dat_bak,
         realw_p dat_glob)
{
    int bx = blockIdx.x + gridDim.x*blockIdx.y;
    if(bx >= nspec) return;
    int ispec = bx;
    int tx = threadIdx.x;
    if(tx >= NGLL3) tx = NGLL3 - 1; 
    int offset = ispec * NGLL3PAD + tx;
    int iglob = ibool[offset] - 1;

    if(is_cpml[ispec] && threadIdx.x < NGLL3) {
        atomicExch(&dat_glob[iglob],dat_bak[offset]);
    }
}
 
/**
 * @brief assemble mpi buffers across different procs
 * @param buf_sd/rv mpi buffers, shape(num_intfs,max_nibool,ndim), it could be host/device memory
 * @param my_neighbors shape(num_intfs) all neigbors need to exchange information
 * @param num_intfs no. of neighbors
 * @param nibool shape(num_intfs) no. of points needto exchange for each neighbor
 * @param req_sd/rv MPI asynchronize requests, shape(num_intfs)
 */
void assemble_asyn_send(int ndim,realw *buf_sd, realw *buf_rv,
                        int num_intfs,int max_nibool,const int* my_neighbors,
                        const int *nibool,MPI_Request *req_sd,
                        MPI_Request *req_rv,int tag)
{
    for(int it = 0; it < num_intfs; it ++) {
        int iloc = it * max_nibool * ndim;
        int npts = ndim * nibool[it];

        MPI_Isend(buf_sd + iloc,npts,MPI_FLOAT,
                  my_neighbors[it],tag,MPI_COMM_WORLD,
                  &req_sd[it]);
        MPI_Irecv(buf_rv + iloc,npts,MPI_FLOAT,
                  my_neighbors[it],tag,MPI_COMM_WORLD,
                  &req_rv[it]);
    }
}


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
                realw_p dat_glob,realw_p ddat_glob,
                const int * irregular_element_number, const realw * h_xix_regular,
                const int * h_nspec_irregular)
{
    // allocate memory on device
    cudaError_t err;
    #define ALLOCEM(type,a,n) type * d_ ## a;\
                            cudaMalloc((void**)(&d_ ## a),sizeof(type)*n*NGLL3PAD);\
                            err =cudaMemcpy2D(d_ ## a,NGLL3PAD*sizeof(type),a,\
                                        NGLL3*sizeof(realw),NGLL3*sizeof(realw),\
                                    n,cudaMemcpyHostToDevice); \
                            if(err != cudaSuccess) {printf("error here in !\n");}
    #define ALLOC(type,a,n) type * d_ ## a; \
                            cudaMalloc((void**)(&d_ ## a),sizeof(type)*n);  \
                            cudaMemcpy(d_ ## a, a,sizeof(type)*n,cudaMemcpyHostToDevice);
    
    #define GET(type,a) type a = * h_ ## a
    GET(int,nspec); GET(int,nglob);
    GET(int,nstep); GET(int,nspec_outer);
    GET(int,nspec_inner);
    GET(realw,cv); GET(realw,ch);
    GET(int,num_intfs); GET(int,max_nibool);
    GET(int,is_sph);
    GET(int,nspec_irregular);
    GET(realw,xix_regular);
    #undef GET
    
    // allocate a new rotate with shape(nspec,3,NGLL3)
    realw * rotate = new realw[nspec*3*NGLL3];
    for(int ispec = 0; ispec < nspec; ispec ++) {
        for(int i = 0; i < 3; i ++) {
        for(int igll = 0; igll < NGLL3; igll ++) {
            rotate[(ispec*3+i)*NGLL3+igll] = rotate0[(ispec*NGLL3+igll)*3+i];
        }}
    }
    ALLOCEM(realw,rotate,nspec*3);
    ALLOCEM(realw,xix,nspec_irregular);   
    ALLOCEM(realw,xiy,nspec_irregular);   ALLOCEM(realw,jaco,nspec);
    ALLOCEM(realw,xiz,nspec_irregular);   ALLOCEM(realw,dat_bak,nspec);
    ALLOCEM(realw,etax,nspec_irregular);  ALLOCEM(realw,gamx,nspec_irregular); 
    ALLOCEM(realw,etay,nspec_irregular);  ALLOCEM(realw,gamy,nspec_irregular); 
    ALLOCEM(realw,etaz,nspec_irregular);  ALLOCEM(realw,gamz,nspec_irregular);
    ALLOCEM(int,ibool,nspec); 
    ALLOC(realw,wgllwgll_xy,NGLL2); ALLOC(realw,hprimeT,NGLL2);
    ALLOC(realw,wgllwgll_xz,NGLL2); ALLOC(realw,hprime_wgll,NGLL2);
    ALLOC(realw,wgllwgll_yz,NGLL2); ALLOC(realw,ddat_glob,nglob);
    ALLOC(realw,rvol,nglob);
    ALLOC(realw,dat_glob,nglob);
    ALLOC(int,is_CPML,nspec);
    ALLOC(int,irregular_element_number,nspec);

    // phase ispec
    ALLOC(int,phase_ispec_inner_elastic,*nspec_max*2);

    // communication arrays
    realw * buf_sd, *buf_rv;
    cudaMallocHost(&buf_rv,sizeof(realw)*num_intfs*max_nibool);
    cudaMallocHost(&buf_sd,sizeof(realw)*num_intfs*max_nibool);
    ALLOC(realw,buf_sd,num_intfs*max_nibool);
    ALLOC(realw,buf_rv,num_intfs*max_nibool);
    ALLOC(int,nibool_intf,num_intfs);
    ALLOC(int,ibool_intf,num_intfs*max_nibool);
    MPI_Request req_sd[num_intfs],req_rv[num_intfs];

    // gpu resource
    int nb_glob = (nglob + NGLL3PAD - 1) / NGLL3PAD;

    // loop time step
    int nelmnts;
    for(int it = 0; it < nstep; it ++) {
        // set ddat_glob to zero
        cudaMemset(d_ddat_glob,0,sizeof(realw)*nglob);
        for(int iphase = 1; iphase <=2; iphase ++) {
            if(iphase == 1) {
                nelmnts = nspec_outer;
            }
            else {
                nelmnts = nspec_inner;
            }

            // launch kernel;
            kernel_smooth_sph_pde <<<nelmnts,NGLL3PAD>>> (
                is_sph,nelmnts,iphase,*nspec_max,d_phase_ispec_inner_elastic,d_xix,
                d_xiy,d_xiz,d_etax,d_etay,d_etaz,d_gamx,d_gamy,d_gamz,
                d_jaco,d_rotate,d_wgllwgll_xy,d_wgllwgll_xz,d_wgllwgll_yz,
                d_hprimeT,d_hprime_wgll,d_ibool,cv,ch,d_dat_glob,d_ddat_glob,
                d_irregular_element_number,xix_regular
            );

            // communication if required
            if(*nprocs != 1) {
                int nblock1 = (max_nibool * num_intfs + NGLL3PAD - 1)/ NGLL3PAD; 
                if(iphase == 1) {
                    kernel_prepare_boundary_matrix <<<nblock1,NGLL3PAD >>> (
                        d_ddat_glob,1,num_intfs,max_nibool,d_nibool_intf,
                        d_ibool_intf,d_buf_sd,false
                    );
                    size_t size_cp = max_nibool*num_intfs*sizeof(realw);
                    cudaMemcpy(buf_sd,d_buf_sd,size_cp,cudaMemcpyDeviceToHost);
                    
                    assemble_asyn_send(1,buf_sd,buf_rv,num_intfs,max_nibool,
                                    my_neighbors,nibool_intf,req_sd,req_rv,
                                    15141);
                }
                else {
                    MPI_Waitall(num_intfs,req_rv,MPI_STATUS_IGNORE);
                    size_t size_cp = max_nibool*num_intfs*sizeof(realw);
                    cudaMemcpy(d_buf_rv,buf_rv,size_cp,cudaMemcpyHostToDevice);
                    kernel_prepare_boundary_matrix <<<nblock1,NGLL3PAD >>> (
                        d_ddat_glob,1,num_intfs,max_nibool,d_nibool_intf,
                        d_ibool_intf,d_buf_rv,true
                    );
                    MPI_Waitall(num_intfs,req_sd,MPI_STATUS_IGNORE);
                }
            } 
        }
        
        // mass matrix
        apply_mass <<<nb_glob,NGLL3PAD >>> (
            nglob,d_rvol,d_dat_glob,d_ddat_glob
        );

        // keep pml
        keep_pml <<<nspec,NGLL3PAD >>> (
            nglob,nspec,d_ibool,d_is_CPML,d_dat_bak,d_dat_glob
        );
    }

    // copy data to ddat
    cudaMemcpy(dat_glob,d_dat_glob,nglob*sizeof(realw),cudaMemcpyDeviceToHost);

    // copy to dat
    for(int ispec = 0; ispec < nspec; ispec ++) {
        for(int igll = 0; igll < NGLL3; igll ++) {
            int idx = ispec*NGLL3+igll;
            int iglob = ibool[idx] - 1;
            dat[idx] = dat_glob[iglob];
        }
    }

    // free space
    cudaFree(d_rotate);

    cudaFree(d_xix);   
    cudaFree(d_xiy);   cudaFree(d_jaco);
    cudaFree(d_xiz);
    cudaFree(d_etax);  cudaFree(d_gamx); 
    cudaFree(d_etay);  cudaFree(d_gamy); 
    cudaFree(d_etaz);  cudaFree(d_gamz);
    cudaFree(d_ibool); cudaFree(d_wgllwgll_xy);
    cudaFree(d_wgllwgll_xz);cudaFree(d_wgllwgll_yz);
    cudaFree(d_hprimeT); cudaFree(d_hprime_wgll);
    cudaFree(d_ddat_glob);
    cudaFree(d_phase_ispec_inner_elastic);
    cudaFreeHost(buf_rv); cudaFreeHost(buf_sd);
    cudaFree(d_buf_sd); cudaFree(d_buf_rv);
    cudaFree(d_nibool_intf);
    cudaFree(d_ibool_intf);
    cudaFree(d_is_CPML);
    cudaFree(d_irregular_element_number);
    cudaFree(d_dat_glob); cudaFree(d_dat_bak);
    cudaFree(d_rvol);
    delete[] rotate;

    #undef ALLOCEM
    #undef ALLOC
}