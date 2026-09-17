// GPL-3.0-or-later. Layer matrices follow SPECFEM's couple_with_injection.f90.
// Double-complex propagation and cuFFT; storage/filter rounding follows CUSTOM_REAL.
#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/complex.h>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>
#ifdef FK_DOUBLE_PRECISION
using Real = double;
#else
using Real = float;
#endif
using Z = thrust::complex<double>;
static_assert(sizeof(Z)==sizeof(cufftDoubleComplex), "complex ABI mismatch");
constexpr double pi = 3.1415926535897932384626433832795;

static void check(cudaError_t e) {
  if (e != cudaSuccess) throw std::runtime_error(cudaGetErrorString(e));
}
static void fftcheck(cufftResult e) {
  if (e != CUFFT_SUCCESS) throw std::runtime_error("cuFFT error " + std::to_string(e));
}
// Restore the solver's CUDA device after FK frees its temporary buffers and FFT plan.
struct DeviceScope {
  int previous;
  explicit DeviceScope(int selected) {
    check(cudaGetDevice(&previous));
    check(cudaSetDevice(selected));
  }
  ~DeviceScope() { cudaSetDevice(previous); }
};
template<class T> struct Buffer {
  T* p=nullptr;
  explicit Buffer(size_t n) { check(cudaMalloc(reinterpret_cast<void**>(&p),n*sizeof(T))); }
  ~Buffer() { if (p) cudaFree(p); }
  Buffer(const Buffer&)=delete;
  Buffer& operator=(const Buffer&)=delete;
  void upload(const T* src, size_t n) { check(cudaMemcpy(p,src,n*sizeof(T),cudaMemcpyHostToDevice)); }
};
struct Plan {
  cufftHandle h=0;
  ~Plan() { if(h) cufftDestroy(h); }
};

__device__ void psv(double om, Z ea, Z eb, double rho, double vs, double height, double p, Z g, Z* m) {
  Z ca=(thrust::exp(om*ea*height)+thrust::exp(-om*ea*height))/2.0;
  Z sa=(thrust::exp(om*ea*height)-thrust::exp(-om*ea*height))/2.0;
  Z cb=(thrust::exp(om*eb*height)+thrust::exp(-om*eb*height))/2.0;
  Z sb=(thrust::exp(om*eb*height)-thrust::exp(-om*eb*height))/2.0;
  Z xa=ea*sa/p, ya=p*sa/ea, xb=eb*sb/p, yb=p*sb/eb;
  double tm=2.0*rho*vs*vs;
  // Row-major local matrix.
  m[0]=ca-g*cb; m[1]=xb-g*ya; m[2]=(ya-xb)/tm; m[3]=(cb-ca)/tm;
  m[4]=xa-g*yb; m[5]=cb-g*ca; m[6]=(ca-cb)/tm; m[7]=(yb-xa)/tm;
  m[8]=tm*(xa-g*g*yb); m[9]=tm*g*(cb-ca); m[10]=m[0]; m[11]=g*yb-xa;
  m[12]=tm*g*(ca-cb); m[13]=tm*(xb-g*g*ya); m[14]=g*ya-xb; m[15]=m[5];
}

__global__ void spectrum(int nl,int np,int nf,int first,int count,const double* layers,
    const Z* ea,const Z* eb,const Z* gamma,const Z* E,const Z* bottom,const Z* cache,
    const double* geom,const int* layer,const int* acoustic,const double* o,Z* out) {
  size_t idx=blockIdx.x*size_t(blockDim.x)+threadIdx.x;
  int npos=nf/2+1;
  if(idx>=size_t(count)*npos) return;
  int q=idx/npos, f=idx%npos, point=first+q;
  // Match the reference's CUSTOM_REAL frequency/omega rounding.
  Real freq=Real(f)*Real(o[6]);
  double om=Real(2.0*pi*freq);
  double p=o[0], phi=o[1];
  const double* x=geom+10*point;
  int l=layer[point];
  Z v[4];
  if(!acoustic[point] && x[2]<=0) {
    Z phase[4]={thrust::exp(om*eb[nl-1]*x[2]),thrust::exp(-om*eb[nl-1]*x[2]),
                thrust::exp(om*ea[nl-1]*x[2]),thrust::exp(-om*ea[nl-1]*x[2])};
    for(int i=0;i<4;i++) {
      v[i]=0;
      for(int j=0;j<4;j++) v[i]+=E[i+4*j]*phase[j]*bottom[j+4*f];
    }
  } else if(!acoustic[point]) {
    Z mat[16];
    psv(om,ea[l],eb[l],layers[3*l],layers[3*l+1],x[9],p,gamma[l],mat);
    double g0=Real(2.0)*Real(layers[3*l+1])*Real(layers[3*l+1])*Real(p)*Real(p);
    for(int i=0;i<4;i++) {
      v[i]=0;
      for(int j=0;j<4;j++) v[i]+=mat[4*i+j]*cache[j+4*(l+nl*f)];
      v[i]*=g0;
    }
  } else {
    Z ca=(thrust::exp(om*ea[l]*x[9])+thrust::exp(-om*ea[l]*x[9]))/2.0;
    Z sa=(thrust::exp(om*ea[l]*x[9])-thrust::exp(-om*ea[l]*x[9]))/2.0;
    const Z* base=cache+4*(min(l,int(o[10])-1)+nl*f);
    v[0]=ca*base[0]-sa*ea[l]*p/layers[3*l]*base[1];
    v[1]=sa*layers[3*l]/(p*ea[l])*base[0]+ca*base[1];
    if(l==int(o[10])) { v[0]=base[0]; v[1]=base[1]; }
    // At a fluid/solid interface the reference chooses the adjacent solid layer
    // for density and applies no partial acoustic propagator.
    // Match the reference's free-surface test.
    if(fabs(x[9]-layers[2]) < 1.e-6*layers[2]) v[1]=0;
  }
  Real delay=Real(Real(p)*Real(x[0]-o[2])*cos(Real(phi)) +
                  Real(p)*Real(x[1]-o[3])*sin(Real(phi)) + o[7]*(-o[4]));
  Real window=exp(-pow(Real(om)*Real(o[5])/Real(2),Real(2)));
  Z stf=double(window)*thrust::exp(Z(0,-om*delay));
  Z fields[5];
  if(!acoustic[point]) {
    fields[0]=stf*v[0]*om;
    fields[1]=stf*v[1]*Z(0,om);
    fields[2]=stf*om*p*(x[3]*v[3]-4.0*x[4]*v[0]);
    fields[3]=stf*om*p*v[2]*Z(0,-1);
    fields[4]=stf*om*p*v[3];
  } else {
    fields[0]=stf*Z(0,-1)*p*p/layers[3*l]*v[1];
    fields[1]=stf*v[0];
    fields[2]=fields[3]=fields[4]=Z(0,1)*v[1]*stf*p;
  }
  for(int j=0;j<5;j++) {
    size_t base=(size_t(q)*5+j)*nf;
    // FFTinv/rspec force DC real and Nyquist zero.
    out[base+f]=f==nf/2 ? Z(0,0) : (f==0 ? Z(fields[j].real(),0) : fields[j]);
    if(f>0 && f<nf/2) out[base+nf-f]=thrust::conj(fields[j]);
  }
}

__device__ Real field_at(const Z* s,int nf,int q,int j,int t,int shift,double scale) {
  int index=(t-shift+nf)%nf;
  return Real(s[(size_t(q)*5+j)*nf+index].real()*scale);
}

__global__ void spline(int nf,int first,int count,const Z* spectral,const double* geom,
    const int* acoustic,const double* o,double* work,double* result) {
  int trace=blockIdx.x*blockDim.x+threadIdx.x;
  if(trace>=count*6) return;
  int q=trace/6, component=trace%6, point=first+q, ns=nf/2;
  const double* x=geom+point*10;
  Real cp=cos(Real(o[1])), sp=sin(Real(o[1]));
  // Adjacent threads filter adjacent traces; keep each time slice contiguous.
  double* c=work+trace;
  size_t stride=size_t(count)*6;
  constexpr double pole=-0.267949192431122706472553658494;
  const double factor=(1.0-pole)*(1.0-1.0/pole);
  for(int t=0;t<nf;t++) {
    Real uz=field_at(spectral,nf,q,1,t,int(o[9]),o[11]);
    Real value=uz;
    if(component<2) {
      Real ur=field_at(spectral,nf,q,0,t,int(o[9]),o[11]);
      value=ur*(component==0 ? cp : sp);
    } else if(component>=3) {
      if(acoustic[point]) {
        value=field_at(spectral,nf,q,2,t,int(o[9]),o[11]);
      } else if(t<ns) {
        Real rr=field_at(spectral,nf,q,2,t,int(o[9]),o[11]);
        Real rz=field_at(spectral,nf,q,3,t,int(o[9]),o[11]);
        Real zz=field_at(spectral,nf,q,4,t,int(o[9]),o[11]);
        Real tt=Real(x[5])*(rr+zz);
        Real xx=rr*cp*cp+tt*sp*sp, xy=cp*sp*(rr-tt), xz=rz*cp;
        Real yy=rr*sp*sp+tt*cp*cp, yz=rz*sp;
        Real nx=Real(x[6]), ny=Real(x[7]), nz=Real(x[8]);
        if(component==3) value=xx*nx+xy*ny+xz*nz;
        if(component==4) value=xy*nx+yy*ny+yz*nz;
        if(component==5) value=xz*nx+yz*ny+zz*nz;
      }
      // For elastic tractions the legacy tmp_t1 tail retains vertical velocity.
      // Preserve it for baseline compatibility instead of silently changing it.
    }
    c[size_t(t)*stride]=double(value)*factor;
  }
  double sum=c[0], power=pole;
  int ninit=min(nf,42);
  for(int t=0;t<ninit;t++) { sum+=power*c[size_t(t)*stride]; power*=pole; }
  c[0]=sum;
  for(int t=1;t<nf;t++) c[size_t(t)*stride]+=pole*c[size_t(t-1)*stride];
  c[size_t(nf-1)*stride]*=pole/(pole-1.0);
  for(int t=nf-2;t>=0;t--) c[size_t(t)*stride]=pole*(c[size_t(t+1)*stride]-c[size_t(t)*stride]);
  for(int t=0;t<ns;t++) {
    Real value=Real(c[size_t(t)*stride]);
    if(ns>40 && t<20) {
      Real taper=Real((1.0-cos(pi*t/20.0))*0.5);
      value=taper*value;
    }
    result[(size_t(t)*count+q)*6+component]=value;
  }
}

extern "C" int fk_cuda_compute(int nl,int np,int nf,int local_rank,
    const double* layers,const Z* ea,const Z* eb,const Z* gamma,const Z* E,const Z* bottom,
    const Z* cache,const double* geom,const int* layer,const int* acoustic,const double* opts,
    double* velocity,double* traction) {
  try {
    int devices=0;
    check(cudaGetDeviceCount(&devices));
    if(devices==0) throw std::runtime_error("No visible CUDA device");
    // Use FWAT's node-local MPI rank within the scheduler-visible devices.
    int selected=local_rank%devices;
    DeviceScope device_scope(selected);
    int npos=nf/2+1, ns=nf/2;
    Buffer<double> dl(3*size_t(nl)), dg(10*size_t(np)), opt(12);
    Buffer<Z> da(nl),db(nl),dgam(nl),de(16),dbo(4*size_t(npos)),dc(4*size_t(nl)*npos);
    Buffer<int> dlayer(np),daco(np);
    dl.upload(layers,3*size_t(nl)); dg.upload(geom,10*size_t(np)); opt.upload(opts,12);
    da.upload(ea,nl); db.upload(eb,nl); dgam.upload(gamma,nl); de.upload(E,16);
    dbo.upload(bottom,4*size_t(npos)); dc.upload(cache,4*size_t(nl)*npos);
    dlayer.upload(layer,np); daco.upload(acoustic,np);
    size_t free_bytes,total_bytes;
    check(cudaMemGetInfo(&free_bytes,&total_bytes));
    // Cap the internal batch at 1024 points, with headroom for cuFFT workspace.
    // Spectra + spline scratch + output determine the per-point memory budget.
    size_t per_point=size_t(nf)*(5*sizeof(Z)+6*sizeof(double)+3*sizeof(double));
    int batch=int(std::min({size_t(1024),size_t(np),free_bytes/2/per_point}));
    if(batch<1) throw std::runtime_error("Insufficient GPU memory for one FK point");
    if(batch>std::numeric_limits<int>::max()/6) throw std::runtime_error("FK batch too large");
    Buffer<Z> spec(size_t(batch)*5*nf);
    Buffer<double> work(size_t(batch)*6*nf), result(size_t(batch)*6*ns);
    std::vector<double> host(size_t(batch)*6*ns);
    Plan plan;
    fftcheck(cufftPlan1d(&plan.h,nf,CUFFT_Z2Z,batch*5));
    std::fprintf(stdout,"FK CUDA: device=%d, batch=%d, points=%d, FFT=%d\n",selected,batch,np,nf);
    for(int first=0;first<np;first+=batch) {
      int count=std::min(batch,np-first);
      // Clear padding for the last partial batch; reuse the same FFT plan.
      check(cudaMemset(spec.p,0,size_t(batch)*5*nf*sizeof(Z)));
      size_t tasks=size_t(count)*npos;
      spectrum<<<(tasks+127)/128,128>>>(nl,np,nf,first,count,dl.p,da.p,db.p,dgam.p,de.p,dbo.p,dc.p,
        dg.p,dlayer.p,daco.p,opt.p,spec.p);
      check(cudaGetLastError());
      fftcheck(cufftExecZ2Z(plan.h,reinterpret_cast<cufftDoubleComplex*>(spec.p),
        reinterpret_cast<cufftDoubleComplex*>(spec.p),CUFFT_INVERSE));
      spline<<<(count*6+127)/128,128>>>(nf,first,count,spec.p,dg.p,daco.p,opt.p,work.p,result.p);
      check(cudaGetLastError());
      check(cudaMemcpy(host.data(),result.p,size_t(count)*6*ns*sizeof(double),cudaMemcpyDeviceToHost));
      for(int t=0;t<ns;t++) for(int q=0;q<count;q++) for(int j=0;j<3;j++) {
        size_t dest=(size_t(t)*np+first+q)*3+j;
        size_t src=(size_t(t)*count+q)*6+j;
        velocity[dest]=host[src]; traction[dest]=host[src+3];
      }
    }
    check(cudaDeviceSynchronize());
    return 0;
  } catch(const std::exception& e) {
    std::fprintf(stderr,"FK CUDA failed: %s\n",e.what());
    return 1;
  }
}
