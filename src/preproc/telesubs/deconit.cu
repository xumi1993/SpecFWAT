#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/complex.h>
#include <cmath>
#include <cstdio>
#include <vector>

#define PI 3.14159265358979323846

// Error checking macro
#define CHECK_CUDA(call) { \
    cudaError_t err = call; \
    if (err != cudaSuccess) { \
        fprintf(stderr, "CUDA Error: %s at line %d\n", cudaGetErrorString(err), __LINE__); \
        exit(1); \
    } \
}

#define CHECK_CUFFT(call) { \
    cufftResult err = call; \
    if (err != CUFFT_SUCCESS) { \
        fprintf(stderr, "CUFFT Error: %d at line %d\n", err, __LINE__); \
        exit(1); \
    } \
}

// Kernel to initialize Gaussian filter
__global__ void init_gauss_kernel(cufftComplex* gauss, int nft, float dt, float f0) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nft) return;

    float df = 1.0f / (nft * dt);
    int nft21 = nft / 2 + 1;
    
    float val = 0.0f;
    if (i < nft21) {
        float freq = df * i;
        val = expf(-0.25f * powf(2.0f * PI * freq / f0, 2.0f)) / dt;
    } else {
        // Symmetric for negative frequencies (stored in upper half)
        // i goes from nft21 to nft-1. 
        // Corresponding positive index is nft - i?
        // Fortran: gauss1(i) = gauss1(2*nft21-i)
        // 2*nft21 = 2*(nft/2 + 1) = nft + 2.
        // So index is nft + 2 - i. Note Fortran is 1-based.
        // C is 0-based.
        // Fortran i=1 corresponds to freq 0. C i=0.
        // Fortran i=nft21 corresponds to Nyquist. C i=nft/2.
        // Fortran loop i=nft21+1 to nft.
        // C loop i=nft/2+1 to nft-1.
        // Mirror index:
        // Fortran: k = 2*nft21 - i.
        // C: k = (nft + 2 - (i+1)) - 1 = nft - i.
        // Let's verify. i=nft-1 -> k=1. Correct.
        
        int k = nft - i;
        float freq = df * k;
        val = expf(-0.25f * powf(2.0f * PI * freq / f0, 2.0f)) / dt;
    }
    gauss[i] = make_cuComplex(val, 0.0f);
}

// Kernel to apply filter: Z = Z * G * dt
__global__ void apply_filter_kernel(cufftComplex* z, cufftComplex* g, int nft, float dt) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nft) return;
    
    cufftComplex val = z[i];
    cufftComplex filter = g[i];
    // Complex mul: (a+bi)(c+di) = (ac-bd) + (ad+bc)i
    float re = val.x * filter.x - val.y * filter.y;
    float im = val.x * filter.y + val.y * filter.x;
    
    z[i] = make_cuComplex(re * dt, im * dt);
}

// Kernel for correlation: RW = R * conj(W)
__global__ void correl_kernel(cufftComplex* rw, cufftComplex* r, cufftComplex* w, int nft) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nft) return;
    
    cufftComplex rv = r[i];
    cufftComplex wv = w[i];
    // conj(w) = x - iy
    // (rx + iry)(wx - iwy) = (rx*wx + ry*wy) + i(ry*wx - rx*wy)
    
    float re = rv.x * wv.x + rv.y * wv.y;
    float im = rv.y * wv.x - rv.x * wv.y;
    
    rw[i] = make_cuComplex(re, im);
}

// Kernel to update R_freq
// R_freq[k] -= amp * exp(-i*2pi*k*(lag)/N) * G[k] * WF[k] * dt^2
// Note: lag is 0-based index here.
__global__ void update_r_freq_kernel(cufftComplex* r_freq, cufftComplex* g, cufftComplex* wf, 
                                     int nft, float dt, float amp, int lag) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nft) return;
    
    float angle = -2.0f * PI * k * lag / nft;
    float cos_a = cosf(angle);
    float sin_a = sinf(angle);
    
    // Spike response in freq: amp * (cos + i sin)
    // Combined filter: G * WF * dt * dt
    
    cufftComplex G = g[k];
    cufftComplex WF = wf[k];
    
    // GWF = G * WF
    float gwf_re = G.x * WF.x - G.y * WF.y;
    float gwf_im = G.x * WF.y + G.y * WF.x;
    
    // Scale by amp * dt^2
    float scale = amp * dt * dt;
    gwf_re *= scale;
    gwf_im *= scale;
    
    // Multiply by exp(angle)
    // (gwf_re + i gwf_im) * (cos + i sin)
    float update_re = gwf_re * cos_a - gwf_im * sin_a;
    float update_im = gwf_re * sin_a + gwf_im * cos_a;
    
    r_freq[k].x -= update_re;
    r_freq[k].y -= update_im;
}

// Kernel to update p0 at a specific index
__global__ void update_p0_kernel(double* p0, int idx, double amp) {
    p0[idx] += amp;
}

// Kernel for final phase shift and filter
// P = P * G * dt * exp(i * 2pi * k * shift_i / N)
__global__ void final_filter_shift_kernel(cufftComplex* p, cufftComplex* g, int nft, float dt, int shift_i) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k >= nft) return;
    
    float angle = -2.0f * PI * k * shift_i / nft;

    float cos_a = cosf(angle);
    float sin_a = sinf(angle);
    
    cufftComplex val = p[k];
    cufftComplex filter = g[k];
    
    // val * filter * dt
    float re = val.x * filter.x - val.y * filter.y;
    float im = val.x * filter.y + val.y * filter.x;
    re *= dt;
    im *= dt;
    
    // * exp(angle)
    float final_re = re * cos_a - im * sin_a;
    float final_im = re * sin_a + im * cos_a;
    
    p[k] = make_cuComplex(final_re, final_im);
}

// Helper to copy double to complex float
__global__ void double2complex(const double* src, cufftComplex* dst, int n, int nft) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= nft) return;
    if (i < n) {
        dst[i] = make_cuComplex((float)src[i], 0.0f);
    } else {
        dst[i] = make_cuComplex(0.0f, 0.0f);
    }
}

// Helper to copy complex float to double (real part)
__global__ void complex2double(const cufftComplex* src, double* dst, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    dst[i] = (double)src[i].x;
}

// Helper to copy complex float to double (real part) with scaling (1/N)
__global__ void complex2double_scaled(const cufftComplex* src, double* dst, int n, float scale) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    dst[i] = (double)(src[i].x * scale);
}

// Helper to copy double to double
__global__ void double2double(const double* src, double* dst, int n) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    dst[i] = src[i];
}

// Functor for abs
struct abs_functor {
    __host__ __device__
    float operator()(const double& x) const {
        return fabsf((float)x);
    }
};

// Functor for square
struct square_functor {
    __host__ __device__
    double operator()(const double& x) const {
        return x * x;
    }
};

// Functor for complex magnitude square
struct complex_norm_sq_functor {
    __host__ __device__
    double operator()(const cufftComplex& x) const {
        return (double)(x.x * x.x + x.y * x.y);
    }
};

// Functor for abs comparison
struct AbsLess {
    __host__ __device__
    bool operator()(const double& a, const double& b) const {
        return fabs(a) < fabs(b);
    }
};

extern "C" void deconit_cuda(
    double* h_utr, double* h_wtr, int nt, int nft, 
    float dt, float tshift, float f0, 
    int maxit, float minderr, int ipart, 
    double* h_rfi
) {
    // Allocate device memory
    double *d_utr_time, *d_wtr_time, *d_rw_time, *d_p0_time, *d_rfi;
    cufftComplex *d_gauss, *d_wf, *d_U_freq, *d_W_freq, *d_R_freq, *d_RW_freq, *d_P_freq;
    
    CHECK_CUDA(cudaMalloc(&d_utr_time, nft * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_wtr_time, nft * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rw_time, nft * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_p0_time, nft * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rfi, nt * sizeof(double)));
    
    CHECK_CUDA(cudaMalloc(&d_gauss, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_wf, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_U_freq, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_W_freq, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_R_freq, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_RW_freq, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_P_freq, nft * sizeof(cufftComplex)));
    
    // Initialize p0 to 0
    CHECK_CUDA(cudaMemset(d_p0_time, 0, nft * sizeof(double)));
    
    // Copy input
    CHECK_CUDA(cudaMemcpy(d_utr_time, h_utr, nt * sizeof(double), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_wtr_time, h_wtr, nt * sizeof(double), cudaMemcpyHostToDevice));
    // Zero pad
    if (nft > nt) {
        CHECK_CUDA(cudaMemset(d_utr_time + nt, 0, (nft - nt) * sizeof(double)));
        CHECK_CUDA(cudaMemset(d_wtr_time + nt, 0, (nft - nt) * sizeof(double)));
    }
    
    // Create plans
    cufftHandle plan;
    CHECK_CUFFT(cufftPlan1d(&plan, nft, CUFFT_C2C, 1));
    
    int blockSize = 256;
    int numBlocks = (nft + blockSize - 1) / blockSize;
    
    // 1. Compute Gauss
    init_gauss_kernel<<<numBlocks, blockSize>>>(d_gauss, nft, dt, f0);
    
    // 2. Compute WF = FFT(wtr)
    // Need to convert double wtr to complex
    double2complex<<<numBlocks, blockSize>>>(d_wtr_time, d_W_freq, nt, nft);
    CHECK_CUFFT(cufftExecC2C(plan, d_W_freq, d_wf, CUFFT_FORWARD)); // d_wf holds FFT(wtr)
    
    // 3. Compute U_freq = FFT(utr) * G * dt
    double2complex<<<numBlocks, blockSize>>>(d_utr_time, d_U_freq, nt, nft);
    CHECK_CUFFT(cufftExecC2C(plan, d_U_freq, d_U_freq, CUFFT_FORWARD));
    apply_filter_kernel<<<numBlocks, blockSize>>>(d_U_freq, d_gauss, nft, dt);
    
    // 4. Compute W_freq = FFT(wtr) * G * dt
    // We already have FFT(wtr) in d_wf. Copy it to d_W_freq then filter.
    CHECK_CUDA(cudaMemcpy(d_W_freq, d_wf, nft * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    apply_filter_kernel<<<numBlocks, blockSize>>>(d_W_freq, d_gauss, nft, dt);
    
    // 5. Compute powerU and powerW
    // Need IFFT of U_freq and W_freq to get filtered time domain signals
    // Use d_RW_freq as temp buffer for IFFT
    
    // U filtered
    CHECK_CUDA(cudaMemcpy(d_RW_freq, d_U_freq, nft * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    CHECK_CUFFT(cufftExecC2C(plan, d_RW_freq, d_RW_freq, CUFFT_INVERSE));
    // Store in d_utr_time (reuse buffer, we don't need raw utr anymore)
    complex2double_scaled<<<numBlocks, blockSize>>>(d_RW_freq, d_utr_time, nft, 1.0f/nft);
    
    // W filtered
    CHECK_CUDA(cudaMemcpy(d_RW_freq, d_W_freq, nft * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    CHECK_CUFFT(cufftExecC2C(plan, d_RW_freq, d_RW_freq, CUFFT_INVERSE));
    // Store in d_wtr_time
    complex2double_scaled<<<numBlocks, blockSize>>>(d_RW_freq, d_wtr_time, nft, 1.0f/nft);
    
    // Compute powers
    thrust::device_ptr<double> t_utr(d_utr_time);
    thrust::device_ptr<double> t_wtr(d_wtr_time);
    
    double powerU = thrust::transform_reduce(t_utr, t_utr + nft, square_functor(), 0.0, thrust::plus<double>());
    double powerW = thrust::transform_reduce(t_wtr, t_wtr + nft, square_functor(), 0.0, thrust::plus<double>());
    
    // 6. Init R_freq = U_freq
    CHECK_CUDA(cudaMemcpy(d_R_freq, d_U_freq, nft * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    
    double sumsq_i = 1.0;
    double d_error = 100.0 * powerU + minderr;
    int maxlag = (ipart == 0) ? (nft / 2) : nft;
    
    int i = 0;
    while (fabs(d_error) > minderr && i < maxit) {
        i++;
        
        // RW = R * conj(W)
        correl_kernel<<<numBlocks, blockSize>>>(d_RW_freq, d_R_freq, d_W_freq, nft);
        
        // IFFT RW
        CHECK_CUFFT(cufftExecC2C(plan, d_RW_freq, d_RW_freq, CUFFT_INVERSE));
        
        // Convert to double and scale by 1/N
        complex2double_scaled<<<numBlocks, blockSize>>>(d_RW_freq, d_rw_time, nft, 1.0f/nft);
        
        // Find max in rw[0..maxlag]
        thrust::device_ptr<double> t_rw(d_rw_time);
        thrust::device_ptr<double> max_ptr = thrust::max_element(t_rw, t_rw + maxlag, AbsLess());
        
        int idx = max_ptr - t_rw;
        double max_val = *max_ptr; // This is rw(idx)
        
        // amp = rw(idx) / powerW / dt
        // Note: rw in code was normalized by sum(wflt**2) which is powerW.
        double amp = (max_val / powerW) / dt;
        
        // Update p0
        update_p0_kernel<<<1, 1>>>(d_p0_time, idx, amp);

        // Update R_freq
        update_r_freq_kernel<<<numBlocks, blockSize>>>(d_R_freq, d_gauss, d_wf, nft, dt, (float)amp, idx);
        
        // Calculate RMS
        // sumsq = sum(|R_freq|^2) / N / powerU
        // Note: Parseval: sum(r^2) = sum(|R|^2)/N
        thrust::device_ptr<cufftComplex> t_R_freq(d_R_freq);
        double sum_R_sq = thrust::transform_reduce(t_R_freq, t_R_freq + nft, complex_norm_sq_functor(), 0.0, thrust::plus<double>());
        
        double sumsq = (sum_R_sq / nft) / powerU;
        
        d_error = 100.0 * (sumsq_i - sumsq);
        sumsq_i = sumsq;
    }
    
    // 7. Final result
    // FFT p0
    double2complex<<<numBlocks, blockSize>>>(d_p0_time, d_P_freq, nft, nft);
    CHECK_CUFFT(cufftExecC2C(plan, d_P_freq, d_P_freq, CUFFT_FORWARD));
    
    // Filter and phase shift
    int shift_i = (int)(tshift / dt + 0.5f);
    final_filter_shift_kernel<<<numBlocks, blockSize>>>(d_P_freq, d_gauss, nft, dt, shift_i);
    
    // IFFT
    CHECK_CUFFT(cufftExecC2C(plan, d_P_freq, d_P_freq, CUFFT_INVERSE));
    
    // Copy to output
    complex2double_scaled<<<numBlocks, blockSize>>>(d_P_freq, d_rfi, nt, 1.0f/nft); // Only copy nt
    CHECK_CUDA(cudaMemcpy(h_rfi, d_rfi, nt * sizeof(double), cudaMemcpyDeviceToHost));
    
    // Cleanup
    cufftDestroy(plan);
    cudaFree(d_utr_time);
    cudaFree(d_wtr_time);
    cudaFree(d_rw_time);
    cudaFree(d_p0_time);
    cudaFree(d_gauss);
    cudaFree(d_wf);
    cudaFree(d_U_freq);
    cudaFree(d_W_freq);
    cudaFree(d_R_freq);
    cudaFree(d_RW_freq);
    cudaFree(d_P_freq);
    cudaFree(d_rfi);
}


// ---------------------------------------------------------------------------
// Batched iterative deconvolution: nb traces of nt samples (trace b at
// offset b*nt) are deconvolved together. Each trace follows the same steps
// as deconit_cuda and stops on its own (per-trace convergence flags); the
// host only checks the flags every BATCH_CHECK_EVERY iterations.
// ---------------------------------------------------------------------------

#define BATCH_BLOCK 256
#define BATCH_CHECK_EVERY 20

// real input (nb x nt, double) -> zero-padded complex (nb x nft)
__global__ void batch_pad_kernel(const double* src, cufftComplex* dst, int nt, int nft, int nb) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    int b = i / nft, t = i % nft;
    dst[i] = make_cuComplex(t < nt ? (float)src[(long)b * nt + t] : 0.0f, 0.0f);
}

// z = z * g * dt (g shared by all traces)
__global__ void batch_filter_kernel(cufftComplex* z, const cufftComplex* g, int nft, int nb, float dt) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    cufftComplex v = z[i], f = g[i % nft];
    z[i] = make_cuComplex((v.x * f.x - v.y * f.y) * dt, (v.x * f.y + v.y * f.x) * dt);
}

// out = a * b (complex), b shared by all traces
__global__ void batch_mul_shared_kernel(cufftComplex* out, const cufftComplex* a, const cufftComplex* b,
                                        int nft, int nb) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    cufftComplex x = a[i], y = b[i % nft];
    out[i] = make_cuComplex(x.x * y.x - x.y * y.y, x.x * y.y + x.y * y.x);
}

// rw = r * conj(w)
__global__ void batch_correl_kernel(cufftComplex* rw, const cufftComplex* r, const cufftComplex* w,
                                    int nft, int nb) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    cufftComplex rv = r[i], wv = w[i];
    rw[i] = make_cuComplex(rv.x * wv.x + rv.y * wv.y, rv.y * wv.x - rv.x * wv.y);
}

// per-trace sum of (scale * Re x)^2 over nft samples; one block per trace
__global__ void batch_power_kernel(const cufftComplex* x, int nft, float scale, double* out) {
    __shared__ double s[BATCH_BLOCK];
    const cufftComplex* xb = x + (long)blockIdx.x * nft;
    double acc = 0.0;
    for (int k = threadIdx.x; k < nft; k += blockDim.x) {
        double v = (double)(xb[k].x * scale);
        acc += v * v;
    }
    s[threadIdx.x] = acc;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if (threadIdx.x < h) s[threadIdx.x] += s[threadIdx.x + h];
        __syncthreads();
    }
    if (threadIdx.x == 0) out[blockIdx.x] = s[0];
}

// per trace: lag of max |rw| in [0, maxlag) (first one on ties), spike update of p0
__global__ void batch_pick_kernel(const cufftComplex* rw, int nft, int maxlag, float scale,
                                  const int* active, const double* powerW, float dt,
                                  double* p0, double* amp, int* lag) {
    __shared__ double sv[BATCH_BLOCK];
    __shared__ int si[BATCH_BLOCK];
    int b = blockIdx.x;
    if (!active[b]) {
        if (threadIdx.x == 0) amp[b] = 0.0;
        return;
    }
    const cufftComplex* rb = rw + (long)b * nft;
    double best = -1.0; int ibest = 0;
    for (int k = threadIdx.x; k < maxlag; k += blockDim.x) {
        double v = fabs((double)(rb[k].x * scale));
        if (v > best) { best = v; ibest = k; }
    }
    sv[threadIdx.x] = best; si[threadIdx.x] = ibest;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if (threadIdx.x < h) {
            double v2 = sv[threadIdx.x + h]; int i2 = si[threadIdx.x + h];
            if (v2 > sv[threadIdx.x] || (v2 == sv[threadIdx.x] && i2 < si[threadIdx.x])) {
                sv[threadIdx.x] = v2; si[threadIdx.x] = i2;
            }
        }
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        int k = si[0];
        double val = (double)(rb[k].x * scale);
        double a = (val / powerW[b]) / dt;
        p0[(long)b * nft + k] += a;
        amp[b] = a;
        lag[b] = k;
    }
}

// R[b,k] -= amp_b * exp(-i 2 pi k lag_b / nft) * G * WF * dt^2
__global__ void batch_update_r_kernel(cufftComplex* r, const cufftComplex* gwf, int nft, int nb, float dt,
                                      const int* active, const double* amp, const int* lag) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    int b = i / nft, k = i % nft;
    if (!active[b]) return;
    float angle = -2.0f * PI * k * lag[b] / nft;
    float cos_a = cosf(angle), sin_a = sinf(angle);
    float scale = (float)amp[b] * dt * dt;
    float gre = gwf[i].x * scale, gim = gwf[i].y * scale;
    r[i].x -= gre * cos_a - gim * sin_a;
    r[i].y -= gre * sin_a + gim * cos_a;
}

// per trace: relative residual power and convergence test; one block per trace
__global__ void batch_converge_kernel(const cufftComplex* r, int nft, const double* powerU, float minderr,
                                      double* sumsq_i, int* active) {
    __shared__ double s[BATCH_BLOCK];
    int b = blockIdx.x;
    if (!active[b]) return;
    const cufftComplex* rb = r + (long)b * nft;
    double acc = 0.0;
    for (int k = threadIdx.x; k < nft; k += blockDim.x)
        acc += (double)(rb[k].x * rb[k].x + rb[k].y * rb[k].y);
    s[threadIdx.x] = acc;
    __syncthreads();
    for (int h = blockDim.x / 2; h > 0; h >>= 1) {
        if (threadIdx.x < h) s[threadIdx.x] += s[threadIdx.x + h];
        __syncthreads();
    }
    if (threadIdx.x == 0) {
        double sumsq = (s[0] / nft) / powerU[b];
        double derr = 100.0 * (sumsq_i[b] - sumsq);
        sumsq_i[b] = sumsq;
        if (!(fabs(derr) > minderr)) active[b] = 0;
    }
}

// p0 (double) -> complex
__global__ void batch_p0_to_complex_kernel(const double* p0, cufftComplex* z, long n) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    z[i] = make_cuComplex((float)p0[i], 0.0f);
}

// P = P * G * dt * exp(-i 2 pi k shift_i / nft)
__global__ void batch_final_filter_shift_kernel(cufftComplex* p, const cufftComplex* g, int nft, int nb,
                                                float dt, int shift_i) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nft) return;
    int k = i % nft;
    float angle = -2.0f * PI * k * shift_i / nft;
    float cos_a = cosf(angle), sin_a = sinf(angle);
    cufftComplex v = p[i], f = g[k];
    float re = (v.x * f.x - v.y * f.y) * dt;
    float im = (v.x * f.y + v.y * f.x) * dt;
    p[i] = make_cuComplex(re * cos_a - im * sin_a, re * sin_a + im * cos_a);
}

// first nt samples of Re(x) / nft -> output (nb x nt, double)
__global__ void batch_unpad_kernel(const cufftComplex* x, double* dst, int nt, int nft, int nb, float scale) {
    long i = (long)blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= (long)nb * nt) return;
    int b = i / nt, t = i % nt;
    dst[i] = (double)(x[(long)b * nft + t].x * scale);
}

extern "C" void deconit_cuda_batch(
    double* h_utr, double* h_wtr, int nt, int nft, int nb,
    float dt, float tshift, float f0,
    int maxit, float minderr, int ipart,
    double* h_rfi
) {
    if (nb <= 0) return;
    long nall = (long)nb * nft;
    int grid_all = (int)((nall + BATCH_BLOCK - 1) / BATCH_BLOCK);
    int grid_out = (int)(((long)nb * nt + BATCH_BLOCK - 1) / BATCH_BLOCK);
    float inv_nft = 1.0f / nft;

    double *d_in, *d_p0, *d_rfi, *d_powerU, *d_powerW, *d_sumsq_i, *d_amp;
    int *d_active, *d_lag;
    cufftComplex *d_gauss, *d_wf, *d_W, *d_R, *d_gwf, *d_tmp;
    CHECK_CUDA(cudaMalloc(&d_in, (long)nb * nt * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_p0, nall * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_rfi, (long)nb * nt * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_powerU, nb * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_powerW, nb * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_sumsq_i, nb * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_amp, nb * sizeof(double)));
    CHECK_CUDA(cudaMalloc(&d_active, nb * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_lag, nb * sizeof(int)));
    CHECK_CUDA(cudaMalloc(&d_gauss, nft * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_wf, nall * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_W, nall * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_R, nall * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_gwf, nall * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMalloc(&d_tmp, nall * sizeof(cufftComplex)));
    CHECK_CUDA(cudaMemset(d_p0, 0, nall * sizeof(double)));

    cufftHandle plan;
    CHECK_CUFFT(cufftPlan1d(&plan, nft, CUFFT_C2C, nb));
    init_gauss_kernel<<<(nft + BATCH_BLOCK - 1) / BATCH_BLOCK, BATCH_BLOCK>>>(d_gauss, nft, dt, f0);

    // WF = FFT(w); W = WF * G * dt; GWF = G * WF
    CHECK_CUDA(cudaMemcpy(d_in, h_wtr, (long)nb * nt * sizeof(double), cudaMemcpyHostToDevice));
    batch_pad_kernel<<<grid_all, BATCH_BLOCK>>>(d_in, d_tmp, nt, nft, nb);
    CHECK_CUFFT(cufftExecC2C(plan, d_tmp, d_wf, CUFFT_FORWARD));
    CHECK_CUDA(cudaMemcpy(d_W, d_wf, nall * sizeof(cufftComplex), cudaMemcpyDeviceToDevice));
    batch_filter_kernel<<<grid_all, BATCH_BLOCK>>>(d_W, d_gauss, nft, nb, dt);
    batch_mul_shared_kernel<<<grid_all, BATCH_BLOCK>>>(d_gwf, d_wf, d_gauss, nft, nb);

    // R = U = FFT(u) * G * dt
    CHECK_CUDA(cudaMemcpy(d_in, h_utr, (long)nb * nt * sizeof(double), cudaMemcpyHostToDevice));
    batch_pad_kernel<<<grid_all, BATCH_BLOCK>>>(d_in, d_tmp, nt, nft, nb);
    CHECK_CUFFT(cufftExecC2C(plan, d_tmp, d_R, CUFFT_FORWARD));
    batch_filter_kernel<<<grid_all, BATCH_BLOCK>>>(d_R, d_gauss, nft, nb, dt);

    // powers of the filtered traces
    CHECK_CUFFT(cufftExecC2C(plan, d_R, d_tmp, CUFFT_INVERSE));
    batch_power_kernel<<<nb, BATCH_BLOCK>>>(d_tmp, nft, inv_nft, d_powerU);
    CHECK_CUFFT(cufftExecC2C(plan, d_W, d_tmp, CUFFT_INVERSE));
    batch_power_kernel<<<nb, BATCH_BLOCK>>>(d_tmp, nft, inv_nft, d_powerW);

    // initial state: active unless a trace has no energy
    std::vector<double> powU(nb), powW(nb), ones(nb, 1.0);
    std::vector<int> act(nb);
    CHECK_CUDA(cudaMemcpy(powU.data(), d_powerU, nb * sizeof(double), cudaMemcpyDeviceToHost));
    CHECK_CUDA(cudaMemcpy(powW.data(), d_powerW, nb * sizeof(double), cudaMemcpyDeviceToHost));
    for (int b = 0; b < nb; b++)
        act[b] = (powU[b] > 0.0 && powW[b] > 0.0 && fabs(100.0 * powU[b] + minderr) > minderr) ? 1 : 0;
    CHECK_CUDA(cudaMemcpy(d_active, act.data(), nb * sizeof(int), cudaMemcpyHostToDevice));
    CHECK_CUDA(cudaMemcpy(d_sumsq_i, ones.data(), nb * sizeof(double), cudaMemcpyHostToDevice));

    int maxlag = (ipart == 0) ? (nft / 2) : nft;
    for (int it = 1; it <= maxit; it++) {
        batch_correl_kernel<<<grid_all, BATCH_BLOCK>>>(d_tmp, d_R, d_W, nft, nb);
        CHECK_CUFFT(cufftExecC2C(plan, d_tmp, d_tmp, CUFFT_INVERSE));
        batch_pick_kernel<<<nb, BATCH_BLOCK>>>(d_tmp, nft, maxlag, inv_nft, d_active, d_powerW, dt,
                                              d_p0, d_amp, d_lag);
        batch_update_r_kernel<<<grid_all, BATCH_BLOCK>>>(d_R, d_gwf, nft, nb, dt, d_active, d_amp, d_lag);
        batch_converge_kernel<<<nb, BATCH_BLOCK>>>(d_R, nft, d_powerU, minderr, d_sumsq_i, d_active);
        if (it % BATCH_CHECK_EVERY == 0) {
            CHECK_CUDA(cudaMemcpy(act.data(), d_active, nb * sizeof(int), cudaMemcpyDeviceToHost));
            int nact = 0;
            for (int b = 0; b < nb; b++) nact += act[b];
            if (nact == 0) break;
        }
    }

    // RF = IFFT(FFT(p0) * G * dt * shift)
    batch_p0_to_complex_kernel<<<grid_all, BATCH_BLOCK>>>(d_p0, d_tmp, nall);
    CHECK_CUFFT(cufftExecC2C(plan, d_tmp, d_tmp, CUFFT_FORWARD));
    int shift_i = (int)(tshift / dt + 0.5f);
    batch_final_filter_shift_kernel<<<grid_all, BATCH_BLOCK>>>(d_tmp, d_gauss, nft, nb, dt, shift_i);
    CHECK_CUFFT(cufftExecC2C(plan, d_tmp, d_tmp, CUFFT_INVERSE));
    batch_unpad_kernel<<<grid_out, BATCH_BLOCK>>>(d_tmp, d_rfi, nt, nft, nb, inv_nft);
    CHECK_CUDA(cudaMemcpy(h_rfi, d_rfi, (long)nb * nt * sizeof(double), cudaMemcpyDeviceToHost));

    cufftDestroy(plan);
    cudaFree(d_in); cudaFree(d_p0); cudaFree(d_rfi);
    cudaFree(d_powerU); cudaFree(d_powerW); cudaFree(d_sumsq_i); cudaFree(d_amp);
    cudaFree(d_active); cudaFree(d_lag);
    cudaFree(d_gauss); cudaFree(d_wf); cudaFree(d_W); cudaFree(d_R); cudaFree(d_gwf); cudaFree(d_tmp);
}
