#include <cuda_runtime.h>
#include <cufft.h>
#include <thrust/device_vector.h>
#include <thrust/extrema.h>
#include <thrust/transform_reduce.h>
#include <thrust/functional.h>
#include <thrust/complex.h>
#include <cmath>
#include <cstdio>

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

