extern "C" void deconit_cuda(
    double* h_utr, double* h_wtr, int nt, int nft, 
    float dt, float tshift, float f0, 
    int maxit, float minderr, int ipart, 
    double* h_rfi
) {}