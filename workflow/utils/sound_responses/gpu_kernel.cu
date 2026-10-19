// Unweighted float64 Bayesian blocks, matching the hepstats 0.9.2 recurrence.
// One persistent CUDA block per fit; Configured penalties share each fitness value.
// Prefix counts are integer-valued doubles (<2^53), hence their differences
// exactly reproduce the reference's reverse cumulative counts.
extern "C" __global__ void bayesian_blocks_batch(
    const double* lengths, const double* prefix, const double* penalties,
    const long long* offsets, const int* sizes, const long long total,
    double* best, int* last, int* invalid)
{
    const int task = blockIdx.x;
    const int tid = threadIdx.x;
    const int lane = tid % 32;
    const int warp = tid / 32;
    const int nw = blockDim.x / 32;
    const int n = sizes[task];
    const long long base = offsets[task];
    const long long edgebase = base + task;
    __shared__ double warpval[N_PENALTIES * 32];
    __shared__ int warparg[N_PENALTIES * 32];
    for (int r = 0; r < n; ++r) {
        double val[N_PENALTIES]; int arg[N_PENALTIES];
        #pragma unroll
        for (int k=0; k<N_PENALTIES; ++k) { val[k] = -1.0/0.0; arg[k] = 2147483647; }
        for (int j=tid; j<=r; j+=blockDim.x) {
            double width = lengths[edgebase+j] - lengths[edgebase+r+1];
            if (width <= 0.0) { atomicExch(invalid+task, 1); continue; }
            double count = prefix[edgebase+r+1] - prefix[edgebase+j];
            double fitness = count * log(count/width);
            #pragma unroll
            for (int k=0; k<N_PENALTIES; ++k) {
                double score = fitness - penalties[k*total+base+r];
                if (j>0) score += best[k*total+base+j-1];
                if (score > val[k] || (score == val[k] && j < arg[k])) {
                    val[k] = score; arg[k] = j;
                }
            }
        }
        #pragma unroll
        for (int k=0; k<N_PENALTIES; ++k) {
            for (int shift=16; shift>0; shift/=2) {
                double v = __shfl_down_sync(0xffffffff, val[k], shift);
                int a = __shfl_down_sync(0xffffffff, arg[k], shift);
                if (lane+shift < 32 && (v > val[k] || (v == val[k] && a < arg[k]))) {
                    val[k] = v; arg[k] = a;
                }
            }
            if (lane==0) { warpval[k*32+warp]=val[k]; warparg[k*32+warp]=arg[k]; }
        }
        __syncthreads();
        if (warp==0) {
            #pragma unroll
            for (int k=0; k<N_PENALTIES; ++k) {
                double v = lane<nw ? warpval[k*32+lane] : -1.0/0.0;
                int a = lane<nw ? warparg[k*32+lane] : 2147483647;
                for (int shift=16; shift>0; shift/=2) {
                    double w = __shfl_down_sync(0xffffffff,v,shift);
                    int b = __shfl_down_sync(0xffffffff,a,shift);
                    if (lane+shift < 32 && (w > v || (w == v && b < a))) {v=w; a=b;}
                }
                if (lane==0) { best[k*total+base+r]=v; last[k*total+base+r]=a; }
            }
        }
        __syncthreads();
        if (invalid[task]) return; // Uniform across block after the barrier.
    }
}
