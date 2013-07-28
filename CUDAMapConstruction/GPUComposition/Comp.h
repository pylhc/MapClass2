#include <stdio.h>
#include <assert.h>
#include <math.h>

// CUDA runtime
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_launch_parameters.h>

#define THREADBLOCK_SIZE 1024

#define BL_SIZE_X    16
#define BL_SIZE_Y    16
#define NRVARS       6
#define MAX_E        100

/*template<int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EPX>
__global__ void multiply_pols_truncate(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil);
template<int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EPX>
__global__ void multiply_pols_truncate2(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil);

template<int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EPX>
__global__ void multiply_pols_truncate_key(int *exp_C,  unsigned long long *keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil, double coeff);

template<int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EPX>
__global__ void multiply_pols_truncate_key2(int *exp_C,  unsigned long long *keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil, double coeff);
*/

void multiply_truncate(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil, cudaStream_t stream);
void multiply_truncate2(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil, cudaStream_t stream);


void multiply_truncate_key(int *exp_C,  unsigned long long *exp_keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil, double coef, cudaStream_t stream);

void multiply_truncate_key2(int *exp_C,  unsigned long long *exp_keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil, double coef, cudaStream_t stream);
 
void get_exponents(int *exp_C, unsigned long long *exp_keys, unsigned int nC);

void get_keys(unsigned long long *exp_keys, int* exp_C, unsigned int nC, cudaStream_t stream);

void compact_pol(int *exp_input,
            double *coeff_input,
            uint *stencil,
            uint *indices,
            int n,
            int final_size,
            int *exp_result,
            double *coeff_result,
        cudaStream_t stream);



void compact_pol_key(unsigned long long *exp_input,
            double *coeff_input,
            uint *stencil,
            uint *indices,
            int n,
            int final_size,
            unsigned long long *exp_result,
            double *coeff_result,
        cudaStream_t stream);

typedef unsigned int uint;

////////////////////////////////////////////////////////////////////////////////
// Implementation limits
////////////////////////////////////////////////////////////////////////////////
extern "C" const uint MAX_BATCH_ELEMENTS;
extern "C" const uint MIN_SHORT_ARRAY_SIZE;
extern "C" const uint MAX_SHORT_ARRAY_SIZE;
extern "C" const uint MIN_LARGE_ARRAY_SIZE;
extern "C" const uint MAX_LARGE_ARRAY_SIZE;

////////////////////////////////////////////////////////////////////////////////
// CUDA scan
////////////////////////////////////////////////////////////////////////////////

void scanExclusiveShort(
    uint *d_Dst,
    uint *d_Src,
    uint arrayLength,
    cudaStream_t stream
);

void scanExclusiveLarge(
    uint *d_Dst,
    uint *d_Src,
    uint arrayLength,
    cudaStream_t stream
);

