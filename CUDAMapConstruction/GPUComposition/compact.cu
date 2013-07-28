/*
 * Diana-Andreea Popescu, EPFL & CERN, Switerland. All rights reserved.
 */

#include "Comp.cuh"

template <int BLOCK_SIZE, int NVARS>
__global__ void stream_compactation(int *exp_input,
			            double *coeff_input,
				    unsigned int *stencil,
                                    unsigned int *indices,
                                    unsigned int n,
				    unsigned int final_size,
				    int *exp_result,
                                    double *coeff_result)
{
    unsigned int tx = threadIdx.x;
    unsigned int tbx = blockIdx.x * BLOCK_SIZE + tx;
    const unsigned int grid_size = gridDim.x * BLOCK_SIZE;
    double coeff = 0; 
    __shared__ int exps[NVARS * BLOCK_SIZE];
   

    while (tbx < n)
    {
	coeff = coeff_input[tbx];
	for (int k = 0; k < NVARS; k++)
		exps[tx + k * BLOCK_SIZE] = exp_input[tbx + k * n];
	
	if(stencil[tbx] == 1) {
	    coeff_result[indices[tbx]] = coeff;
	    for (int k = 0; k < NVARS; k ++)
		exp_result[indices[tbx] + k * final_size] = /*exp_input[tbx + k * n];*/ exps[tx + k * BLOCK_SIZE];
			
	}
	tbx += grid_size;
    }
}


template <int NVARS>
__global__ void stream_compactation_key(unsigned long long *exp_key_input,
                                    double *coeff_input,
                                    unsigned int *stencil,
                                    unsigned int *indices,
                                    unsigned int n,
                                    unsigned int final_size,
                                    unsigned long long *exp_result,
                                    double *coeff_result)
{
    unsigned int tbx = blockIdx.x * blockDim.x + threadIdx.x;
    const unsigned int grid_size = gridDim.x * blockDim.x;
    double coeff = 0;
    unsigned long long key = 0;

    while (tbx < n)
    {
        coeff = coeff_input[tbx];
	key = exp_key_input[tbx];
        if(stencil[tbx] == 1) {
            coeff_result[indices[tbx]] = coeff;
            exp_result[indices[tbx]] = key;
        }
        tbx += grid_size;
    }
}


void compact_pol(int *exp_input,
            double *coeff_input,
            uint *stencil,
            uint *indices,
            int n,
            int final_size,
            int *exp_result,
            double *coeff_result,
        cudaStream_t stream){

        uint dimGrid = n/THREADBLOCK_SIZE;
        if (n % THREADBLOCK_SIZE != 0)
                dimGrid ++;
        stream_compactation<THREADBLOCK_SIZE, NRVARS><<<dimGrid, THREADBLOCK_SIZE, 0, stream>>>(exp_input, coeff_input, stencil, indices, n, final_size, exp_result, coeff_result);
}

void compact_pol_key(unsigned long long *exp_input,
            double *coeff_input,
            uint *stencil,
            uint *indices,
            int n,
            int final_size,
            unsigned long long *exp_result,
            double *coeff_result,
        cudaStream_t stream){

        uint dimGrid = n/THREADBLOCK_SIZE;
        if (n % THREADBLOCK_SIZE != 0)
                dimGrid ++;
        stream_compactation_key<NRVARS><<<dimGrid, THREADBLOCK_SIZE, 0, stream>>>(exp_input, coeff_input, stencil, indices, n, final_size, exp_result, coeff_result);
}


