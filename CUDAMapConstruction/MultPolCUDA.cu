/**
 * Copyright 2013 Diana-Andreea Popescu, EPFL & CERN, Switzerland.  All rights reserved.
 *
 */

// System includes
#include <stdio.h>
#include <assert.h>
#include <math.h>

// CUDA runtime
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_launch_parameters.h>

// Helper functions and utilities to work with CUDA
#include <helper_functions.h>
#if defined(__CUDA_ARCH__) && (__CUDA_ARCH__ < 200)
	#define printf(f, ...) ((void)(f, __VA_ARGS__),0)
#endif
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/remove.h>

inline
void checkCuda(cudaError_t result)
{
    if (result != cudaSuccess) {
	    fprintf(stderr, "CUDA Runtime Error: %s\n", 
	cudaGetErrorString(result));
    exit(EXIT_FAILURE);
  }
}


struct is_order_less
{
   __host__ __device__
   bool operator() (const int x)
   {
       return (x < 0);
   }
};

template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
computeResultTerms(int *exp_C, unsigned long long *exp_keys, int *exp_A, int *exp_B,
					double *coeff_C, double *coeff_A, double *coeff_B, 					
					int nC, int nA, int nB,
				       	int order, int* stencil)
{
    // Block index
    int bx = blockIdx.x;
    int by = blockIdx.y;

    // Thread index
    int tx = threadIdx.x;
    int ty = threadIdx.y;
	
    int tbx0 = tx + bx * BLOCK_SIZE_X;
    int tby0 = ty + by * BLOCK_SIZE_Y;

    int tbx = tbx0;
    int tby = tby0;

    int offsetx = BLOCK_SIZE_X * gridDim.x;
    int offsety = BLOCK_SIZE_Y * gridDim.y;
	
    __shared__ int Aes[BLOCK_SIZE_X * NVARS];
    __shared__ double Acs[BLOCK_SIZE_X];

    __shared__ int Bes[BLOCK_SIZE_Y * NVARS];
    __shared__ double Bcs[BLOCK_SIZE_Y];

    int Cexp[NVARS];
    double Ccoeff = 0;
    int c = 0;
    int sum = 0;
    unsigned long long ekey = 0;

    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		  
		    Ccoeff = Acs[tx] * Bcs[ty];
		    sum = 0;
		    ekey = 0;
		    c = nA * tby + tbx;
		    coeff_C[c] = Ccoeff;
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) 
		    {
			Cexp[k] = Aes[tx + k * BLOCK_SIZE_X] + Bes[ty + k * BLOCK_SIZE_Y];
			exp_C[c + k * nC] = Cexp[k]; 
			ekey = MAX_EXP * ekey + Cexp[k];
			sum += Cexp[k];
		    }
		    //see it has to be removed or not
		    stencil[c] = order - sum;
		    exp_keys[c] = ekey;	
		    //update index
		    tbx += offsetx; 
	    }
	    //update index
	    tby += offsety;
	    //reset 
	    tbx = tbx0;
    }  
}

template <int NVARS, int MAX_EXP> __global__ void
computeExponents(int *exp_C, unsigned long long *exp_keys, int nC)
{
    // Block index
    int bx = blockIdx.x;
    // Thread index
    int tx = threadIdx.x;
    int tbx = tx + bx * blockDim.x;
    int offset = blockDim.x * gridDim.x;
    unsigned long long key = 0, kd = 0;
    while (tbx < nC) {
	key = exp_keys[tbx];
#pragma unroll
	for (int k = NVARS - 1; k >= 0; k--) {
		kd = key/MAX_EXP;
		exp_C[tbx + k * nC] = key - kd * MAX_EXP; 
		key = kd;
	}
	tbx += offset;
    }
}

extern "C"
int* allocExponentsMemory(int dimA, int nvars){
    int size_A = dimA * nvars;
    int mem_size_exp_A = sizeof(int) * size_A;
    int *exp_A = NULL;
    checkCuda(cudaMallocHost((void **)&exp_A, mem_size_exp_A));
    return exp_A;
}

extern "C"
double* allocCoefficientsMemory(int dimA){
    int mem_size_coeff_A = sizeof(double) * dimA;
    double *coeff_A = NULL;
    checkCuda(cudaMallocHost((void **)&coeff_A, mem_size_coeff_A));
    return coeff_A;
}

extern "C"
void freeMemory(int* exp_C, double* coeff_C){
  cudaFree(exp_C);
  cudaFree(coeff_C);
}

extern "C"
int multPol(int* exp_A, double* coeff_A, int* exp_B, double* coeff_B, int* exp_C, double* coeff_C, 	 
	     int &dimA, int &dimB, 
	     int &order, int &nvars){

    int block_size_x = 16;
    int block_size_y = 16;
    int size_A = dimA * nvars;
    int mem_size_exp_A = sizeof(int) * size_A;
    int mem_size_coeff_A = sizeof(double) * dimA;

    int size_B = dimB * nvars;
    int mem_size_exp_B = sizeof(int) * size_B;
    int mem_size_coeff_B = sizeof(double) * dimB;

    // Allocate device memory
    int *e_A, *e_B, *e_C;
    double *c_A, *c_B, *c_C;
    double *final_coeff_C;
    unsigned long long *e_keys_C;
    unsigned long long *final_keys_C;

    int dimC = dimA * dimB;
    int size_C = dimA * dimB * nvars;
    int mem_size_exp_C = size_C * sizeof(int);
    int mem_size_keys_C = dimC * sizeof(unsigned long long); 
    int mem_size_coeff_C = sizeof(double) * dimC;

    checkCuda(cudaMalloc((void **) &e_A, mem_size_exp_A));
    checkCuda(cudaMalloc((void **) &c_A, mem_size_coeff_A));

    checkCuda(cudaMalloc((void **) &e_B, mem_size_exp_B));
    checkCuda(cudaMalloc((void **) &c_B, mem_size_coeff_B));

    checkCuda(cudaMalloc((void **) &e_C, mem_size_exp_C));
    checkCuda(cudaMalloc((void **) &c_C, mem_size_coeff_C));

    checkCuda(cudaMalloc((void **) &final_coeff_C, mem_size_coeff_C));
    checkCuda(cudaMalloc((void **) &e_keys_C, mem_size_keys_C));
    checkCuda(cudaMalloc((void **) &final_keys_C, mem_size_keys_C));

    //stencil used for truncation
    int *stencil = NULL;
    int mem_size_stencil = sizeof(int) * dimC;
    checkCuda(cudaMalloc((void **) &stencil, mem_size_stencil));

     // copy host memory to device
    checkCuda(cudaMemcpy(e_A, exp_A, mem_size_exp_A, cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(c_A, coeff_A, mem_size_coeff_A, cudaMemcpyHostToDevice));

    checkCuda(cudaMemcpy(e_B, exp_B, mem_size_exp_B, cudaMemcpyHostToDevice));
    checkCuda(cudaMemcpy(c_B, coeff_B, mem_size_coeff_B, cudaMemcpyHostToDevice));

    dim3 threads(block_size_x, block_size_y);
    int limitA = dimA % threads.x;
    int x1 = 1;
    if (limitA == 0)
      	x1 = 0;
    int limitB = dimB % threads.y;
    int y1 = 1;
    if (limitB == 0)
	y1 = 0;
    dim3 grid(dimA / threads.x + x1, dimB / threads.y + y1);
    //printf("Dims %d %d\n", dimA, dimB);
    //printf("Grid  %d %d order %d\n", grid.x, grid.y, order); 	
    //printf("Computing result using CUDA Kernel...\n");    

    computeResultTerms<16, 16, 6, 100><<< grid, threads >>>(e_C, e_keys_C, 
			e_A, e_B, c_C, c_A, c_B, dimC, dimA, dimB, order, stencil);
   
    thrust::device_ptr<unsigned long long> keys_C_dev(e_keys_C);
    thrust::device_ptr<double> coeffs_C_dev(c_C);
    thrust::device_ptr<unsigned long long> final_keys_C_dev(final_keys_C);
    thrust::device_ptr<double> final_coeff_C_dev(final_coeff_C);
    //remove terms which have order higher than the input order
    thrust::device_ptr<int> stencil_dev(stencil);
    thrust::device_ptr<double> end_coeffs_dev = thrust::remove_if(coeffs_C_dev, coeffs_C_dev + dimC, stencil_dev, is_order_less());
    thrust::device_ptr<unsigned long long> end_keys_dev = thrust::remove_if(keys_C_dev, keys_C_dev + dimC, stencil_dev, is_order_less());
    //sort
    thrust::sort_by_key(keys_C_dev, end_keys_dev, coeffs_C_dev);
    thrust::pair<thrust::device_ptr<unsigned long long>, thrust::device_ptr<double> > end;
    //reduce
    end = thrust::reduce_by_key(keys_C_dev, end_keys_dev, coeffs_C_dev, final_keys_C_dev, final_coeff_C_dev); 	
    int sizeC = end.first - final_keys_C_dev; 
    dim3 threads_exp(block_size_x * block_size_x);
    if (sizeC % threads_exp.x == 0)
	x1 = 0;
    else 
	x1 = 1;
    dim3 grid_exp(sizeC/threads_exp.x + x1);

    computeExponents<6,100><<< grid_exp, threads_exp >>>(e_C, final_keys_C, sizeC);

    // Copy result from device to host
    checkCuda(cudaMemcpy(exp_C, e_C, sizeC * nvars * sizeof(int), cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(coeff_C, final_coeff_C, sizeC * sizeof(double), cudaMemcpyDeviceToHost));

    // Clean up memory
   
    cudaFree(e_A);
    cudaFree(e_B);
    cudaFree(e_C);
    cudaFree(c_A);
    cudaFree(c_B);
    cudaFree(c_C);

    cudaFree(e_keys_C);
    cudaFree(final_keys_C);
    cudaFree(final_coeff_C);
    cudaFree(stencil);

    return sizeC;

}
