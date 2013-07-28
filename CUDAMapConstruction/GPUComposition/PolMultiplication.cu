/**
 * Copyright 2013 Diana-Andreea Popescu, EPFL & CERN, Switzerland.  All rights reserved.
 *
 */

// System includes
#include <stdio.h>
#include <assert.h>
#include <math.h>
#include "Comp.cuh"

// CUDA runtime
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#include <helper_cuda.h>

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


template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
computeResultTersmArbitrarySized(int *exp_C, unsigned long long *exp_keys, int *exp_A, int *exp_B,
	double *coeff_C, double *coeff_A, double *coeff_B,
	unsigned int nC, unsigned int nA, unsigned int nB)
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
    unsigned long long ekey = 0;
    double Ccoeff = 0;
    int c = 0;

    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		    // Ccoeff is used to store the coefficient of the term product
		    // that is computed by the thread
		    Ccoeff = Acs[tx] * Bcs[ty];
		    ekey = 0;
		    c = nA * tby + tbx;
		    coeff_C[c] = Ccoeff;
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) 
		    {
			Cexp[k] = Aes[tx + k * BLOCK_SIZE_X] + Bes[ty + k * BLOCK_SIZE_Y];
			exp_C[c + k * nC] = Cexp[k]; 
			ekey = MAX_EXP * ekey + Cexp[k];
		    }
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

template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
multiply_pols_truncate2(int *exp_C,  int *exp_A, int *exp_B,
					double *coeff_C, double *coeff_A, double *coeff_B, 					
					unsigned int nC, unsigned int nA, unsigned int nB,
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
    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		    // Ccoeff is used to store the coefficient of the term product
		    // that is computed by the thread
		    Ccoeff = Acs[tx] * Bcs[ty];
		    sum = 0;
		    c = nA * tby + tbx;
		    coeff_C[c] = Ccoeff;
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) 
		    {
			Cexp[k] = Aes[tx + k * BLOCK_SIZE_X] + Bes[ty + k * BLOCK_SIZE_Y];
			exp_C[c + k * nC] = Cexp[k]; 
			//printf("%d %d %d %d\n", k, Cexp[k], Aes[tx + k * BLOCK_SIZE_X], Bes[ty + k * BLOCK_SIZE_Y]);
			sum += Cexp[k];
		    }
		    stencil[c] = order - sum;
			
		    //update index
		    tbx += offsetx; 
	    }
	    //update index
	    tby += offsety;
	    //reset 
	    tbx = tbx0;
    }  
}



template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
multiply_pols_truncate(int *exp_C,  int *exp_A, int *exp_B,
					double *coeff_C, double *coeff_A, double *coeff_B, 					
					unsigned int nC, unsigned int nA, unsigned int nB,
					int order, unsigned int* stencil)
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
    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		    // Ccoeff is used to store the coefficient of the term product
		    // that is computed by the thread
		    Ccoeff = Acs[tx] * Bcs[ty];
		    sum = 0;
		    c = nA * tby + tbx;
		    coeff_C[c] = Ccoeff;
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) 
		    {
			Cexp[k] = Aes[tx + k * BLOCK_SIZE_X] + Bes[ty + k * BLOCK_SIZE_Y];
			exp_C[c + k * nC] = Cexp[k]; 
			sum += Cexp[k];
		    }
		    if (sum <= order) {
		    	stencil[c] = 1;
		    }
		    else stencil[c] = 0;
			
		    //update index
		    tbx += offsetx; 
	    }
	    //update index
	    tby += offsety;
	    //reset 
	    tbx = tbx0;
    }  
}

template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
multiply_pols_truncate_key(int *exp_C, unsigned long long *exp_keys, int *exp_A, int *exp_B,
			   double *coeff_C, double *coeff_A, double *coeff_B, 					
			   unsigned int nC, unsigned int nA, unsigned int nB,
			   int order, unsigned int* stencil, double transform)
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
    unsigned long long ekey = 0;
    double Ccoeff = 0;
    int c = 0;
    int sum = 0;
    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		    // Ccoeff is used to store the coefficient of the term product
		    // that is computed by the thread
		    Ccoeff = Acs[tx] * Bcs[ty] * transform;
		    ekey = 0;
		    sum = 0;
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
		    if (sum <= order)
		    	stencil[c] = 1;
		    else stencil[c] = 0;
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

template <int BLOCK_SIZE_X, int BLOCK_SIZE_Y, int NVARS, int MAX_EXP> __global__ void
multiply_pols_truncate_key2(int *exp_C, unsigned long long *exp_keys, int *exp_A, int *exp_B,
			   double *coeff_C, double *coeff_A, double *coeff_B, 					
			   unsigned int nC, unsigned int nA, unsigned int nB,
			   int order, int* stencil, double transform)
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
    unsigned long long ekey = 0;
    double Ccoeff = 0;
    int c = 0;
    int sum = 0;
    while (tby < nB) {
	    while (tbx < nA) {
		    Acs[tx] = coeff_A[tbx];
		    Bcs[ty] = coeff_B[tby];
#pragma unroll
		    for (int k = 0; k < NVARS; ++k) {
			    Aes[tx + k * BLOCK_SIZE_X] = exp_A[tbx + k * nA]; 
			    Bes[ty + k * BLOCK_SIZE_Y] = exp_B[tby + k * nB];
		    }
		    // Ccoeff is used to store the coefficient of the term product
		    // that is computed by the thread
		    Ccoeff = Acs[tx] * Bcs[ty] * transform;
		    ekey = 0;
		    sum = 0;
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
get_keys_from_exponents(unsigned long long *exp_keys, int *exp_A, unsigned int nA)
{
    // Block index
    int bx = blockIdx.x;

    // Thread index
    int tx = threadIdx.x;
	
    int tbx = tx + bx * blockDim.x;
    int offset = blockDim.x * gridDim.x;	
    
    unsigned long long ekey = 0;
    
    while (tbx < nA) {
	ekey = 0;	  
#pragma unroll
	for (int k = 0; k < NVARS; ++k) {
		ekey = MAX_EXP * ekey + exp_A[tbx + k * nA]; 
		exp_keys[tbx] = ekey;
	}
	//update index
        tbx += offset;
     }
}

template <int NVARS, int MAX_EXP> __global__ void
get_exponents_from_key(int *exp_C, unsigned long long *exp_keys, unsigned int nC)
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
		//	printf("%llu \n", key);
#pragma unroll
		for (int k = NVARS - 1; k >= 0; k--) {
		 	kd = key/MAX_EXP;
			exp_C[tbx + k * nC] = key - kd * MAX_EXP;
			key = kd;
		  	//exp_C[tbx + k * nC] = key % MAX_EXP;
		  	//key /= MAX_EXP;
			
		}
		tbx += offset;
	}
}


/*void initPol(unsigned int *exps, unsigned int dim, double *coeffs, unsigned int nvars)
{
    for (unsigned int i = 0; i < dim; ++i)
    {
        for (unsigned int k = 0; k < nvars; ++k)
			exps[i + k * dim] = rand() % 20;
		coeffs[i] = 1;
    }
//	for (unsigned int i = 0; i < dim * nvars; ++i)
//		  printf("%d ", exps[i]);
}*/


void multiply_truncate(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil, cudaStream_t stream){
	dim3 threads(BL_SIZE_X, BL_SIZE_Y);
	int x1 = (nA % threads.x == 0) ? 0 : 1;
        int y1 = (nB % threads.y == 0) ? 0 : 1;

        dim3 grid(nA/BL_SIZE_X + x1, nB/BL_SIZE_Y + y1);
	multiply_pols_truncate<BL_SIZE_X, BL_SIZE_Y, NRVARS, MAX_E> <<< grid, threads, 0, stream >>>(exp_C, exp_A, exp_B, coeff_C, coeff_A, coeff_B, nC, nA, nB, order, stencil);
	getLastCudaError("multiply truncate execution FAILED\n");
}

void multiply_truncate2(int *exp_C,  int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil, cudaStream_t stream){
        dim3 threads(BL_SIZE_X, BL_SIZE_Y);
        int x1 = (nA % threads.x == 0) ? 0 : 1;
        int y1 = (nB % threads.y == 0) ? 0 : 1;

        dim3 grid(nA/BL_SIZE_X + x1, nB/BL_SIZE_Y + y1);
        multiply_pols_truncate2<BL_SIZE_X, BL_SIZE_Y, NRVARS, MAX_E> <<< grid, threads, 0, stream >>>(exp_C, exp_A, exp_B, coeff_C, coeff_A, coeff_B, nC, nA, nB, order, stencil);
        getLastCudaError("multiply truncate execution FAILED\n");
}

void multiply_truncate_key(int *exp_C,  unsigned long long *exp_keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, unsigned int* stencil, double coef, cudaStream_t stream){
        dim3 threads(BL_SIZE_X, BL_SIZE_Y);
	int x1 = (nA % threads.x == 0) ? 0 : 1;
    	int y1 = (nB % threads.y == 0) ? 0 : 1;

        dim3 grid(nA/BL_SIZE_X + x1, nB/BL_SIZE_Y + y1);
        multiply_pols_truncate_key<BL_SIZE_X, BL_SIZE_Y, NRVARS, MAX_E> <<< grid, threads, 0, stream >>>(exp_C, exp_keys,exp_A, exp_B, coeff_C, coeff_A, coeff_B, nC, nA, nB, order, stencil, coef);
	getLastCudaError("multiply truncate key execution FAILED\n");
}


void multiply_truncate_key2(int *exp_C,  unsigned long long *exp_keys, int *exp_A, int *exp_B,
              double *coeff_C, double *coeff_A, double *coeff_B,
              unsigned int nC, unsigned int nA, unsigned int nB,
              int order, int* stencil, double coef, cudaStream_t stream){
        dim3 threads(BL_SIZE_X, BL_SIZE_Y);
        int x1 = (nA % threads.x == 0) ? 0 : 1;
        int y1 = (nB % threads.y == 0) ? 0 : 1;

        dim3 grid(nA/BL_SIZE_X + x1, nB/BL_SIZE_Y + y1);
        multiply_pols_truncate_key2<BL_SIZE_X, BL_SIZE_Y, NRVARS, MAX_E> <<< grid, threads, 0, stream >>>(exp_C, exp_keys,exp_A, exp_B, coeff_C, coeff_A, coeff_B, nC, nA, nB, order, stencil, coef);
        getLastCudaError("multiply truncate key execution FAILED\n");
}

void get_exponents(int *exp_C, unsigned long long *exp_keys, unsigned int nC){
	dim3 threads(THREADBLOCK_SIZE);
        int x1 = (nC % threads.x == 0) ? 0 : 1;
	dim3 grid(nC/THREADBLOCK_SIZE + x1);
	get_exponents_from_key<NRVARS, MAX_E><<<grid, threads>>>(exp_C, exp_keys, nC);
}

void get_keys(unsigned long long *exp_keys, int* exp_C, unsigned int nC, cudaStream_t stream){
        dim3 threads(THREADBLOCK_SIZE);
        int x1 = (nC % threads.x == 0) ? 0 : 1;
        dim3 grid(nC/THREADBLOCK_SIZE + x1);
        get_keys_from_exponents<NRVARS, MAX_E><<<grid, threads, 0, stream>>>(exp_keys, exp_C, nC);
}

