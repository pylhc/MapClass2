/*
 * Diana-Andreea Popescu, EPFL & CERN, Switerland. All rights reserved.
 */


#include <stdio.h>
#include <assert.h>
#include <math.h>

// CUDA runtime
#include <cuda_runtime.h>
#include <cuda.h>
#include <device_launch_parameters.h>
#include "Comp.cuh"
// Helper functions and utilities to work with CUDA
#include <helper_functions.h>
#include <thrust/device_ptr.h>
#include <thrust/reduce.h>
#include <thrust/sort.h>
#include <thrust/remove.h>
#include <omp.h>
#include <vector>
#include <list>
#include <algorithm>
#include <numeric>

using namespace std;

inline
void checkCuda(cudaError_t result)
{
  if (result != cudaSuccess) {
    fprintf(stderr, "CUDA Runtime Error: %s\n",
            cudaGetErrorString(result));
        exit(EXIT_FAILURE);
  }
}


unsigned int nextPow2(unsigned int x)
{
    --x;
    x |= x >> 1;
    x |= x >> 2;
    x |= x >> 4;
    x |= x >> 8;
    x |= x >> 16;
    return ++x;
}


struct is_order_less
{
    __host__ __device__
    bool operator() (const int x)
    {
	    return (x < 0);
    }
};

extern "C"
int composeOnGPU(vector<double> input_coeff, vector<list<int> > terms, int inputSize, 
		 vector<int*> other_exp, vector<double*> other_coeff, vector<uint> otherSize, int order,
		int* final_exponents, double* final_coeffs)
{
  double inits = omp_get_wtime();
    int nr_terms = terms.size();
    int nr_functions = otherSize.size();
    cudaStream_t streams[nr_terms];
    for (int i = 0; i < nr_terms; ++i) 
	    cudaStreamCreate(&streams[i]); 
    //duplicate the input polynoms for each stream
    //allocate device memory for input
    ////////////////////////////////////////////
    int *e_input;
    double *c_input;
    unsigned int mem_size_exp_input = 0;
    unsigned int mem_size_coeff_input = 0;
    int all_size = accumulate(otherSize.begin(), otherSize.end(), 0);
    mem_size_exp_input = NRVARS * all_size * nr_terms * sizeof(int);
    mem_size_coeff_input = all_size * nr_terms * sizeof(double);

    checkCuda(cudaMalloc((void **) &e_input, mem_size_exp_input));
    checkCuda(cudaMalloc((void **) &c_input, mem_size_coeff_input));
	
    ////////////////////////////////////////////
    //determine maximum size to use for allocation 
    //allocate memory for multiplication result, aux result, stencil, indices on each stream
    unsigned int max_size = *(max_element(otherSize.begin(), otherSize.end()));
    unsigned int result_size = 8 * max_size * max_size;
    unsigned int mem_size_result_exp = result_size * NRVARS * sizeof(int) * nr_terms;
    unsigned int mem_size_result_coeff = result_size  * sizeof(double) * nr_terms;
    int *e_result, *e_aux;
    double *c_result, *c_aux;
    checkCuda(cudaMalloc((void **) &e_result, mem_size_result_exp));
    checkCuda(cudaMalloc((void **) &c_result, mem_size_result_coeff));
    checkCuda(cudaMalloc((void **) &e_aux, mem_size_result_exp));
    checkCuda(cudaMalloc((void **) &c_aux, mem_size_result_coeff));
   
    //determine next power of 2 for stencil for stream compactation
    unsigned int max_size_pow2 = nextPow2(max_size * max_size);
    unsigned int stencil_size = 8 * max_size_pow2;
    unsigned int mem_size_stencil = stencil_size * sizeof(unsigned int) * nr_terms;
    unsigned int* stencil, *indices;
    checkCuda(cudaMalloc((void **) &stencil, mem_size_stencil));
    checkCuda(cudaMalloc((void **) &indices, mem_size_stencil));
    //alloc memoy for keys
    unsigned long long *e_keys, *keys;
    unsigned int mem_size_keys = result_size * sizeof(unsigned long long) * nr_terms;
    checkCuda(cudaMalloc((void **) &e_keys, mem_size_keys));
    checkCuda(cudaMalloc((void **) &keys, mem_size_keys));
    ///////////////////////////////////////////////
    //copy data to GPU
    for (int i = 0; i < nr_terms; ++i){
	    //compute offset
	    int offset_exp = i * all_size * NRVARS;
	    int offset_coeff = i * all_size;
	    for (int j = 0; j < nr_functions; ++j){ 
		    //compute offset
		    if (j != 0){
			    offset_exp += otherSize[j - 1] * NRVARS; 
			    offset_coeff += otherSize[j - 1];
		    }  
		    int mem_size_exp_j = otherSize[j] * NRVARS * sizeof(int);
		    int mem_size_coeff_j = otherSize[j] * sizeof(double);
		    checkCuda(cudaMemcpyAsync(e_input + offset_exp, other_exp[j], 
			    mem_size_exp_j, cudaMemcpyHostToDevice, streams[i]));
		    checkCuda(cudaMemcpyAsync(c_input + offset_coeff, other_coeff[j], mem_size_coeff_j, 
			    cudaMemcpyHostToDevice, streams[i]));
	    }	
    }
    //cout << result_size << " " << stencil_size << endl;
    ////////////////////////////////////////////
    //compute number of multiplications 		
    int iterations = -nr_terms;
    for (vector<list<int> >:: const_iterator it = terms.begin(); it != terms.end(); ++it)
	    iterations += it->size();
    unsigned int index_result_exp, index_result_coeff, index_input_exp, index_input_coeff;
    vector<pair<int*, double*> > p_result;
    vector<pair<int*, double*> > p_aux;
    for (int i = 0; i < nr_terms; i ++){
	index_result_exp = i * result_size * NRVARS;
	index_result_coeff = i * result_size;
      	p_result.push_back(make_pair(e_result + index_result_exp, c_result + index_result_coeff));
    } 
     
    vector<int> single_terms;
    uint* sizes_aux;
    checkCuda(cudaMallocHost((void **) &sizes_aux, sizeof(uint) * nr_terms));
    for (int i = 0; i < nr_terms; i ++){
    	index_input_exp = i * all_size * NRVARS;
	index_input_coeff = i * all_size;
	int index = terms[i].front();
	terms[i].pop_front();
	if (terms[i].empty())
	  	single_terms.push_back(i);
	for (int k = 1; k <= index; k ++){
		 index_input_exp += otherSize[k - 1] * NRVARS;
		 index_input_coeff += otherSize[k - 1];
	}
	sizes_aux[i] = otherSize[index];
	p_aux.push_back(make_pair(e_input + index_input_exp, c_input + index_input_coeff));
    }
    
    vector<bool> firstmult(nr_terms, true);
    vector<uint> dimRes(nr_terms, 0);
    vector<uint> is_zero(nr_terms, 0);
    double inite = omp_get_wtime();
    cout << "init + copy " << 1000 * (inite - inits) << endl;
    double start = omp_get_wtime();
    for (int i = 0; i < single_terms.size(); i ++) 
      	get_keys(keys + single_terms[i] * result_size, p_aux[single_terms[i]].first, sizes_aux[single_terms[i]], streams[single_terms[i]]);
	
    unsigned int dimResult, dimInput, dimAux;
   
    while (iterations > 0){
	    //launch kernels from each term
           #pragma omp parallel for private(index_input_exp, index_input_coeff, index_result_exp, index_result_coeff, dimInput, dimAux, dimResult) shared(iterations, indices, terms, dimRes, p_aux, p_result, stencil, order, streams, input_coeff) num_threads(32)  
	for (int i = 0; i < nr_terms; i ++){
		    if (!terms[i].empty()) {
			int index = terms[i].front();
			index_input_exp = i * all_size * NRVARS;
			index_input_coeff = i * all_size;
			for (int k = 1; k <= index; k ++){
				index_input_exp += otherSize[k - 1] * NRVARS;
				index_input_coeff += otherSize[k - 1];
			}
		       
			//cudaStreamSynchronize(streams[i]);
			dimInput = otherSize[index];
			dimAux = sizes_aux[i];
			dimResult = dimInput * dimAux;
			dimRes[i] = dimResult;
			if (terms[i].size() == 1){ 
			 	 multiply_truncate_key(p_result[i].first, e_keys + i * result_size,
			    		p_aux[i].first, e_input + index_input_exp,
			    		p_result[i].second, p_aux[i].second, c_input + index_input_coeff, 
				  	dimResult, dimAux, dimInput, order, stencil + i * stencil_size, input_coeff[i], streams[i]); 
			} 
			else {	
				multiply_truncate(p_result[i].first, p_aux[i].first, e_input + index_input_exp,
			 	p_result[i].second, p_aux[i].second, c_input + index_input_coeff, 
				 dimResult, dimAux, dimInput, order, stencil + i * stencil_size, streams[i]); 	 			    
			}
			if (stencil_size <= 4096)
				scanExclusiveShort(indices + i * stencil_size, stencil + i * stencil_size, stencil_size, streams[i]);
			else  scanExclusiveLarge(indices + i * stencil_size, stencil + i * stencil_size, stencil_size, streams[i]);
			
			checkCuda(cudaMemcpyAsync(&sizes_aux[i], indices + i * stencil_size + dimResult, sizeof(uint), cudaMemcpyDeviceToHost, streams[i]));
				cudaStreamSynchronize(streams[i]);	
		
                        if (firstmult[i]){
				firstmult[i] = false;
				index_result_exp = i * result_size * NRVARS;
				index_result_coeff = i * result_size;
				p_aux[i] = make_pair(e_aux + index_result_exp, c_aux + index_result_coeff);	
			} 
			//cout << sizes_aux[i] << endl;
			//kernel for truncation
			//	uint dimGrid = dimRes[i]/THREADBLOCK_SIZE;
                        //if (dimRes[i] % THREADBLOCK_SIZE != 0)
                        //        dimGrid ++;
			
			if (terms[i].size() == 1) 
			  compact_pol_key(e_keys + i * result_size, p_result[i].second, stencil + i * stencil_size, indices + i * stencil_size, dimResult, sizes_aux[i], keys + i * result_size, p_aux[i].second, streams[i]); 
			  //stream_compactation_key<NRVARS><<<dimGrid, THREADBLOCK_SIZE, 0, streams[i]>>>(e_keys + i * result_size, p_result[i].second, stencil + i * stencil_size, indices + i * stencil_size, dimResult, sizes_aux[i], keys + i * result_size, p_aux[i].second); 
			else 
			  	compact_pol(p_result[i].first, p_result[i].second, stencil + i * stencil_size, indices + i * stencil_size, 
			  		dimResult, sizes_aux[i], p_aux[i].first, p_aux[i].second, streams[i]);
			  //stream_compactation<THREADBLOCK_SIZE, NRVARS><<<dimGrid, THREADBLOCK_SIZE, 0, streams[i]>>>(p_result[i].first, p_result[i].second, stencil + i * stencil_size, indices + i * stencil_size, dimResult, sizes_aux[i], p_aux[i].first, p_aux[i].second);
			terms[i].pop_front();
			      	#pragma omp critical
			iterations --;
		
	     	}	
	    }	    
    }
    double endt = omp_get_wtime();
    cout << "iter " << 1000 * (endt - start) << endl;
    cudaThreadSynchronize();
    double startcopy = omp_get_wtime();
    uint size_keys = 0;
    for (int i = 0; i < nr_terms; i ++)
    	size_keys += sizes_aux[i];
  	
    cout << "total " << size_keys << endl;
    

    thrust::device_ptr<unsigned long long> e_keys_ptr(e_keys);
    thrust::device_ptr<unsigned long long> keys_ptr(keys);
    int index_keys = 0;
    thrust::device_ptr<double> c_result_ptr(c_result);
    thrust::device_ptr<double> c_aux_ptr(c_aux);
   
    for (int i = 0; i < nr_terms; i ++){
      	thrust::copy(keys_ptr + i * result_size, keys_ptr + i * result_size + sizes_aux[i], e_keys_ptr + index_keys);
    	thrust::copy(c_aux_ptr + i * result_size, c_aux_ptr + i * result_size + sizes_aux[i], c_result_ptr + index_keys);
	index_keys += sizes_aux[i];
    }
    
    double endcopy = omp_get_wtime();
    cout << "copy " << 1000 * (endcopy - startcopy) << endl;
    double startsort = omp_get_wtime();
    thrust::sort_by_key(e_keys_ptr, e_keys_ptr + size_keys, c_result_ptr);
    double endsort = omp_get_wtime();
    cout << "sort " << 1000 * (endsort - startsort) << endl;
    //reduce by key 
    double startreduce = omp_get_wtime();
    thrust::pair<thrust::device_ptr<unsigned long long>, thrust::device_ptr<double> > end;
    end = thrust::reduce_by_key(e_keys_ptr, e_keys_ptr + size_keys, c_result_ptr, keys_ptr, c_aux_ptr); 
    size_keys = end.first - keys_ptr;
    double endreduce = omp_get_wtime();
    cout << "reduce " << 1000 * (endreduce - startreduce) << endl;
    cout << "size " << size_keys << endl;
    get_exponents(e_result, keys, size_keys);		
    //coeffs are in c_aux_ptr
    for (int i = 0; i < nr_terms; ++i)
	    cudaStreamDestroy(streams[i]);

    checkCuda(cudaMemcpy(final_exponents, e_result,
                            size_keys * sizeof(int) * NRVARS, cudaMemcpyDeviceToHost));
    checkCuda(cudaMemcpy(final_coeffs, c_aux, size_keys * sizeof(double),
                            cudaMemcpyDeviceToHost));

 	
    //free memory
    cudaFree(e_input);
    cudaFree(c_input);
    cudaFree(e_result);
    cudaFree(c_result);
    cudaFree(e_aux);
    cudaFree(c_aux);
    cudaFree(e_keys);
    cudaFree(keys);
    cudaFree(stencil);
    cudaFree(indices);
 	
    return size_keys; 
}
