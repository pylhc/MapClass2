#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include <complex>
#include <sstream>
#include <omp.h>
#include "Pol.h"
#include <cuda_runtime.h>
#include <cuda_profiler_api.h>
#include "Comp.h"
using namespace std;
using namespace boost;

#define MAX_SIZE	1024

extern "C"
int composeOnGPU(vector<double> input_coeff, vector<list<int> > terms, int inputSize,
                 vector<int*> other_exp, vector<double*> other_coeff, vector<uint> otherSizes, int order,
		int* final_exponents, double* final_coeffs);


template <typename Type>
struct pair_hash {
    std::size_t operator()(pair<Type, Type> const& c) const {
        std::size_t seed = 0;
        boost::hash_combine(seed, get<0>(c));
        boost::hash_combine(seed, get<1>(c));
        return seed;
    }
};

template  <class T> 
class Polmap {
	public:
		unordered_map<string, Polynom<T>> pols;
	public:
		Polmap(unordered_map<string, Polynom<T>> ps);
		Polmap(const Polmap &other);
		Polmap();
		Polmap& operator= (const Polmap &other);
		Polmap operator +(Polmap other);
		Polmap operator -(Polmap other);
		Polmap operator *(Polmap other);
		Polmap eval(string var, Polynom<T> other);
		Polmap parallelmult(Polmap other);
		Polmap operator ^(int n);
		Polynom<T> operator [](string s);
		Polmap& truncate();
		Polmap& truncate(int order);
		int getOrder();
		double getEps();
		void printpolmap();
		string getMap();
		Polmap compositionWithoutHash(Polmap other);
		Polmap compositionTermsParallel(Polmap other, int nb_threads);
};

template <class T> Polmap<T>::Polmap(unordered_map<string, Polynom<T>> ps) {
	pols = ps;
}

template <class T> Polmap<T>::Polmap() {
	
}

template <class T> Polmap<T>::Polmap(const Polmap &other) {
	pols = other.pols;
}

template <class T> Polmap<T>&:: Polmap<T>::operator=( const Polmap& other ) {
     if (this == &other)
        return *this;
     pols = other.pols;
     return *this;
}

template <class T> int Polmap<T> ::getOrder() {
	typename unordered_map<string, Polynom<T>>::const_iterator i = pols.begin();
	return i->second.order;
}


template <class T> double Polmap<T> ::getEps() {
	typename unordered_map<string, Polynom<T>>::const_iterator i = pols.begin();
	return i->second.eps;
}

template <class T> Polmap<T> Polmap<T>::operator+(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	for (typename unordered_map<string, Polynom<T>>::const_iterator i = other.pols.begin(); i != other.pols.end(); ++i) {
		result_pols[i->first] = this[i->first] + i->second;
	}
	return Polmap<T>(result_pols);
}

template <class T> Polmap<T> Polmap<T>::operator-(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	for (typename unordered_map<string, Polynom<T>>::const_iterator i = other.pols.begin(); i != other.pols.end(); ++i) {
		result_pols[i->first] = this[i->first] - i->second;
	}
	return Polmap<T>(result_pols);
}

template <class T> Polmap<T> Polmap<T>::eval(string var, Polynom<T> other) {
	
	for (typename unordered_map<string, Polynom<T>>::const_iterator i = pols.begin(); i != pols.end(); ++i) {
		pols[i->first] = pols[i->first].eval(var, other);
	}
	return *this;
}

/*template <class T> Polmap<T> Polmap<T>::operator*(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
		Polynom<T> p = itmap->second;
		Polynom<T> result = Polynom<T>(p.order, p.eps, p.vars, 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
			//construct pol for each var
			Polynom<T> res = Polynom<T>(p.order, p.eps, p.vars[0], i->first[0]);
			res = res.eval(p.vars[0], other.pols[p.vars[0]]); 
			for (int k = 1; k < p.vars.size(); k ++) {
				Polynom<T> subterm = Polynom<T>(p.order, p.eps, p.vars[k], i->first[k]); 
				subterm = subterm.eval(p.vars[k], other.pols[p.vars[k]]);
				res = res * subterm;
			}
			result = result + res * i->second; 
		} 
		result_pols[itmap->first] = result;
	}
	return Polmap<T>(result_pols);
} */

template <class T> Polmap<T> Polmap<T>::operator*(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	unordered_map<pair<string, string>, string, pair_hash<string>> temp;
	unordered_map<string, Polynom<T>> tempvalues;
	string lm = "0";
	tempvalues = other.pols;
	list<string> m;
	string t;
	for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
		Polynom<T> p = itmap->second;
		Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
		vector<int> zero(p.vars.size(), 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
			m.clear();
			Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, i->second);
			for (int k = 0; k < p.vars.size(); k ++) {
				typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
				if (found == other.pols.end())
					coef = coef * Polynom<T> (p.order, p.eps, p.vars[k], i->first[k]);  
				else
					for (int a = 0; a < i->first[k]; a ++)
						m.push_back(p.vars[k]);
			}
			if (m.size() > 0) {
				while (m.size() > 1) {
					string first = m.back();
					m.pop_back();
					string second = m.back();
					m.pop_back();
					pair<string, string> pp = make_pair (first, second);	
					unordered_map<pair<string, string>, string, pair_hash<string>>::const_iterator found = temp.find (pp);
					if (found != temp.end())
						t = found->second;
					else {
						int ilm = lexical_cast<int>(lm) + 1;
						lm = lexical_cast<string>(ilm);
						temp[pp] = lm;
						t = lm;
						Polynom<T> newv = tempvalues[get<0>(pp)] * tempvalues[get<1>(pp)]; 
						tempvalues[t] = newv;
					}	
					m.push_back(t);
				}
				newpol = newpol + coef * tempvalues[m.front()];
			}
			else 
				newpol = newpol + coef;
			 
		} 
		if (fabs(newpol.terms[zero]) < p.eps)
			newpol.terms.erase(zero);
		result_pols[itmap->first] = newpol;
	}
	return Polmap<T>(result_pols);
}

template <class T> Polmap<T> Polmap<T>::compositionWithoutHash(Polmap<T> other){
 	unordered_map<string, Polynom<T>> result_pols;
      	vector<string> m;
        for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
                Polynom<T> p = itmap->second;
		//do not use s and d in composition, as they do not change 
		if (itmap->first.compare("s") == 0 || itmap->first.compare("d") == 0){
		  	result_pols[itmap->first] = p;
		 	continue;
		}                          
                Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
		vector<int> zero(p.vars.size(), 0);
                for (typename unordered_map<vector<int>, double, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
                        m.clear();
                        Polynom<T> result = Polynom<T>(p.order, p.eps, p.vars, i->second);
                        for (int k = 0; k < p.vars.size(); k ++) {
                                typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
                                if (found == other.pols.end())
                                        result = result * Polynom<T> (p.order, p.eps, p.vars[k], i->first[k]);  
                                else
                                        for (int a = 0; a < i->first[k]; a ++)
					       m.push_back(p.vars[k]);
			} 	
			//cout << m.size() << endl;
			//	double start = omp_get_wtime();
			for (vector<string>::iterator k = m.begin(); k != m.end(); k ++)	        
				result = result * other.pols[*k];
			//	double end = omp_get_wtime();
			//	cout << end - start << " s" << endl;
                        newpol = newpol + result;
                } 
		if (fabs(newpol.terms[zero]) < p.eps)
			 newpol.terms.erase(zero);

                result_pols[itmap->first] = newpol;
        }

        return Polmap<T>(result_pols);

}

template <class T> Polmap<T> Polmap<T>::compositionTermsParallel(Polmap<T> other, int nb_threads){
 	unordered_map<string, Polynom<T>> result_pols;
	unordered_map<string, Polynom<T>> temp;
	temp = other.pols;
      	vector<string> m;
	//multimap<int, vector<int>> estimator;
        for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
                Polynom<T> p = itmap->second;
		//do not use s and d in composition, as they do not change 
		if (itmap->first.compare("s") == 0 || itmap->first.compare("d") == 0){
		  	result_pols[itmap->first] = p;
		 	continue;
		}                          
		//distribute the work to threads */
		Polynom<T>* newpols = new Polynom<T>[nb_threads];
		for (int i = 0; i < nb_threads; i ++)
			newpols[i] = Polynom<T>(p.order, p.eps, p.vars, 0);
		typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator it[p.terms.size()];
		int size = 0;
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) 
		  	it[size ++] = i;
		
#pragma omp parallel for firstprivate(temp, m) shared(newpols) schedule(dynamic)
		for (int i = 0; i < size; ++i) {
			
		  	int index = omp_get_thread_num();
                        m.clear();
                        Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, it[i]->second);
                        for (int k = 0; k < p.vars.size(); k ++) {
                                typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
                                if (found == other.pols.end())
                                        coef = coef * Polynom<T> (p.order, p.eps, p.vars[k], it[i]->first[k]);  
                                else
                                        for (int a = 0; a < it[i]->first[k]; a ++)
					       m.push_back(p.vars[k]);
			} 	
		 	if (m.size() > 0) {
                                while (m.size() > 1) {
                                        string first = m.back();
                                        m.pop_back();
                                        string second = m.back();
                                        m.pop_back();
                                        string pp = first + " " + second;       
                                        typename unordered_map<string, Polynom<T>>::const_iterator found = temp.find (pp);
                                        if (found == temp.end()) {
						Polynom<T> newv = temp[first] * temp[second]; 
					  	temp.insert(make_pair(pp, newv));
                                        }       
                                        m.push_back(pp);
                                }
				newpols[index] = newpols[index] + coef * temp[m.front()];
                        }
                        else 
			  	newpols[index] = newpols[index] + coef;
		
		}
		vector<int> zero(p.vars.size(), 0);
		Polynom<T> newpol = newpols[0];
		for (int i = 1; i < nb_threads; i ++)
		  	newpol = newpol + newpols[i];
		if (fabs(newpol.terms[zero]) < p.eps)
			newpol.terms.erase(zero);
		result_pols[itmap->first] = newpol;
		delete[] newpols;
        }

        return Polmap<T>(result_pols);

}

template <class T> Polmap<T> Polmap<T>::parallelmult(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	unordered_map<pair<string, string>, string, pair_hash<string>> temp;
	unordered_map<string, Polynom<T>> tempvalues;
	vector<pair<string, Polynom<T>>> vpols;
	string lm = "0";
	tempvalues = other.pols;
	for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) 
		vpols.push_back(make_pair(itmap->first, itmap->second));
	list<string> m;
	string t;
	#pragma omp parallel for firstprivate(lm, tempvalues, temp, m, t) shared(result_pols) schedule(static)
	for (int itmap = 0; itmap < vpols.size(); itmap ++) {
		Polynom<T> p = get<1>(vpols[itmap]);
		Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
		vector<int> zero(p.vars.size(), 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
			m.clear();
			Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, i->second);
			for (int k = 0; k < p.vars.size(); k ++) {
				typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
				if (found == other.pols.end())
					coef = coef * Polynom<T> (p.order, p.eps, p.vars[k], i->first[k]);  
				else
					for (int a = 0; a < i->first[k]; a ++)
						m.push_back(p.vars[k]);
			}
			if (m.size() > 0) {
				while (m.size() > 1) {
					string first = m.back();
					m.pop_back();
					string second = m.back();
					m.pop_back();
					pair<string, string> pp = make_pair (first, second);	
					unordered_map<pair<string, string>, string, pair_hash<string>>::const_iterator found = temp.find (pp);
					if (found != temp.end())
						t = found->second;
					else {
						int ilm = lexical_cast<int>(lm) + 1;
						lm = lexical_cast<string>(ilm);
						temp[pp] = lm;
						t = lm;
						Polynom<T> newv = tempvalues[get<0>(pp)] * tempvalues[get<1>(pp)]; 
						tempvalues[t] = newv;
					}	
					m.push_back(t);
				}
				newpol = newpol + coef * tempvalues[m.front()];
				
			}
			else 
				newpol = newpol + coef;
			 
		} 
		if (fabs(newpol.terms[zero]) < p.eps)
			newpol.terms.erase(zero);
		result_pols[get<0>(vpols[itmap])] = newpol;
	}
	return Polmap<T>(result_pols);
}


template <class T> Polmap<T> Polmap<T>::operator^(int n) {
	if (n == 0) {
		unordered_map<string, Polynom<T>> result_pols;
		vector<string> vars;
		for (typename unordered_map<string, Polynom<T>>::const_iterator i = pols.begin(); i != pols.end(); ++i) {
			vars.clear();
			vars.push_back(i->first);
			unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
			vector<int> v;
			v.push_back(1);
			result_terms[v] = 1;
			result_pols[i->first] = new Polynom<T>(getOrder(), getEps(), vars, result_terms);  
		}
		return Polmap<T>(result_pols);
	
	}
	Polmap<T> res = *this;
	Polmap<T> p = *this;
	for (int i = 1; i < n; i ++)
		res = res * p;
	return res;
}

template <class T> Polynom<T> Polmap<T>::operator[](string s) {
	return pols[s];
}

template <class T> Polmap<T>& Polmap<T>::truncate() {
	for (typename unordered_map<string, Polynom<T>>::iterator it=pols.begin(); it !=pols.end(); ++it)
		pols[it->first] = it->second.truncate();
	return *this;
}

template <class T> Polmap<T>& Polmap<T>::truncate(int o) {
	for (typename unordered_map<string, Polynom<T>>::iterator it=pols.begin(); it !=pols.end(); ++it)
		pols[it->first] = it->second.truncate(o);
	return *this;
}

template <class T> void Polmap<T>::printpolmap() {
	for (typename unordered_map<string, Polynom<T>>::iterator it=pols.begin(); it !=pols.end(); ++it) {
		cout << it->first << "=";
		it->second.printpol();
		cout<<endl;
	}
}

template <class T> string Polmap<T>::getMap() {
	string resultMap = "";	
	for (typename unordered_map<string, Polynom<T>>::iterator it=pols.begin(); it !=pols.end(); ++it) {
		resultMap += it->first + "|" + it->second.getPol() + "|";
	}
	return resultMap;
}

Polmap<double> composeGPU(Polmap<double> first, Polmap<double> other) {

  	//prepare second polmap 
	int nr_functions = first.pols.size();
	vector<int*> other_exp(nr_functions);
  	vector<double*> other_coeff(nr_functions);
  	vector<unsigned int> other_sizes(nr_functions);
	int *exp;
	double *coeff;
	vector<double> coeffs;
	list<int> term;
	vector<list<int>> input_terms;
       
  	for (typename unordered_map<string, Polynom<double>>::const_iterator itmap = other.pols.begin(); itmap != other.pols.end(); ++itmap) {
		//for each pol alloc memory and get raw representation
	//     	if (itmap->first.compare("s") == 0 || itmap->first.compare("d") == 0)
	//		continue;
		double ts = omp_get_wtime();
		Polynom<double> p = itmap->second;
		unsigned int size = p.terms.size();
		cudaMallocHost((void **)&exp, sizeof(int) * size * NRVARS);
		cudaMallocHost((void **)&coeff, sizeof(double) * size);
		p.getPolRaw(exp, coeff);
		if (itmap->first.compare("s") == 0){
		other_exp[3] = exp;
		other_coeff[3] = coeff;
		other_sizes[3] = size;
		}
		else  if (itmap->first.compare("d") == 0){
                other_exp[0] = exp;
                other_coeff[0] = coeff;
                other_sizes[0] = size;
                }
		else  if (itmap->first.compare("px") == 0){
                other_exp[1] = exp;
                other_coeff[1] = coeff;
                other_sizes[1] = size;
                }
		else  if (itmap->first.compare("py") == 0){
                other_exp[2] = exp;
                other_coeff[2] = coeff;
                other_sizes[2] = size;
                }
		else  if (itmap->first.compare("x") == 0){
                other_exp[4] = exp;
                other_coeff[4] = coeff;
                other_sizes[4] = size;
                }
		else  if (itmap->first.compare("y") == 0){
                other_exp[5] = exp;
                other_coeff[5] = coeff;
                other_sizes[5] = size;
                }
		double te = omp_get_wtime();
		cout << 1000 * (te - ts) << endl;

  	}
        unordered_map<string, Polynom<double>> result_pols;
	int *final_exponents;
	double *final_coeffs;
	
 	cudaMallocHost((void **)&final_exponents, sizeof(int) * MAX_SIZE * first.pols.size());
    	cudaMallocHost((void **)&final_coeffs, sizeof(double) * MAX_SIZE);
       //bool prof = true;
        for (typename unordered_map<string, Polynom<double>>::const_iterator itmap = first.pols.begin(); itmap != first.pols.end(); ++itmap) {

		Polynom<double> p = itmap->second;                
		//do not use s and d in composition, as they do not change 
		if (itmap->first.compare("s") == 0 || itmap->first.compare("d") == 0){
		  	result_pols[itmap->first] = p;
		 	continue;
		}  
		coeffs.clear();
		input_terms.clear();
		for (typename unordered_map<vector<int>, double, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i){
			coeffs.push_back(i->second);
			//prepare list of terms
			term.clear();
			for (int k = 0; k < p.vars.size(); k ++){
				for (int j = 0; j < i->first[k]; j ++)
                                	term.push_back(k);
			}
			input_terms.push_back(term);	
		}
		//double start = omp_get_wtime();
		//if (prof) cudaProfilerStart();

		int nb_terms = composeOnGPU(coeffs, input_terms, 4, other_exp, other_coeff, other_sizes, p.order,
			final_exponents, final_coeffs);                     
              	//if (prof) cudaProfilerStop();
		//prof = false;
		//double end = omp_get_wtime();
		//cout << 1000 * (end - start) << endl;
		Polynom<double> newpol = Polynom<double>(p.order, p.eps, p.vars, final_exponents, final_coeffs, nb_terms);
		result_pols[itmap->first] = newpol;
        }
	cudaFree(exp);
	cudaFree(coeff);
	cudaFree(final_exponents);
	cudaFree(final_coeffs);
	return Polmap<double>(result_pols);
} 
