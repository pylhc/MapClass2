#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include "boost/functional/hash.hpp"
#include <complex>
#include <numeric>
#include <iomanip> 

#ifndef POL_H
#define POL_H

using namespace std;

template <typename Container> 
struct container_hash {
    std::size_t operator()(Container const& c) const {
        return boost::hash_range(c.begin(), c.end());
    }
};


template  <class T> 
class Polynom {
	public:
		int order;
		double eps;
		vector<string> vars;
		unordered_map<vector<int>, T, container_hash<vector<int>>> terms;
	public:
		Polynom(int, double, vector<string>, unordered_map<vector<int>, T, container_hash<vector<int>>>);
		Polynom(int, double, string, int);
		Polynom(int, double, string);
		Polynom(int, double, T);
		Polynom(int, double, vector<string>, T);
		Polynom();
		Polynom(const Polynom &other);
		Polynom& operator= (Polynom other);
		Polynom(Polynom&& other);
		Polynom& truncate();
		Polynom& truncate(int order);
		Polynom operator +(Polynom other);
		Polynom operator +(T val);
		Polynom operator -(Polynom other);
		Polynom operator -(T val);
		Polynom operator *(Polynom other);
		Polynom operator *(T val);
		Polynom operator /(Polynom other);
		Polynom operator /(T val);
		Polynom operator ^(int n);
		Polynom eval(string var, Polynom val);
		Polynom eval(string var, T val);	
		friend void swap<>(Polynom<T>& first, Polynom<T>& second);
		T getZeroTerm();
		Polynom getRestOfTerms();
		void printpol();
		Polynom pinv();
		static Polynom phorner (vector<T> lst, Polynom p);
		void setOrder(int o);
		string getPol();
		template <class U> friend Polynom<U> operator*(U val, Polynom<U> other);
		template <class U> friend Polynom<U> operator+(U val, Polynom<U> other);
};

template <class T> Polynom<T>::Polynom (int o, double e, vector<string> v, unordered_map<vector<int>, T, container_hash<vector<int>>> t) {
	order = o;
	eps = e;
	vars = v;
	sort(vars.begin(), vars.end());
	terms = t;
}

template <class T> Polynom<T>::Polynom (int o, double e, T init) {
	order = o;
	eps = e;
	vector<int> v;
	terms[v] = init;
}


template <class T> Polynom<T>::Polynom (int o, double e, vector<string> v, T init) {
	order = o;
	eps = e;
	vars = v;
	vector<int> vv(v.size(), 0);
	terms[vv] = init;
}


template <class T> Polynom<T>::Polynom (int o, double e, string pol ) {
	order = o;
	eps = e;
	//start parsing of string pol
	//TODO:
}

template <class T> Polynom<T>::Polynom (int o, double e, string var, int exp) {
	order = o;
	eps = e;
	vars.push_back(var);
	vector<int> v;
	v.push_back(exp);
	terms[v] = 1;
}

template <class T> Polynom<T>::Polynom(const Polynom &other) {
	order = other.order;
     eps = other.eps;
     vars = other.vars;
     terms = other.terms;
}

template <class T> Polynom<T>::Polynom () {
	order = 6;
	eps = 1E-18;
}

template <class T> void swap(Polynom<T>& first, Polynom<T>& second) {
	swap(first.order, second.order);
	swap(first.eps, second.eps);
	swap(first.vars, second.vars);
	swap(first.terms, second.terms);
}

template <class T> Polynom<T>& Polynom<T>::operator=(Polynom other ) {
	swap(*this, other);
    	return *this;
}

template <class T> Polynom<T>:: Polynom(Polynom&& other) {
	order = other.order;
     eps = other.eps;
	vars = move(other.vars);
	terms = move(other.terms);
   
    	other.order = 0;
    	other.eps = 0; 
}

template <class T> Polynom<T>& Polynom<T>::truncate() {
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator it=terms.begin(); it != terms.end(); ++it)
		if (accumulate(it->first.begin(), it->first.end(), 0) > order || fabs(it->second) < eps)
			terms.erase (it);
	return *this;
}

template <class T> Polynom<T>& Polynom<T>::truncate(int o) {
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator it=terms.begin(); it != terms.end(); ++it)
		if (accumulate(it->first.begin(), it->first.end(), 0) > o || fabs(it->second) < eps)
			terms.erase (it);
	return *this;
}

template <class T> Polynom<T> Polynom<T>::operator*(Polynom<T> other) {
	int neworder = min(order, other.order);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	
	vector<string> newvars;
	if (vars == other.vars) {
		newvars = vars;			
		vector<int> newexp (newvars.size(), 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator i = terms.begin(); i != terms.end(); ++i) 
			for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator j = other.terms.begin(); j != other.terms.end(); ++j) {
				for (int k = 0; k < i->first.size(); k++)  
					newexp[k] = i->first[k] + j->first[k];
				if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder && fabs(i->second * j->second) >= eps)
					result_terms[newexp] += i->second * j->second;  
			}
	}
	else {
		newvars.resize (vars.size() + other.vars.size());
		vector<string> :: iterator it = set_union (vars.begin(), vars.end(), other.vars.begin(), other.vars.end(), newvars.begin()); 
		newvars.resize(it - newvars.begin());
		vector<int> newexp (newvars.size(), 0);
		unordered_map<string, int> expi, expj;
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator i = terms.begin(); i != terms.end(); ++i) {
			for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator j = other.terms.begin(); j != other.terms.end(); ++j) {
				
				for (int k = 0; k < vars.size(); k ++)
					expi[vars[k]] = i->first[k];
				for (int k = 0; k < other.vars.size(); k ++)
					expj[other.vars[k]] = j->first[k];
							
				for (int k = 0; k < newvars.size(); k++) {
					unordered_map<string,int>::iterator it1 = expi.find(newvars[k]);
					unordered_map<string,int>::iterator it2 = expj.find(newvars[k]);
					newexp[k] = 0;
					if (it1 != expi.end())
						newexp[k] += it1->second;
					if (it2 != expj.end())
						newexp[k] += it2->second;
				} 
				
				if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder && fabs(i->second * j->second) >= eps)
					result_terms[newexp] += i->second * j->second;	
			}
		}
			
	}
	Polynom<T> res = Polynom<T>(neworder, eps, newvars, result_terms);
	return move(res);
}

template <class T> Polynom<T> Polynom<T>::operator/(Polynom<T> other) {
	return move(*this * other.pinv());

}

template <class T> Polynom<T> Polynom<T>::operator+(Polynom<T> other) {
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	int neworder = min(order, other.order);
	vector<string> newvars (vars.size() + other.vars.size());
	vector<string> :: iterator it = set_union (vars.begin(), vars.end(), other.vars.begin(), other.vars.end(), newvars.begin());
	newvars.resize(it - newvars.begin());
	vector<int> newexp (newvars.size(), 0);
	unordered_map<string, int> expi, expj;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator i = terms.begin(); i != terms.end(); ++i) 
	{
		for (int k = 0; k < vars.size(); k ++)
			expi[vars[k]] = i->first[k];
		for (int k = 0; k < newvars.size(); k ++) {
			unordered_map<string, int>::iterator it1 = expi.find(newvars[k]);
			if (it1 != expi.end())
				newexp[k] = it1->second;
			else newexp[k] = 0;
		}
		if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder)
			result_terms[newexp] += i->second;
		
	}
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator j = other.terms.begin(); j != other.terms.end(); ++j) 
	{		
		for (int k = 0; k < other.vars.size(); k ++)
			expj[other.vars[k]] = j->first[k];
		for (int k = 0; k < newvars.size(); k ++) {
			unordered_map<string, int>::iterator it1 = expj.find(newvars[k]);
			if (it1 != expj.end())
				newexp[k] = it1->second;
			else newexp[k] = 0;
		}
		if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder)
			result_terms[newexp] += j->second;
	}
	Polynom<T> res = Polynom<T>(neworder, eps, newvars, result_terms);
	return move(res);
}

template <class T> Polynom<T> Polynom<T>::operator-(Polynom<T> other) {
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	int neworder = min(order, other.order);
	vector<string> newvars (vars.size() + other.vars.size());
	vector<string> :: iterator it = set_union (vars.begin(), vars.end(), other.vars.begin(), other.vars.end(), newvars.begin());
	newvars.resize(it - newvars.begin());
	vector<int> newexp(newvars.size(), 0);
	unordered_map<string, int> expi, expj;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator i = terms.begin(); i != terms.end(); ++i) 
	{
		for (int k = 0; k < vars.size(); k ++)
			expi[vars[k]] = i->first[k];
		for (int k = 0; k < newvars.size(); k ++) {
			unordered_map<string, int>::iterator it1 = expi.find(newvars[k]);
			if (it1 != expi.end())
				newexp[k] = it1->second;
			else newexp[k] = 0;
		}
		if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder)
			result_terms[newexp] += i->second;
		
	}
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator j = other.terms.begin(); j != other.terms.end(); ++j) 
	{		
		for (int k = 0; k < other.vars.size(); k ++)
			expj[other.vars[k]] = j->first[k];
		for (int k = 0; k < newvars.size(); k ++) {
			unordered_map<string, int>::iterator it1 = expj.find(newvars[k]);
			if (it1 != expj.end())
				newexp[k] = it1->second;
			else newexp[k] = 0;

		}
		if (accumulate(newexp.begin(), newexp.end(), 0) <= neworder)
			result_terms[newexp] -= j->second;
	}
	Polynom<T> res = Polynom<T>(neworder, eps, newvars, result_terms);
	return move(res);
}


template <class T> Polynom<T> Polynom<T>::operator^(int n) {
	
	if (n == 0) 
		return Polynom<T>(order, eps, vars, 1);
	else if (n < 0) 
		return pinv() ^ n;
	Polynom<T> res = *this;
	Polynom<T> p = *this;
	for (int i = 0; i < n - 1; i ++) 
		res = res * p;
	return move(res);

}

template <class T> void Polynom<T>::printpol() {
	bool first = true;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = terms.begin(); i != terms.end(); ++i) {
		if (i->second != 0) {
			double val = abs(i->second);
			if (i->second < 0)
				cout << " - ";
			else if(!first)
				cout <<" + ";
			cout << setprecision(14) << val;
			for (int k = 0; k < i->first.size(); k ++) {
				if (i->first[k] >= 2)
					cout << "*" << vars[k] << "**" << i->first[k];
				else if (i->first[k] == 1)
					cout << "*" << vars[k];
			if (first)
				first = false;	
			}
			
		}
	}
	cout << endl;
}

template <class T> string Polynom<T>::getPol() {
	bool first = true;
	ostringstream temp;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = terms.begin(); i != terms.end(); ++i) {
		if (i->second != 0) {
			double val = abs(i->second);
			if (i->second < 0)
				temp << " - ";
			else if(!first)
				temp <<" + ";
			temp << setprecision(14) << val;
			for (int k = 0; k < i->first.size(); k ++) {
				if (i->first[k] >= 2)
					temp << "*" << vars[k] << "**" << i->first[k];
				else if (i->first[k] == 1)
					temp << "*" << vars[k];
			if (first)
				first = false;	
			}
			
		}
	}
	return temp.str();
}

template <class T> Polynom<T> Polynom<T>:: eval(string var, Polynom<T> other) {
	vector<string>::iterator itv = find (vars.begin(), vars.end(), var);
	if(itv == vars.end())
		return *this; 
	int index = itv - vars.begin();
    	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	vector<string> valvar;
	valvar.push_back(var);
	//delete var from vars and add new vars from pol val
	vector<string> newvars (vars.size() + other.vars.size());
	vector<string> lessvars (vars.size() - 1);
	set_difference (vars.begin(), vars.end(), valvar.begin(), valvar.end(), lessvars.begin());
	vector<string> :: iterator it = set_union (lessvars.begin(), lessvars.end(), other.vars.begin(), other.vars.end(), newvars.begin());
	newvars.resize(it - newvars.begin());
	vector<int> newexp (newvars.size(), 0);
	unordered_map<string, int> expi;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = terms.begin(); i != terms.end(); ++i) {
		for (int k = 0; k < vars.size(); k ++)
			expi[vars[k]] = i->first[k];
		if (i->first[index] != 0) {		
			Polynom resexp = other ^ i->first[index];
			vector<int> current_term_exp = i->first;
			current_term_exp.erase(current_term_exp.begin() + index);
			unordered_map<vector<int>, T, container_hash<vector<int>>> current_term_map;
			current_term_map[current_term_exp] = i->second; 
			Polynom<T> current_term = Polynom<T>(order, eps, lessvars, current_term_map);
			current_term = current_term * resexp;
			result_terms.insert (current_term.terms.begin(), current_term.terms.end()); 
		}
		else {
			for (int k = 0; k < newvars.size(); k++) {
				unordered_map<string,int>::iterator it1 = expi.find(newvars[k]);	
				newexp[k] = 0;
				if (it1 != expi.end())
					newexp[k] += it1->second;	
			}
			if (accumulate(newexp.begin(), newexp.end(), 0) <= order)
				result_terms[newexp] += i->second;
		}
	
	}
	Polynom<T> res = Polynom<T>(order, eps, newvars, result_terms);
	return move(res);
}

template <class T> Polynom<T> Polynom<T>:: eval(string var, T val) {
	vector<string>::iterator it = find (vars.begin(), vars.end(), var);
	if(it == vars.end())
		return *this; 
	int index = it - vars.begin();
    
	vector<string> newvars (vars.size());
	newvars = vars;
	newvars.erase(newvars.begin() + index);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	vector<int> newexp;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = terms.begin(); i != terms.end(); ++i) {	
		newexp = i->first;
		newexp.erase(newexp.begin() + index);
		result_terms[newexp] += i->second * pow(val, i->first[index]);
	}
	Polynom<T> res = Polynom<T>(order, eps, newvars, result_terms);
	return move(res);

}

template <class T> T Polynom<T>::getZeroTerm () {
	vector<int> zero (vars.size(), 0);
	typename unordered_map<vector<int>, T, container_hash<vector<int>>>::const_iterator i = terms.find(zero);
	if (i == terms.end()) 
		return 0;
	else return i->second;
}


template <class T> Polynom<T> Polynom<T>::getRestOfTerms() {
	vector<int> zero;
	for (int i = 0; i < vars.size(); i ++)
		zero.push_back(0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = terms;
	result_terms.erase(zero);
	Polynom<T> res = Polynom<T>(order, eps, vars, result_terms);
	return move(res);
}


template <class T> Polynom<T> Polynom<T>::operator+(T val) {
	vector<int> zero (vars.size(), 0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = terms;
	result_terms[zero] += val;
	Polynom<T> res = Polynom<T>(order, eps, vars, result_terms);
	return move(res);
}

template <class T> Polynom<T> Polynom<T>::operator-(T val) {
	vector<int> zero (vars.size(), 0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = terms;
	result_terms[zero] -= val;
	Polynom<T> res = Polynom<T>(order, eps, vars, result_terms);
	return move(res);
}

template <class T> Polynom<T> Polynom<T>::operator*(T val) {
	vector<int> zero (vars.size());
	for (int i = 0; i < vars.size(); i ++)
		zero.push_back(0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = terms;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = terms.begin(); i != terms.end(); ++i) {
		result_terms[i->first] = i->second * val;
	}
	Polynom<T> res = Polynom<T>(order, eps, vars, result_terms);
	return move(res);
}

template <class T> Polynom<T> operator+(T val, Polynom<T> other) {
	vector<int> zero (other.vars.size(), 0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = other.terms;
	result_terms[zero] += other.val;
	Polynom<T> res = Polynom<T>(other.order, other.eps, other.vars, result_terms);
	return move(res);
}

template <class T> Polynom<T> operator*(T val, Polynom<T> other) {
	vector<int> zero (other.vars.size());
	for (int i = 0; i < other.vars.size(); i ++)
		zero.push_back(0);
	unordered_map<vector<int>, T, container_hash<vector<int>>> result_terms;
	result_terms = other.terms;
	for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = other.terms.begin(); i != other.terms.end(); ++i) {
		result_terms[i->first] = i->second * val;
	}
	Polynom<T> res = Polynom<T>(other.order, other.eps, other.vars, result_terms);
	return move(res);
}


template <class T> Polynom<T> Polynom<T>::operator/(T val) {
	return move(*this * (1/val));
}

template <class T> Polynom<T> Polynom<T>::pinv() {
	T a0 = getZeroTerm();
	Polynom<T> p = getRestOfTerms();
	p = p/a0;
	vector<T> lst;
	lst.push_back(1/a0);
	for (int i = 1; i < order + 1; i ++) 
		lst.push_back(-lst[i - 1]);
	return move(phorner(lst, p));
}

template <class T> Polynom<T> Polynom<T>::phorner (vector<T> lst, Polynom<T> p) {
	Polynom<T> out = Polynom(p.order, p.eps, lst[lst.size() - 1]);
	for (int i = lst.size() - 2; i >= 0; i --) {
		out = out * p;
		out = out + lst[i];
	}
	return move(out);	
}

template <class T> void Polynom<T>::setOrder(int o) {
	order = o;
}

#endif
