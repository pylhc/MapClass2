#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <boost/functional/hash.hpp>
#include <complex>
#include "Pol.h"

using namespace std;

template <class T> Polynom<T> sqrtp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	p = p/a0;
	vector<T> lst;
	lst.push_back(sqrt(a0));
	for (int i = 1; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 1]/2/i*(2*i - 3));
	return Polynom<T>::phorner(lst, p);
}

template <class T> Polynom<T> expp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	vector<T> lst;
	lst.push_back(exp(a0));
	for (int i = 1; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 1]/i);
	return Polynom<T>::phorner(lst, p);
}

template <class T> Polynom<T> logp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	p = p/a0;
	vector<T> lst;
	lst.push_back(log(a0));
	lst.push_back(1);
	for (int i = 2; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 1]/i * (i - 1));
	return Polynom<T>::phorner(lst, p);
}


template <class T> Polynom<T> sinp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	vector<T> lst;
	lst.push_back(sin(a0));
	lst.push_back(cos(a0));
	for (int i = 2; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 2]/i/(i - 1));
	return Polynom<T>::phorner(lst, p);
}

template <class T> Polynom<T> cosp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	vector<T> lst;
	lst.push_back(cos(a0));
	lst.push_back(-sin(a0));
	for (int i = 2; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 2]/i/(i - 1));
	return Polynom<T>::phorner(lst, p);
}

template <class T> Polynom<T> tanp(Polynom<T> c) {
	return sinp(c)/cosp(c);
}

template <class T> Polynom<T> sinhp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	vector<T> lst;
	lst.push_back(sinh(a0));
	lst.push_back(cosh(a0));
	for (int i = 2; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 2]/i/(i - 1));
	return Polynom<T>::phorner(lst, p);
}

template <class T> Polynom<T> coshp(Polynom<T> c) {
	T a0 = c.getZeroTerm();
	Polynom<T> p = c.getRestOfTerms();
	vector<T> lst;
	lst.push_back(cosh(a0));
	lst.push_back(sinh(a0));
	for (int i = 2; i < c.order + 1; i ++) 
		lst.push_back(-lst[i - 2]/i/(i - 1));
	return Polynom<T>::phorner(lst, p);
}
