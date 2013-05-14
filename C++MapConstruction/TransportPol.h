#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <vector>
#include <string>
#include <cmath>
#include <complex>
#include "Polmap.h"
#include "Funset.h"

using namespace std;
using namespace boost;

static string xstring = "x";
static string pxstring = "px";
static string ystring = "y";
static string pystring = "py";
static string dstring = "d";
static string sstring = "s";


template <class T> Polynom<T> X (int order) {
	return Polynom<T>(order, 1E-18, xstring, 1);	
}
template <class T> Polynom<T> PX (int order) {
	return Polynom<T>(order, 1E-18, pxstring, 1);	
}
template <class T> Polynom<T> Y (int order) {
	return Polynom<T>(order, 1E-18, ystring, 1);	
}
template <class T> Polynom<T> PY (int order) {
	return Polynom<T>(order, 1E-18, pystring, 1);	
}
template <class T> Polynom<T> D (int order) {
	return Polynom<T>(order, 1E-18, dstring, 1);	
}
template <class T> Polynom<T> S (int order) {
	return Polynom<T>(order, 1E-18, sstring, 1);	
}

template <class T> Polynom<T> L(Polynom<T> length, Polynom<T> d) {
	return (d + 1).pinv() * length;
}

template <class T> Polmap<T> DRIFTMap(Polynom<T> length, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	unordered_map<string, Polynom<T>> mp;
  	Polynom<T> l = L(length, d);
     	 
  	mp[xstring] = l * px + x;
  	mp[pxstring] = px;
  	mp[ystring] = l * py + y;
	mp[pystring] = py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polynom<T>  Q11(Polynom<T> L, Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K = (d + 1).pinv() * abs(K1L/L); 
	return cosp(sqrtp(K) * L);
} 

template <class T> Polynom<T> Q12(Polynom<T> L,Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return (sqrtp(K).pinv() * sinp(sqrtp(K) * L))/(d + 1);
}

template <class T> Polynom<T> Q21(Polynom<T> L,Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K) * sinp(sqrtp(K) * L) * (d + 1) * (-1);
}

template <class T> Polynom<T> Q33(Polynom<T> L,Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return coshp(sqrtp(K) * L);
}

template <class T> Polynom<T> Q34(Polynom<T> L,Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K).pinv() * sinhp(sqrtp(K) * L)/(d + 1);
}

template <class T> Polynom<T> Q43(Polynom<T> L,Polynom<T> K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K)*sinhp(sqrtp(K) * L) * (d + 1);
 	
}

template <class T> Polmap<T> QFMap(Polynom<T> L, Polynom<T> K1L, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
    unordered_map<string, Polynom<T>> mp;
  	Polynom<T> q11 = Q11p(L, K1L, d);
  	Polynom<T> q12 = Q12p(L, K1L, d);
  	Polynom<T> q21 = Q21p(L, K1L, d);
  	Polynom<T> q33 = Q33p(L, K1L, d);
  	Polynom<T> q34 = Q34p(L, K1L, d);
  	Polynom<T> q43 = Q43p(L, K1L, d);
  	mp[xstring] = q11 * x + q12 * px;
  	mp[pxstring] = q21 * x  + q11 * px;
  	mp[ystring] = q33 * y + q34 * py;
  	mp[pystring] = q43 * y + q33 * py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polmap<T> QDMap(Polynom<T> L,Polynom<T> K1L, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
  	unordered_map<string, Polynom<T>> mp;
  	Polynom<T> q11 = Q11p(L, K1L, d);
  	Polynom<T> q12 = Q12p(L, K1L, d);
  	Polynom<T> q21 = Q21p(L, K1L, d);
  	Polynom<T> q33 = Q33p(L, K1L, d);
  	Polynom<T> q34 = Q34p(L, K1L, d);
  	Polynom<T> q43 = Q43p(L, K1L, d);
  	mp[xstring] = q33 * x + q34 * px;
  	mp[pxstring] = q43 * x  + q33 * px;
  	mp[ystring] = q11 * y + q12 * py;
  	mp[pystring] = q21 * y + q11 * py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polynom<T>  D11(Polynom<T> L, Polynom<T> ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	return cosp(THETA);
} 

template <class T> Polynom<T>  D12(Polynom<T> L, Polynom<T> ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	Polynom<T> P = sqrtp(d + 1).pinv() * (L/ANGLE);
	return P * sinp(THETA);
}

template <class T> Polynom<T>  D15(Polynom<T> L, Polynom<T> ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	Polynom<T> P = L/ANGLE;
	return (cosp(THETA) * (-1) + 1) * P;
}

template <class T> Polynom<T>  D21(Polynom<T> L, Polynom<T> ANGLE, Polynom<T> d) {
	
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE;  
	Polynom<T> P = sqrtp(d + 1).pinv() * (L/ANGLE);
	return P.pinv() * sinp(THETA) * (-1);
}

template <class T> Polynom<T>  D25(Polynom<T> L, Polynom<T> ANGLE, Polynom<T> d) {
	
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	return sinp(THETA) * sqrtp(d + 1);
}

template <class T> Polynom<T>  D34(Polynom<T> L, Polynom<T> d) {
	
	return (d + 1).pinv() * L;
}

template <class T> Polmap<T> DIMap(Polynom<T> L,Polynom<T> ANGLE, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	unordered_map<string, Polynom<T>> mp;
  	Polynom<T> d11 = D11(L, ANGLE, d);
  	Polynom<T> d12 = D12(L, ANGLE, d);
  	Polynom<T> d21 = D21(L, ANGLE, d);
  	Polynom<T> d15 = D15(L, ANGLE, d);
  	Polynom<T> d25 = D25(L, ANGLE, d);
  	Polynom<T> d34 = D34(L, d);
  	mp[xstring] = d11 * x  + d12  * px + d15 * d;
  	mp[pxstring] = d21 * x  + d11  * px + d25 * d;
  	mp[ystring] = y + d34 * py;
  	mp[pystring] = py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polmap<T> generateDefaultMap(Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	unordered_map<string, Polynom<T>> mp;
	mp[xstring] = x;
  	mp[pxstring] = px;
  	mp[ystring] = y;
  	mp[pystring] = py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

static vector<Polynom<double>> separateComplex(Polynom<complex<double>> p) {
	unordered_map<vector<int>, double, container_hash<vector<int>>> real_terms;
	unordered_map<vector<int>, double, container_hash<vector<int>>> imag_terms;
	for (unordered_map<vector<int>, complex<double>, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
		if (i->second.real() != 0)
			real_terms[i->first] = i->second.real();
		if (i->second.imag() != 0)
			imag_terms[i->first] = i->second.imag();	
	}
	vector<Polynom<double>> v;
	v.push_back(Polynom<double>(p.order, p.eps, p.vars, real_terms));
	v.push_back(Polynom<double>(p.order, p.eps, p.vars, imag_terms));
	return v;
}

static int factorial(int n){
  return (n == 1 || n == 0) ? 1 : factorial(n - 1) * n;
}

static vector<Polynom<complex<double>>> EQ(int n, int order) {
	Polynom<complex<double>> x = X<complex<double>>(order);
	Polynom<complex<double>> y = Y<complex<double>>(order);
	vector<Polynom<complex<double>>> lst;
	complex<double> j (0, 1);
	Polynom<complex<double>> init = x + y * j;
	
	lst.push_back(init);
	for (int i = 2; i <= n; i ++) {
		lst.push_back((lst[i - 2] * init) * (1./factorial(i)));
	}
	return lst;
}

static vector<vector<Polynom<double>>> separateComplexList(vector<Polynom<complex<double>>> lst) {
	vector<vector<Polynom<double>>> res;
	vector<Polynom<double>> real;
	vector<Polynom<double>> imag;
	for (vector<Polynom<complex<double>>>::iterator it = lst.begin(); it != lst.end(); it ++)  {
		vector<Polynom<double>> r = separateComplex(*it);
		real.push_back(r[0]);
		imag.push_back(r[1]);
	}
	res.push_back(real);
	res.push_back(imag);
	return res;
}
template <class T> Polynom<T> Bx(Polynom<T> r, Polynom<T> KnL) {  
  	return r * KnL * (-1);
}
template <class T> Polynom<T> By(Polynom<T> i, Polynom<T> KnL) { 
  	return i * KnL;
}

template <class T> Polmap<T> MUL(Polynom<T> K1L, Polynom<T> K2L, Polynom<T> K3L, Polynom<T> K4L, vector<vector<Polynom<T>>> v,
	Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	unordered_map<string, Polynom<T>> mp; 
	mp[xstring] = x;
  	mp[pxstring] = px;
  	mp[ystring] = y;
  	mp[pystring] = py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	 
  	if (abs(K1L) >= 1e-18) { 
    		
    		mp[pxstring] = mp[pxstring] + Bx(v[0][0], K1L);
    		mp[pystring] = mp[pystring] + By(v[1][0], K1L);	
  	}
  	if (abs(K2L) >= 1e-18) { 
    		
    		mp[pxstring] = mp[pxstring] + Bx(v[0][1], K2L);
    		mp[pystring] = mp[pystring] + By(v[1][1], K2L);	
  	}
  	if (abs(K3L) >= 1e-18) { 
    		
    		mp[pxstring] = mp[pxstring] + Bx(v[0][2], K3L);
    		mp[pystring] = mp[pystring] + By(v[1][2], K3L);	
  	}
  	if (abs(K4L) >= 1e-18) { 
    		
    		mp[pxstring] = mp[pxstring] + Bx(v[0][3], K4L);
    		mp[pystring] = mp[pystring] + By(v[1][3], K4L);	
  	}
  	
  	//m[pxstring] = m[pxstring] + Bx(1,K1L,order) + Bx(2,K2L,order) + Bx(3,K3L,order) + Bx(4,K4L,order)
  	//m[pystring] = m[pystring] + By(1,K1L,order) + By(2,K2L,order) + By(3,K3L,order) + By(4,K4L,order)
	return mp;
}

static double CiK[7] = {0.04880952380952381, 0.2571428571428571, 0.03214285714285714, 0.3238095238095238, 0.03214285714285714, 0.2571428571428571, 0.04880952380952381};

template <class T> Polmap<T> MULTHICK(Polynom<T> K1L,Polynom<T> K2L,Polynom<T> K3L,Polynom<T> K4L,Polynom<T> L, vector<vector<Polynom<T>>> v, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	Polynom<T> Li = L * (1/6);
  	Polmap<T> m = generateDefaultMap(x, px, y, py, d, s);
  	Polmap<T> dr = DRIFTMap(Li, x, px, y, py, d, s);
	
  	for (int i = 0; i < 6; i ++) {
    		m = MUL(K1L*CiK[i], K2L*CiK[i], K3L*CiK[i], K4L*CiK[i], v, x, px, y, py, d, s) * m;
    		m =  dr * m;
	}
  	m = MUL(K1L*CiK[6], K2L*CiK[6], K3L*CiK[6], K4L*CiK[6], v, x, px, y, py, d, s) * m;
  	return m;
} 

template <class T> Polmap<T> mapForElementWithIndex(int index, unordered_map<string, string> e, vector<vector<Polynom<T>>> v,
		Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	string keyword = e["KEYWORD"];
//	cout << keyword << endl;
	int order = x.order;
	double eps = x.eps;
	Polynom<T> l = Polynom<T>(order, eps, "l", atof(e["L"].c_str()));
	Polynom<T> k1l = Polynom<T>(order, eps, "k1l", atof(e["K1L"].c_str()));
	Polynom<T> k2l = Polynom<T>(order, eps, "k2l", atof(e["K2L"].c_str()));
	Polynom<T> k3l = Polynom<T>(order, eps, "k3l", atof(e["K3L"].c_str()));
	Polynom<T> k4l = Polynom<T>(order, eps, "k4l", atof(e["K4L"].c_str()));
	Polynom<T> angle = Polynom<T>(order, eps, "angle", atof(e["ANGLE"].c_str()));
	Polynom<T> delta = Polynom<T>(order, eps, "deltap", atof(e["DELTAP"].c_str()));
	string drift = "DRIFT";
	string quadrupole = "QUADRUPOLE";
	string sbend = "SBEND";
	string multipole = "MULTIPOLE";
	string decapole = "DECAPOLE";
	string sextupole = "SEXTUPOLE";
	string octupole = "OCTUPOLE";
    	
    	if (keyword.compare(drift) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
    	
    	if (keyword.compare(quadrupole) == 0)
		if (l != 0) { 
        		if (k1l > 0)
          		return QFMap(l, k1l, x, px, y, py, d, s); //or delta instead of d
        		else
          		return QDMap(l, k1l, x, px, y, py, d, s); //or delta instead of d
          } 
          
    	if (keyword.compare(sbend) == 0) 
     	return DIMap(l, angle, x, px, y, py, d, s);        
    	if (keyword.compare(quadrupole) == 0)
     	if (l == 0)
        		return MUL(k1l, k2l, k3l, k4l, v, x, px, y, py, d, s);
    	if (keyword.compare(multipole) == 0)
     	if (l == 0)
        		return MUL(k1l, k2l, k3l, k4l, v, x, px, y, py, d, s);
    	if (keyword.compare(sextupole) == 0 || keyword.compare(octupole) == 0 || keyword.compare(decapole) == 0) {
    		
     	if (l == 0) {

        		return MUL(k1l,  k2l, k3l, k4l, v, x, px, y, py, d, s);
        	}
      	else {
        		return MULTHICK(k1l,  k2l, k3l, k4l, l, v, x, px, y, py, d, s);
        	}
	}
	return Polmap<T>();
	//return generateDefaultMap( x, px, y, py, d, s);
}


