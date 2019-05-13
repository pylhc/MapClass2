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
static std::vector<string> vars = {xstring, pxstring, ystring, pystring, dstring, sstring};

template <class T> Polynom<T> X (int order) {
	return Polynom<T>(order, 1E-18, vars, xstring, 1);	
}
template <class T> Polynom<T> PX (int order) {
	return Polynom<T>(order, 1E-18, vars, pxstring, 1);	
}
template <class T> Polynom<T> Y (int order) {
	return Polynom<T>(order, 1E-18, vars, ystring, 1);	
}
template <class T> Polynom<T> PY (int order) {
	return Polynom<T>(order, 1E-18, vars, pystring, 1);	
}
template <class T> Polynom<T> D (int order) {
	return Polynom<T>(order, 1E-18, vars, dstring, 1);	
}
template <class T> Polynom<T> S (int order) {
	return Polynom<T>(order, 1E-18, vars, sstring, 1);	
}

template <class T> Polynom<T> L(double length, Polynom<T> d) {
	return (d + 1).pinv() * length;
}

template <class T> Polmap<T> DRIFTMap(double length, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
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

template <class T> Polynom<T>  Q11(double L, double K1L, Polynom<T> d) {
	Polynom<T> K = (d + 1).pinv() * abs(K1L/L); 
	return cosp(sqrtp(K) * L);
} 

template <class T> Polynom<T> Q12(double L,double K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return (sqrtp(K).pinv() * sinp(sqrtp(K) * L))/(d + 1);
}

template <class T> Polynom<T> Q21(double L,double K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K) * sinp(sqrtp(K) * L) * (d + 1) * (-1);
}

template <class T> Polynom<T> Q33(double L,double K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return coshp(sqrtp(K) * L);
}

template <class T> Polynom<T> Q34(double L,double K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K).pinv() * sinhp(sqrtp(K) * L)/(d + 1);
}

template <class T> Polynom<T> Q43(double L,double K1L, Polynom<T> d) {
	Polynom<T> K =  (d + 1).pinv() * abs(K1L/L);
 	return sqrtp(K)*sinhp(sqrtp(K) * L) * (d + 1);
 	
}

template <class T> Polmap<T> QFMap(double L, double K1L, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
    unordered_map<string, Polynom<T>> mp;
  	Polynom<T> q11 = Q11(L, K1L, d);
  	Polynom<T> q12 = Q12(L, K1L, d);
  	Polynom<T> q21 = Q21(L, K1L, d);
  	Polynom<T> q33 = Q33(L, K1L, d);
  	Polynom<T> q34 = Q34(L, K1L, d);
  	Polynom<T> q43 = Q43(L, K1L, d);
  	mp[xstring] = q11 * x + q12 * px;
  	mp[pxstring] = q21 * x  + q11 * px;
  	mp[ystring] = q33 * y + q34 * py;
  	mp[pystring] = q43 * y + q33 * py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polmap<T> QDMap(double L,double K1L, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
  	unordered_map<string, Polynom<T>> mp;
  	Polynom<T> q11 = Q11(L, K1L, d);
  	Polynom<T> q12 = Q12(L, K1L, d);
  	Polynom<T> q21 = Q21(L, K1L, d);
  	Polynom<T> q33 = Q33(L, K1L, d);
  	Polynom<T> q34 = Q34(L, K1L, d);
  	Polynom<T> q43 = Q43(L, K1L, d);
  	mp[xstring] = q33 * x + q34 * px;
  	mp[pxstring] = q43 * x  + q33 * px;
  	mp[ystring] = q11 * y + q12 * py;
  	mp[pystring] = q21 * y + q11 * py;
  	mp[dstring] = d;
  	mp[sstring] = s;
  	return Polmap<T>(mp);
}

template <class T> Polynom<T>  D11(double L, double ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	return cosp(THETA);
} 

template <class T> Polynom<T>  D12(double L, double ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	Polynom<T> P = sqrtp(d + 1).pinv() * (L/ANGLE);
	return P * sinp(THETA);
}

template <class T> Polynom<T>  D15(double L, double ANGLE, Polynom<T> d) {
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	double P = L/ANGLE;
	return (cosp(THETA) * (-1) + 1) * P;
}

template <class T> Polynom<T>  D21(double L, double ANGLE, Polynom<T> d) {
	
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE;  
	Polynom<T> P = sqrtp(d + 1).pinv() * (L/ANGLE);
	return P.pinv() * sinp(THETA) * (-1);
}

template <class T> Polynom<T>  D25(double L, double ANGLE, Polynom<T> d) {
	
	Polynom<T> THETA = sqrtp(d + 1).pinv() * ANGLE; 
	return sinp(THETA) * sqrtp(d + 1);
}

template <class T> Polynom<T>  D34(double L, Polynom<T> d) {
	
	return (d + 1).pinv() * L;
}

template <class T> Polmap<T> RotMap(double e, double L, double ANGLE, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
  unordered_map<string, Polynom<T>> mp;
  double phi = 0;// correction term. Handbook of accelerators sec2.2 Pole face rotation matrix
  
  double d21  = ANGLE/L * tan(e) ;
  double d43  = -1*ANGLE/L * tan ( e - phi ) ;

  mp[xstring]  =  x ;
  mp[pxstring] = d21 * x  +  px ;
  mp[ystring]  =  y ;
  mp[pystring] = d43 * y  +  py ;
  mp[dstring]  = d;
  mp[sstring]  = s;

  return Polmap<T>(mp);
}

template <class T> Polmap<T> DIMap(double L,double ANGLE, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
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
	for (int i = 2; i <= n; ++i) {
		lst.push_back((lst[i - 2] * init) * (1./(i)));
	}
	return lst;
}

static vector<vector<Polynom<double>>> separateComplexList(vector<Polynom<complex<double>>> lst) {
	vector<vector<Polynom<double>>> res;
	vector<Polynom<double>> real;
	vector<Polynom<double>> imag;
	for (vector<Polynom<complex<double>>>::iterator it = lst.begin(); it != lst.end(); ++it)  {
		vector<Polynom<double>> r = separateComplex(*it);
		real.push_back(r[0]);
		imag.push_back(r[1]);
	}
	res.push_back(real);
	res.push_back(imag);
	return res;
}

template <class T> Polynom<T> Bx(Polynom<T> r, double KnL) {  
  	return r * (-KnL);
}
template <class T> Polynom<T> By(Polynom<T> i, double KnL) { 
  	return i * KnL;
}

template <class T> Polmap<T> MUL(double K1L, double K2L, double K3L, double K4L, vector<vector<Polynom<T>>> v,
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
	return mp;
}

static double CiK[7] = {0.04880952380952381, 0.2571428571428571, 0.03214285714285714, 0.3238095238095238, 0.03214285714285714, 0.2571428571428571, 0.04880952380952381};

template <class T> Polmap<T> MULTHICK(double K1L,double K2L,double K3L,double K4L,double L, vector<vector<Polynom<T>>> v, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	double Li = L / 6;
  	Polmap<T> m = generateDefaultMap(x, px, y, py, d, s);
  	Polmap<T> dr = DRIFTMap(Li, x, px, y, py, d, s);
	
  	for (int i = 0; i < 6; ++i) {
    		m = MUL(K1L*CiK[i], K2L*CiK[i], K3L*CiK[i], K4L*CiK[i], v, x, px, y, py, d, s) * m;
    		m =  dr * m;
	}
  	m = MUL(K1L*CiK[6], K2L*CiK[6], K3L*CiK[6], K4L*CiK[6], v, x, px, y, py, d, s) * m;
  	return m;
} 


template <class T> Polmap<T> MULTHICK2(double K1L,double K2L,double K3L,double K4L,double L, vector<vector<Polynom<T>>> v, Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s) {
	double Li = L / 2;
  	Polmap<T> m = DRIFTMap(Li, x, px, y, py, d, s);
	Polmap<T> dr = m;
    	m = MUL(K1L, K2L, K3L, K4L, v, x, px, y, py, d, s) * m;
    	m =  dr * m;
  	return m;
} 

template <class T> Polmap<T> mapForElement(unordered_map<string, string> e, vector<vector<Polynom<T>>> v,
		Polynom<T> x, Polynom<T> px, Polynom<T> y, Polynom<T> py, Polynom<T> d, Polynom<T> s,
		int fmultipole) {
	string keyword = e["KEYWORD"];
	//cout << keyword << endl;
	double l = atof(e["L"].c_str());
	double k1l = atof(e["K1L"].c_str());
	double k2l = atof(e["K2L"].c_str());
	double k3l = atof(e["K3L"].c_str());
	double k4l = atof(e["K4L"].c_str());
	double angle = atof(e["ANGLE"].c_str());
	double e1 = atof(e["E1"].c_str());
	double e2 = atof(e["E2"].c_str());

	string drift       = "DRIFT";
	string rcollimator = "RCOLLIMATOR";
	string ecollimator = "ECOLLIMATOR";
	string monitor     = "MONITOR";
	string matrix      = "MATRIX";
	string quadrupole  = "QUADRUPOLE";
	string sbend       = "SBEND";
	string multipole   = "MULTIPOLE";
	string decapole    = "DECAPOLE";
	string sextupole   = "SEXTUPOLE";
	string octupole    = "OCTUPOLE";
	
    	
    	if (keyword.compare(drift) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
	if (keyword.compare(rcollimator) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
    	if (keyword.compare(ecollimator) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
    	if (keyword.compare(monitor) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
    	if (keyword.compare(matrix) == 0) 
		return DRIFTMap(l, x, px, y, py, d, s);
    	if (keyword.compare(quadrupole) == 0)
		if (l != 0) {
		        if (k1l == 0)
			        return DRIFTMap(l, x, px, y, py, d, s);
        		if (k1l > 0)
          			return QFMap(l, k1l, x, px, y, py, d, s);
        		else
          			return QDMap(l, k1l, x, px, y, py, d, s);
          } 
	if (keyword.compare(sbend) == 0) {
	  Polmap<T> tmpmap;
	  if (angle == 0 ) tmpmap = DRIFTMap(l, x, px, y, py, d, s);
	  else tmpmap = DIMap(l, angle, x, px, y, py, d, s);
	  if ( e1 == 0 && e2 == 0 )
	    return tmpmap;
	  else
	    if ( e1 != 0 && e2 == 0)
	      return tmpmap * RotMap(e1, l, angle, x, px, y, py, d, s);
	  if ( e1 == 0 && e2 != 0)
	    return RotMap(e2, l, angle, x, px, y, py, d, s) * tmpmap ;
	  else
	    return RotMap(e2, l, angle, x, px, y, py, d, s) * tmpmap * RotMap(e1, l, angle, x, px, y, py, d, s);
	}  
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
			if (fmultipole == 0)
        			return MULTHICK(k1l,  k2l, k3l, k4l, l, v, x, px, y, py, d, s);
        		else return MULTHICK2(k1l,  k2l, k3l, k4l, l, v, x, px, y, py, d, s);
        	}
	}
	return Polmap<T>();
	//return generateDefaultMap( x, px, y, py, d, s);
}


