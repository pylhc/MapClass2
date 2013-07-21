#include <iostream>
#include <algorithm>
#include <unordered_map>
#include <map>
#include <vector>
#include <string>
#include <cmath>
#include "boost/functional/hash.hpp"
#include "boost/lexical_cast.hpp"
#include <boost/algorithm/string.hpp>
#include <complex>
#include <sstream>
#include <omp.h>
#include "Pol.h"

using namespace std;
using namespace boost;

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
		Polmap parallel_composition(Polmap other);
		Polmap operator ^(int n);
		Polynom<T> operator [](string s);
		Polmap& truncate();
		Polmap& truncate(int order);
		int getOrder();
		double getEps();
		void printpolmap();
		string getMap();
		Polmap compositionWithoutHash(Polmap other);
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
			for (unsigned int k = 1; k < p.vars.size(); k ++) {
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
	//cout << "START" << endl; 
	int nrfound = 0;
	int nrmult = 0;		
	for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
		Polynom<T> p = itmap->second;
		Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
		vector<int> zero(p.vars.size(), 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
			m.clear();
			Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, i->second);
			for (unsigned int k = 0; k < p.vars.size(); k ++) {
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
					if (found != temp.end()) {
						t = found->second;
						nrfound ++;
					}
					else {
						int ilm = lexical_cast<int>(lm) + 1;
						lm = lexical_cast<string>(ilm);
						temp[pp] = lm;
						t = lm;
						Polynom<T> newv = tempvalues[get<0>(pp)] * tempvalues[get<1>(pp)]; 
						tempvalues[t] = newv;
						nrmult ++;
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
	//cout << nrfound << " " << nrmult << endl;
	return Polmap<T>(result_pols);
}

/*template <class T> Polmap<T> Polmap<T>::operator*(Polmap<T> other) {
        unordered_map<string, Polynom<T>> result_pols;
        unordered_map<string, Polynom<T>> temp;
        temp = other.pols;      
        list<string> m;
	string mstring = ""; 
        for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) {
                Polynom<T> p = itmap->second;                           
                Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
                for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
                        m.clear();
			mstring = "";
                        Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, i->second);
                        for (unsigned int k = 0; k < p.vars.size(); k ++) {
                                typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
                                if (found == other.pols.end())
                                        coef = coef * Polynom<T> (p.order, p.eps, p.vars[k], i->first[k]);  
                                else
                                        for (int a = 0; a < i->first[k]; a ++) {
                                                m.push_back(p.vars[k]);
						mstring += p.vars[k] + " ";
					}
                        }
                        if (m.size() > 0) {
                       		trim(mstring);
				typename unordered_map<string, Polynom<T>>::const_iterator found = temp.find (mstring);         
				if (found != temp.end())
					newpol = newpol + coef * found->second;
				else {
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
                                newpol = newpol + coef * temp[m.front()];
                		}
                        }
                        else 
                                newpol = newpol + coef;
                         
                } 
                result_pols[itmap->first] = newpol;
        }
        return Polmap<T>(result_pols);
}*/

template <class T> Polmap<T> Polmap<T>::parallel_composition(Polmap<T> other) {
	unordered_map<string, Polynom<T>> result_pols;
	unordered_map<pair<string, string>, string, pair_hash<string>> temp;
	unordered_map<string, Polynom<T>> tempvalues;
	vector<pair<string, Polynom<T>>> vpols;
	string lm = "0";
	tempvalues = other.pols;
	result_pols = pols; 
	for (typename unordered_map<string, Polynom<T>>::const_iterator itmap = pols.begin(); itmap != pols.end(); ++itmap) { 
		if (itmap->first.compare("s") != 0 && itmap->first.compare("d") != 0)
			vpols.push_back(make_pair(itmap->first, itmap->second));
	}
	list<string> m;
	string t;
	int nb_threads = vpols.size();
	#pragma omp parallel for firstprivate(lm, tempvalues, temp, m, t) shared(result_pols) schedule(static) num_threads(nb_threads)
	for (int itmap = 0; itmap < vpols.size(); itmap ++) {
		Polynom<T> p = vpols[itmap].second;
		Polynom<T> newpol = Polynom<T>(p.order, p.eps, p.vars, 0);
		for (typename unordered_map<vector<int>, T, container_hash<vector<int>>>::iterator i = p.terms.begin(); i != p.terms.end(); ++i) {
			m.clear();
			Polynom<T> coef = Polynom<T>(p.order, p.eps, p.vars, i->second);
			for (unsigned int k = 0; k < p.vars.size(); k ++) {
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
		result_pols[vpols[itmap].first] = newpol;
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
			for (unsigned int k = 0; k < p.vars.size(); ++k) {
				typename unordered_map<string, Polynom<T>>::const_iterator found = other.pols.find (p.vars[k]);
				if (found == other.pols.end()) 
					result = result * Polynom<T> (p.order, p.eps, p.vars[k], i->first[k]);
				else	for (int a = 0; a < i->first[k]; a ++)
						m.push_back(p.vars[k]);
			}
			for (vector<string>::iterator k = m.begin(); k != m.end(); k ++)
				result = result * other.pols[*k];
			newpol = newpol + result;
		}
		if (fabs(newpol.terms[zero]) < p.eps)
	          	newpol.terms.erase(zero);
		result_pols[itmap->first] = newpol;
	}
	return Polmap<T>(result_pols);
}                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             
