#ifndef __POLYCOMMITCOMMON_H__
#define __POLYCOMMITCOMMON_H__

#include <iostream>
#include <string>
#include <vector>
#include <NTL/ZZ_pX.h>
#include <NTL/ZZX.h>
#include <NTL/vec_ZZ.h>
#include <stdexcept>
#include <stdint.h>
#include "PBC.h"

#define POLYCOMMIT_MODE_POLY 0
#define POLYCOMMIT_MODE_VEC 1

#define POLYCOMMIT_WINDOW_DEFAULT 12

NTL_CLIENT

class PolyCommitParams {
friend istream& operator>>(istream &is, PolyCommitParams &params);
friend ostream& operator<<(ostream &os, PolyCommitParams &params);
friend class PolyCommitment;

    Pairing pairing;
    vector<G1> galphai;
    vector<G2> ghatalphai;
    vector<G1> glambdai;
    ZZ order;
    unsigned int t, t_hat;
    G1 h;					//Need a random h
    G1 halpha;				//h^alpha
    bool compute_Lagrange;
    int window_width;

    vector< vector<G1> > precomp_products;
    vector< vector<G1> > lambda_precomp_products;

public:
    PolyCommitParams(bool compute_Lagrange = false, int window_width = POLYCOMMIT_WINDOW_DEFAULT)
        : compute_Lagrange(compute_Lagrange), window_width(window_width) { }
    static void create(ostream &os, const string pairingparams,
			    unsigned int t);

    const Pairing &get_pairing() const { return pairing; }
    unsigned int get_t() const { return t; }
    const G1 &get_galphai(int i) const { return galphai[i]; }
    const G2 &get_ghatalphai(int i) const { return ghatalphai[i]; }
    const G1 &get_glambdai(int i) const { return glambdai[i]; }
    const G1 &get_h() const { return h; }
    const G1 &get_halpha() const { return halpha; }
    const ZZ &get_order() const { return order; }
};

class PolyCommitment {
friend class PoK_point_Prover;
friend class PoK_poly_Prover;

    const PolyCommitParams *paramsp;
    ZZ_pX f;
    vec_ZZ_p v;
    bool mode;
    G1 C_internal(const ZZX &expons, const vector< vector<G1> > &precomp_products) const;

public:
    // Commit to a polynomial of degree at most t.
    PolyCommitment(const PolyCommitParams *paramsp, const ZZ_pX &f) :
	paramsp(paramsp), f(f), mode(POLYCOMMIT_MODE_POLY) {
	if (deg(f) > (int)paramsp->get_t()) {
	    throw runtime_error("degree of f is too large");
	}}
    // Commit to a vector of length at most t+1.
    PolyCommitment(const PolyCommitParams *paramsp, const vec_ZZ_p &v) :
	paramsp(paramsp), v(v), mode(POLYCOMMIT_MODE_VEC) {
	int t_max = (int)paramsp->get_t() + 1;
	if (v.length() > t_max) {
	    throw runtime_error("length of v is too large");
	}
	f = poly_from_vec(v, t_max);
    }
    G1 get_C();
    G1 get_C_fast() const;
    G1 createwitness(const ZZ_p &i, const ZZ_p &blinding_factor = to_ZZ_p(1)) const;
    static vector<G1> createwitness(const vector<PolyCommitment> &commits, const ZZ_p &i);
    static ZZ_pX poly_from_vec(const vec_ZZ_p &v, int length);
};

Zr to_Zr(const Pairing &e, const ZZ_p x);
Zr to_Zr(const Pairing &e, const uint64_t x);

#endif
