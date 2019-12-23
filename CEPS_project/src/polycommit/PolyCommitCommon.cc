#include "PolyCommitCommon.h"

#include <NTL/ZZ_p.h>
#include <arpa/inet.h>

void PolyCommitParams::create(ostream &os, const string pairingparams,
	unsigned int t)
{
    // Create the Pairing
    Pairing p(pairingparams);
    size_t len = pairingparams.length();
    unsigned char lenbuf[2];
    lenbuf[0] = len >> 8;
    lenbuf[1] = len;
    os.write((const char *)lenbuf, 2);
    os << pairingparams;

    // Create a random generator for G1
    unsigned int t_nbo;  // network byte order
    t_nbo = htonl(t);
    os.write((const char *)&t_nbo, 4);
    G1 g(p, false);
    os << g;

    // Pick a random alpha
    Zr alpha(p, true);

    for (unsigned int i=0; i<t; ++i) {
	g ^= alpha;
	os << g;
    }

    // Create a random generator for G1 -- for h
    G1 h(p, false);
    os << h;

    //  -- for h^alpha

    h ^= alpha;
    os << h;

    // This controls how many entries can be in a single batch opening
    // of values of a single polynomial
    unsigned int t_hat = 1;
    unsigned int t_hat_nbo;  // network byte order

    // Create a random generator for G2
    t_hat_nbo = htonl(t_hat);
    os.write((const char *)&t_hat_nbo, 4);
    G2 ghat(p, false);
    os << ghat;

    for (unsigned int i=0; i<t_hat; ++i) {
	ghat ^= alpha;
	os << ghat;
    }
}

Zr to_Zr(const Pairing &e, const ZZ_p x)
{
    long len = NumBytes(rep(x));
    unsigned char xbytes[len];
    BytesFromZZ(xbytes, rep(x), len);
    unsigned char revbytes[len];
    for (long i=0;i<len;++i) {
	revbytes[i] = xbytes[len-1-i];
    }
    Zr zr(e, revbytes, len, 0);
    return zr;
}

Zr to_Zr(const Pairing &e, const uint64_t x)
{
    long len = sizeof(uint64_t);
    unsigned char revbytes[len];
    for (long i=0;i<len;++i) {
	revbytes[i] = (char)(x >> (8 * (len-1-i)));
    }
    Zr zr(e, revbytes, len, 0);
    return zr;
}


ostream &operator<<(ostream &os, PolyCommitParams &params)
{
    // Write the pairing params
    size_t len = params.pairing.get_pbc_param_t().length();
    unsigned char lenbuf[2];
    lenbuf[0] = len >> 8;
    lenbuf[1] = len;
    os.write((const char *)lenbuf, 2);
    os << params.pairing.get_pbc_param_t();

    // Write the g^\alpha^i
    unsigned int t_nbo;  // network byte order
    t_nbo = htonl(params.t);
    os.write((const char *)&t_nbo, 4);
    for (unsigned int i=0; i<=params.t; ++i) {
	os << params.galphai[i];
    }

    // Write h,h^alpha
    os << params.h << params.halpha;

    // Write the ghat^\alpha^i
    unsigned int t_hat_nbo;  // network byte order
    t_hat_nbo = htonl(params.t_hat);
    os.write((const char *)&t_hat_nbo, 4);
    for (unsigned int i=0; i<=params.t_hat; ++i) {
	os << params.ghatalphai[i];
    }

 return os;
}

istream &operator>>(istream &is, PolyCommitParams &params)
{
    unsigned char lenbuf[2];
    is.read((char *)lenbuf, 2);
    size_t len = (lenbuf[0] << 8) + lenbuf[1];
    char pairingparams[len+1];
    is.read(pairingparams, len);
    pairingparams[len] = '\0';
    params.pairing.init(pairingparams);

    // Read the maximum exponent of g^\alpha^i
    unsigned int t;
    unsigned int t_nbo;  // network byte order
    is.read((char *)&t_nbo, 4);
    t = ntohl(t_nbo);

    // Read the g^\alpha^i
    for (unsigned int i=0; i<=t; ++i) {
	G1 g(params.pairing, true);
	is >> g;
	params.galphai.push_back(g);
    }
    params.t = t;

    // Create a random generator for G1 -- for h
    G1 h(params.pairing, true);
    is >> h;
    params.h = h;
    //  -- for h^alpha
    G1 halpha(params.pairing, true);
    is >> halpha;
    params.halpha = halpha;

    // Read the maximum exponent of ghat^\alpha^i
    unsigned int t_hat;
    unsigned int t_hat_nbo;  // network byte order
    is.read((char *)&t_hat_nbo, 4);
    t_hat = ntohl(t_hat_nbo);

    // Read the ghat^\alpha^i
    for (unsigned int i=0; i<=t_hat; ++i) {
	G2 g(params.pairing, true);
	is >> g;
	params.ghatalphai.push_back(g);
    }
    params.t_hat = t_hat;

    // Figure out the order of the pairing groups.  There doesn't seem
    // to be an API call for this?
    Zr zero(params.pairing);
    Zr one(zero, 1);
    zero -= one;
    string negone = zero.toString();
    // Convert it to little-endian
    int l = negone.length();
    const unsigned char *negonestr = (const unsigned char *)negone.c_str();
    unsigned char orderm1[l+1];
    for (int i=0;i<l;++i) {
	orderm1[i] = negonestr[l-1-i];
    }
    orderm1[l] = '\0';
    // Convert it to a ZZ
    params.order = ZZFromBytes(orderm1, l);
    params.order += 1;

    //Save ZZ_p context
    ZZ_pContext savectx;
    savectx.save();

    ZZ_p::init(params.order);

    // Compute the precomputation table
    G1 gzero(params.pairing, true);
    int windows = t / params.window_width + 1;
    params.precomp_products.resize(windows);
    for (int k=0;k<windows;++k)
    {
    	params.precomp_products[k].push_back(gzero);
	unsigned int bound = (k==windows-1) ? (t+1) % params.window_width : params.window_width;
	for(unsigned int i=0;i<bound;++i) {
	for(int j=0;j<(1<<i);++j) {
	    params.precomp_products[k].push_back(
		params.precomp_products[k][j] *
		params.galphai[k*params.window_width+i]);
	}
    }
    }
    // Optionally, compute the Lagrange basis
    if (params.compute_Lagrange)
    {
	vec_ZZ_p indices;
	for (unsigned int i = 0; i <= t; ++i) { append(indices, to_ZZ_p(i)); }
	ZZ_pX lambda = BuildFromRoots(indices);

	for (unsigned int i = 0; i <= t; ++i)
	{
	    ZZ_pX lambdai = lambda / (ZZ_pX(1,1) - to_ZZ_p(i));
	    lambdai /= eval(lambdai, to_ZZ_p(i));
	    PolyCommitment C(&params, lambdai);
	    params.glambdai.push_back(C.get_C_fast());
	}

	// Compute the precomputation table
	params.lambda_precomp_products.resize(windows);
	for (int k=0;k<windows;k++)
	{
	    params.lambda_precomp_products[k].push_back(gzero);
	    unsigned int bound = (k==windows-1) ? (t+1) % params.window_width : params.window_width;
	    for(unsigned int i=0;i<bound;++i) {
		for(int j=0;j<(1<<i);++j) {
		    params.lambda_precomp_products[k].push_back(
			params.lambda_precomp_products[k][j] *
			params.glambdai[k*params.window_width+i]);
		}
	    }
	}
    }

    // Restore context
    savectx.restore();

    return is;
}

ZZ_pX PolyCommitment::poly_from_vec(const vec_ZZ_p &v, int length)
{
    vec_ZZ_p indices;
    for (int i = 0 ; i < length; ++i) { append(indices, to_ZZ_p(i)); }

    return interpolate(indices, VectorCopy(v, length));
}


G1 PolyCommitment::get_C() {
    //Save ZZ_p context
    ZZ_pContext savectx;
    savectx.save();

    ZZ_p::init(paramsp->order);
    G1 g(paramsp->pairing);

    if ((mode == POLYCOMMIT_MODE_VEC) && (paramsp->compute_Lagrange))
    {
	int len = v.length();
	for(int i=0;i<len;++i) {
	    ZZ_p cfp = v[i];
	    Zr cf = to_Zr(paramsp->pairing, cfp);
	    g *= ((paramsp->glambdai[i])^cf);
	}
    }
    else
    {
	int d = deg(f);
	for(int i=0;i<=d;++i) {
	    ZZ_p cfp = coeff(f,i);
	    Zr cf = to_Zr(paramsp->pairing, cfp);
	    g *= ((paramsp->galphai[i])^cf);
	}
    }

    // Restore context
    savectx.restore();

    return g;
}

static void divrem2(ZZX &q, ZZX &r, const ZZX &ex)
{
    int d = deg(ex);
    for (int i=0;i<=d;++i) {
	SetCoeff(q, i, coeff(ex, i)/2);
	SetCoeff(r, i, coeff(ex, i)%2);
    }
}

G1 PolyCommitment::C_internal(const ZZX &expons, const vector< vector<G1> > &precomp_products) const
{
    if (IsZero(expons)) {
	return G1(paramsp->pairing, true);
    }
    ZZX q, r;
    divrem2(q, r, expons);
    G1 C = C_internal(q, precomp_products);
    int d = deg(r);

    C = C.square();
    int windows = d / paramsp->window_width + 1;
    for (int k=0;k<windows;++k)
    {
	unsigned int index=0;
	for (int i=0;i<paramsp->window_width;++i) {
	    if (!IsZero(coeff(r,k*paramsp->window_width+i))) {
		index += (1<<i);
	    }
	}
	C *= precomp_products[k][index];
    }
    return C;
}

G1 PolyCommitment::get_C_fast() const {
    if ((mode == POLYCOMMIT_MODE_VEC) && (paramsp->compute_Lagrange))
    {
	int len = v.length();
	ZZX expons;
	for (int i=0;i<len;++i) {
	    SetCoeff(expons, i, rep(v[i]));
	}
	return C_internal(expons, paramsp->lambda_precomp_products);
    }
    else
    {
	int d = deg(f);
	ZZX expons;
	for(int i=0;i<=d;++i) {
	    SetCoeff(expons, i, rep(coeff(f, i)));
	}

	return C_internal(expons, paramsp->precomp_products);
    }
}

G1 PolyCommitment::createwitness(const ZZ_p &i, const ZZ_p &blinding_factor) const {
    ZZ_p fi = eval(f, i);

    // Compute the polynomial psi(x) = (f(x)-f(i))/(x-i)
    ZZ_pX denom(1,1); // x
    denom -= i;
    ZZ_pX psi;
    if (divide(psi, f-fi, denom) == 0) {
	throw runtime_error("x-i doesn't divide f(x)-f(i)?!");
    }

    // Compute w = g^{psi(alpha)}, which is just a commitment to psi
    PolyCommitment psicom(paramsp, blinding_factor * psi);

    return psicom.get_C_fast();
}

vector<G1> PolyCommitment::createwitness(const vector<PolyCommitment> &commits, const ZZ_p &i)
{
    vector<G1> witnesses;
    for (size_t j = 0; j < commits.size(); j++)
    {
	witnesses.push_back(commits[j].createwitness(i));
    }
    return witnesses;
}

#define POLYCOMMIT_SYM_DEFPARAMS \
"type a\n"\
"q 2149279669255358467807031928884602064965849607417268878942578020380711849468507704854756169974266156283454009766052915551928876758979436718043801993006923\n"\
"h 2941193476968633928915514572756480069153912946528476089236551260848359597655536307119448058353958150131532\n"\
"r 730750862221594424981965739670091261094297337857\n"\
"exp2 159\n"\
"exp1 135\n"\
"sign1 1\n"\
"sign0 1\n"

#define POLYCOMMIT_ASYM_DEFPARAMS \
"type d\n"\
"q 34155376728415643485253199952223120817480971620809317\n"\
"n 34155376728415643485253199767411387617406135340735963\n"\
"h 14079\n"\
"r 2425980306017163398341728799446792216592523285797\n"\
"a 16998822741009851476858453412953856547198847906791629\n"\
"b 32465094389155591141149001016286335587250188710863434\n"\
"k 6\n"\
"nk 1587648945903454419702294240696669428988257567004066667807817123624472417955040581724182593942606411788830338483111739923094122455965212711863453680770344432635484301107088581083365986465995771210066008136141674961980627026805959407325799760861946153887760541539723935858058246026804124025609796222220660565141631296\n"\
"hk 269761481129543709905880358277500023395314481104084120305961073572355019695601034021874289343986546858568479155991309118249358857985438172361312160865276114790621806863462094233703391155493997120107408194663831963226944\n"\
"coeff0 26826430320383533291086636892391482089885019715604831\n"\
"coeff1 10041240416257104004698194348998922958451952509010108\n"\
"coeff2 6606091442536002616272024969858190106364268669962931\n"\
"nqr 4719190100394795070479603407875910135476373878231286\n"


#ifdef CREATE_PARAMS
#include <iostream>

int main(int argc, char **argv)
{
    int t = argc > 1 ? atoi(argv[1]) : 5;

    PolyCommitParams::create(cout, string(POLYCOMMIT_SYM_DEFPARAMS), t);
    return 0;
}
#endif

#ifdef READ_PARAMS
#include <iostream>

int main(int argc, char **argv)
{
    PolyCommitParams params;

    cin >> params;
    return 0;
}
#endif
