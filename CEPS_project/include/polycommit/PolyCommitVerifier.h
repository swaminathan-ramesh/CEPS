#ifndef __POLYCOMMITVERIFIER_H__
#define __POLYCOMMITVERIFIER_H__

#include "PolyCommitCommon.h"

bool verifypoly(const PolyCommitParams &p, const G1 &C, const ZZ_pX &f);
bool verifyvec(const PolyCommitParams &p, const G1 &C, const vec_ZZ_p &v);
bool verifyeval(const PolyCommitParams &p, const G1 &C, const ZZ_p &i,
    const ZZ_p &fi, const G1 &witness);
// verify a batch of evaluations (different polys, same evaluation point)
bool verifyeval(const PolyCommitParams &p, const vector<G1> &C, const ZZ_p &i,
    const vec_ZZ_p &fi, const vector<G1> &witness);

#endif
