#ifndef MONT_CURVE_H_
#define MONT_CURVE_H_

#include "fp.h"
#include <assert.h>

typedef struct {
  fp X, Z;
} MP;

void MP_init(MP *Res);

void MP_set(MP *Res, const fp x, const fp z);

void MP_set_mpz(MP *Res, const mpz_t x, const mpz_t z);

void MP_set_si(MP *Res, const long int x, const long int z);

void MP_set_x(MP *Res, const fp x);

void MP_copy(MP *Res, const MP P);

void MP_clear(MP *Res);

void MP_norm(MP *Res, const MP P);

void xADD(MP *Res, const MP P, const MP Q, const MP R);

void xDBL(MP *Res, const MP P, const MP Ap24_C24);

void xDBLADD(MP *Dub, MP *Add, const MP P, const MP Q, const MP R,
             const MP Ap24_1);

void xMUL(MP *Res, const MP P, mpz_t k, const MP A_1);

void mont_rhs(fp *Res, const MP A, const fp x);

void xISOG(MP *Aout, MP *Pout, const MP A, const MP P, const MP K, int k);

#endif
