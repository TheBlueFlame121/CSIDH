#ifndef FP_H_
#define FP_H_

#include <gmp.h>
#include <stdlib.h>
#include <time.h>

typedef struct {
  mpz_t a;
} fp;

struct fp_params;
typedef struct fp_params fp_params;

struct fp_params {
  mpz_t p, tmp1, tmp2, tmp3;
  gmp_randstate_t state;
  fp fp_temp1, t0, t1, t2, t3, t4, t5, t6, fp_res;
};

extern fp_params global_params;

void fp_init(fp *Res);

void fp_setup(mpz_t prime);

void fp_clear(fp *Res);

void fp_set_str(fp *Res, const char *a);

void fp_set_mpz(fp *Res, const mpz_t a);

void fp_set_si(fp *Res, const long int a);

void fp_to_string(const fp x, char *a);

int fp_iszero(const fp x);

void fp_copy(fp *Res, const fp x);

int fp_cmp(const fp x, const fp y);

int fp_isone(const fp x);

void fp_random(fp *Res);

int fp_isequal(const fp x, const fp y);

void fp_add3(fp *Res, const fp x, const fp y);

void fp_add2(fp *Res, const fp x);

void fp_add3_si(fp *Res, const fp x, const long int a);

void fp_add2_si(fp *Res, const long int a);

void fp_sub3(fp *Res, const fp x, const fp y);

void fp_sub2(fp *Res, const fp x);

void fp_sub_si(fp *Res, const fp x, const long int a);

void fp_sub2_si(fp *Res, const long int a);

void fp_neg(fp *Res, const fp x);

void fp_mul3(fp *Res, const fp x, const fp y);

void fp_mul2(fp *Res, const fp x);

void fp_sqr2(fp *Res, const fp x);

void fp_sqr1(fp *Res);

void fp_inv2(fp *Res, const fp x);

void fp_inv1(fp *Res);

void fp_div3(fp *Res, const fp x, const fp y);

void fp_div2(fp *Res, const fp x);

#endif // FP_H_
