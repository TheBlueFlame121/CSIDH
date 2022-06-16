#include "fp.h"

fp_params global_params;

void fp_init(fp *Res) { mpz_init(Res->a); }

void fp_setup(mpz_t prime) {
  mpz_init(global_params.p);
  mpz_init(global_params.tmp1);
  mpz_init(global_params.tmp2);
  mpz_init(global_params.tmp3);

  mpz_set(global_params.p, prime);

  // init random
  srand(time(0));
  gmp_randinit_mt(global_params.state);
  gmp_randseed_ui(global_params.state, time(0));

  // init some fp2 registers
  fp_init(&(global_params.fp_temp1));
  fp_init(&(global_params.t0));
  fp_init(&(global_params.t1));
  fp_init(&(global_params.t2));
  fp_init(&(global_params.t3));
  fp_init(&(global_params.t4));
  fp_init(&(global_params.t5));
  fp_init(&(global_params.t6));
  fp_init(&(global_params.fp_res));
}

void fp_clear(fp *Res) { mpz_clear(Res->a); }

void fp_set_str(fp *Res, const char *a) {
  mpz_set_str(global_params.tmp1, a, 0);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_set_mpz(fp *Res, const mpz_t a) { mpz_mod(Res->a, a, global_params.p); }

void fp_set_si(fp *Res, const long int a) {
  mpz_set_si(global_params.tmp1, a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_to_string(const fp x, char *a) { gmp_sprintf(a, "%#Zx", x.a); }

int fp_iszero(const fp x) { return (mpz_cmp_si(x.a, 0) == 0); }

void fp_copy(fp *Res, const fp x) { mpz_set(Res->a, x.a); }

int fp_cmp(const fp x, const fp y) { return mpz_cmp(x.a, y.a); }

int fp_isone(const fp x) { return (mpz_cmp_si(x.a, 1) == 0); }

void fp_random(fp *Res) {
  mpz_urandomm(Res->a, global_params.state, global_params.p);
}

int fp_isequal(const fp x, const fp y) { return (fp_cmp(x, y) == 0); }

void fp_add3(fp *Res, const fp x, const fp y) {
  mpz_add(global_params.tmp1, x.a, y.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_add2(fp *Res, const fp x) {
  mpz_add(global_params.tmp1, Res->a, x.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_add3_si(fp *Res, const fp x, const long int a) {
  mpz_add_ui(global_params.tmp1, x.a, a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_add2_si(fp *Res, const long int a) {
  mpz_add_ui(global_params.tmp1, Res->a, a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sub3(fp *Res, const fp x, const fp y) {
  mpz_sub(global_params.tmp1, x.a, y.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sub2(fp *Res, const fp x) {
  mpz_sub(global_params.tmp1, Res->a, x.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sub3_si(fp *Res, const fp x, const long int a) {
  mpz_sub_ui(global_params.tmp1, x.a, a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sub2_si(fp *Res, const long int a) {
  mpz_sub_ui(global_params.tmp1, Res->a, a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_neg(fp *Res, const fp x) {
  if (mpz_cmp_si(x.a, 0) == 0)
    mpz_set(Res->a, x.a);
  else
    mpz_sub(Res->a, global_params.p, x.a);
}

void fp_mul3(fp *Res, const fp x, const fp y) {
  mpz_mul(global_params.tmp1, x.a, y.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_mul2(fp *Res, const fp x) {
  mpz_mul(global_params.tmp1, Res->a, x.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sqr2(fp *Res, const fp x) {
  mpz_mul(global_params.tmp1, x.a, x.a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_sqr1(fp *Res) {
  mpz_mul(global_params.tmp1, Res->a, Res->a);
  mpz_mod(Res->a, global_params.tmp1, global_params.p);
}

void fp_inv2(fp *Res, const fp x) { mpz_invert(Res->a, x.a, global_params.p); }

void fp_inv1(fp *Res) { mpz_invert(Res->a, Res->a, global_params.p); }

void fp_div3(fp *Res, const fp x, const fp y) {
  fp_inv2(&(global_params.fp_temp1), y);
  fp_mul3(Res, x, global_params.fp_temp1);
}

void fp_div2(fp *Res, const fp x) {
  fp_inv2(&(global_params.fp_temp1), x);
  fp_mul2(Res, global_params.fp_temp1);
}
