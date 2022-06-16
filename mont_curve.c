#include "mont_curve.h"

void MP_init(MP *Res) {
  fp_init(&(Res->X));
  fp_init(&(Res->Z));
}

void MP_set(MP *Res, const fp x, const fp z) {
  fp_copy(&(Res->X), x);
  fp_copy(&(Res->Z), z);
}

void MP_set_mpz(MP *Res, const mpz_t x, const mpz_t z) {
  fp_set_mpz(&(Res->X), x);
  fp_set_mpz(&(Res->Z), z);
}

void MP_set_si(MP *Res, const long int x, const long int z) {
  fp_set_si(&(Res->X), x);
  fp_set_si(&(Res->Z), z);
}

void MP_set_x(MP *Res, const fp x) {
  fp_copy(&(Res->X), x);
  fp_set_si(&(Res->Z), 1);
}

void MP_copy(MP *Res, const MP P) {
  fp_copy(&(Res->X), P.X);
  fp_copy(&(Res->Z), P.Z);
}

void MP_clear(MP *Res) {
  fp_clear(&(Res->X));
  fp_clear(&(Res->Z));
}

void MP_norm(MP *Res, const MP P) {
  fp_div3(&(Res->X), P.X, P.Z);
  fp_set_si(&(Res->Z), 1);
}

void xADD(MP *Res, const MP P, const MP Q, const MP R) {
  fp_add3(&(global_params.t0), Q.X, Q.Z);
  fp_sub3(&(global_params.t1), Q.X, Q.Z);
  fp_add3(&(global_params.t2), P.X, P.Z);
  fp_sub3(&(global_params.t3), P.X, P.Z);
  fp_mul2(&(global_params.t3), global_params.t0);
  fp_mul2(&(global_params.t2), global_params.t1);
  fp_add3(&(global_params.t0), global_params.t3, global_params.t2);
  fp_sub3(&(global_params.t1), global_params.t3, global_params.t2);
  fp_sqr1(&(global_params.t0));
  fp_sqr1(&(global_params.t1));
  fp_mul3(&(Res->X), R.Z, global_params.t0);
  fp_mul3(&(Res->Z), R.X, global_params.t1);
  return;
}

void xDBL(MP *Res, const MP P, const MP Ap24_C24) {
  fp_sub3(&(global_params.t0), P.X, P.Z);
  fp_add3(&(global_params.t1), P.X, P.Z);
  fp_sqr1(&(global_params.t0));

  fp_sqr1(&(global_params.t1));
  fp_mul3(&(Res->Z), Ap24_C24.Z, global_params.t0);
  fp_mul3(&(Res->X), Res->Z, global_params.t1);

  fp_sub2(&(global_params.t1), global_params.t0);
  fp_mul3(&(global_params.t0), Ap24_C24.X, global_params.t1);
  fp_add2(&(Res->Z), global_params.t0);

  fp_mul2(&(Res->Z), global_params.t1);
  return;
}

void xDBLADD(MP *Dub, MP *Add, const MP P, const MP Q, const MP R,
             const MP Ap24_1) {
  fp_add3(&(global_params.t0), P.X, P.Z);
  fp_sub3(&(global_params.t1), P.X, P.Z);
  fp_sqr2(&(Dub->X), global_params.t0);
  fp_sub3(&(global_params.t2), Q.X, Q.Z);
  fp_add3(&(Add->X), Q.X, Q.Z);
  fp_mul2(&(global_params.t0), global_params.t2);
  fp_sqr2(&(Dub->Z), global_params.t1);

  fp_mul2(&(global_params.t1), Add->X);
  fp_sub3(&(global_params.t2), Dub->X, Dub->Z);
  fp_mul2(&(Dub->X), Dub->Z);
  fp_mul3(&(Add->X), Ap24_1.X, global_params.t2);
  fp_sub3(&(Add->Z), global_params.t0, global_params.t1);
  fp_add2(&(Dub->Z), Add->X);
  fp_add3(&(Add->X), global_params.t0, global_params.t1);

  fp_mul2(&(Dub->Z), global_params.t2);
  fp_sqr1(&(Add->Z));
  fp_sqr1(&(Add->X));
  fp_mul2(&(Add->Z), R.X);
  fp_mul2(&(Add->X), R.Z);
  return;
}

void xMUL(MP *Res, const MP P, mpz_t k, const MP A_1) {
  MP x0, x1, Ap24_1, MP_temp1, MP_temp2;
  MP_init(&x0);
  MP_init(&x1);
  MP_init(&Ap24_1);
  MP_init(&MP_temp1);
  MP_init(&MP_temp2);

  MP_copy(&x0, P);
  fp_set_si(&(Ap24_1.Z), 1);
  fp_add3_si(&(Ap24_1.X), A_1.X, 2);
  fp_set_si(&(global_params.t0), 4);
  fp_div2(&(Ap24_1.X), global_params.t0);

  xDBL(&x1, P, Ap24_1);

  for (int i = mpz_sizeinbase(k, 2) - 2; i >= 0; i--) {
    if (mpz_tstbit(k, i) == 0) {
      xDBLADD(&MP_temp1, &MP_temp2, x0, x1, P, Ap24_1);
      MP_copy(&x0, MP_temp1);
      MP_copy(&x1, MP_temp2);
    } else {
      xDBLADD(&MP_temp1, &MP_temp2, x1, x0, P, Ap24_1);
      MP_copy(&x0, MP_temp2);
      MP_copy(&x1, MP_temp1);
    }
  }

  MP_copy(Res, x0);

  MP_clear(&x0);
  MP_clear(&x1);
  MP_clear(&Ap24_1);
  MP_clear(&MP_temp1);
  MP_clear(&MP_temp2);

  return;
}

void mont_rhs(fp *Res, const MP A, const fp x) {
  fp_copy(Res, x);    // x
  fp_add2(Res, A.X);  // x+A
  fp_mul2(Res, x);    // x^2 + Ax
  fp_add2_si(Res, 1); // x^2 + Ax + 1
  fp_mul2(Res, x);    // x^3 + Ax^2 + x
}

void xISOG(MP *Aout, MP *Pout, const MP A, const MP P, const MP K, int k) {
  assert(k >= 3);
  assert(k % 2 == 1);

  MP_copy(Aout, A);
  MP_copy(Pout, P);

  fp tmp0, tmp1, tmp2, Psum, Pdif;
  fp_init(&tmp0);
  fp_init(&tmp1);
  fp_init(&tmp2);
  fp_init(&Psum);
  fp_init(&Pdif);

  MP Q, Aed, prod, Ap24_C24;
  MP_init(&Q);
  MP_init(&Aed);
  MP_init(&prod);
  MP_init(&Ap24_C24);

  fp_add3_si(&(Ap24_C24.X), Aout->X, 2);
  fp_set_si(&(Ap24_C24.Z), 4);

  // compute twisted Edwards curve coefficients
  fp_add3(&(Aed.Z), Aout->Z, Aout->Z);
  fp_add3(&(Aed.X), Aout->X, Aed.Z);
  fp_sub3(&(Aed.Z), Aout->X, Aed.Z);

  fp_add3(&Psum, P.X, P.Z); // precomputations
  fp_sub3(&Pdif, P.X, P.Z);

  fp_sub3(&(prod.X), K.X, K.Z);
  fp_add3(&(prod.Z), K.X, K.Z);

  fp_mul3(&tmp1, prod.X, Psum);
  fp_mul3(&tmp0, prod.Z, Pdif);
  fp_add3(&(Q.X), tmp0, tmp1);
  fp_sub3(&(Q.Z), tmp0, tmp1);

  // MP M[3] = {K};
  MP M[3];
  for (int i = 0; i < 3; i++)
    MP_init(&M[i]);
  MP_copy(&M[0], K);

  xDBL(&M[1], K, Ap24_C24);

  for (int i = 1; i < k / 2; ++i) {

    if (i >= 2)
      xADD(&M[i % 3], M[(i - 1) % 3], K, M[(i - 2) % 3]);

    fp_sub3(&tmp1, M[i % 3].X, M[i % 3].Z);
    fp_add3(&tmp0, M[i % 3].X, M[i % 3].Z);
    fp_mul2(&(prod.X), tmp1);
    fp_mul2(&(prod.Z), tmp0);
    fp_mul2(&tmp1, Psum);
    fp_mul2(&tmp0, Pdif);
    fp_add3(&tmp2, tmp0, tmp1);
    fp_mul2(&(Q.X), tmp2);
    fp_sub3(&tmp2, tmp0, tmp1);
    fp_mul2(&(Q.Z), tmp2);
  }

  // point evaluation
  fp_sqr1(&(Q.X));
  fp_sqr1(&(Q.Z));
  fp_mul3(&(Pout->X), P.X, Q.X);
  fp_mul3(&(Pout->Z), P.Z, Q.Z);

  // compute Aed.x^k, Aed.z^k
  /* exp_by_squaring_(&Aed.x, &Aed.z, k); */
  mpz_powm_ui(Aed.X.a, Aed.X.a, k, global_params.p);
  mpz_powm_ui(Aed.Z.a, Aed.Z.a, k, global_params.p);

  // compute prod.x^8, prod.z^8
  fp_sqr1(&(prod.X));
  fp_sqr1(&(prod.X));
  fp_sqr1(&(prod.X));
  fp_sqr1(&(prod.Z));
  fp_sqr1(&(prod.Z));
  fp_sqr1(&(prod.Z));

  // compute image curve parameters
  fp_mul2(&(Aed.Z), prod.X);
  fp_mul2(&(Aed.X), prod.Z);

  // compute Montgomery params
  fp_add3(&(Aout->X), Aed.X, Aed.Z);
  fp_sub3(&(Aout->Z), Aed.X, Aed.Z);
  fp_add2(&(Aout->X), Aout->X);

  fp_clear(&tmp0);
  fp_clear(&tmp1);
  fp_clear(&tmp2);
  fp_clear(&Psum);
  fp_clear(&Pdif);

  MP_clear(&Q);
  MP_clear(&Aed);
  MP_clear(&prod);
  MP_clear(&Ap24_C24);
}
