#include "csidh.h"

/* const int primes[NUM_PRIMES] = {3,  5,  7,  11, 13, 17, 19, */
/*                                 23, 29, 31, 37, 41, 43, 97}; */

/* const char pStr[] = "0x2338fc35dff3029b"; */

/* const int max_exp = 2; */

const unsigned primes[NUM_PRIMES] = {
    3,   5,   7,   11,  13,  17,  19,  23,  29,  31,  37,  41,  43,  47,  53,
    59,  61,  67,  71,  73,  79,  83,  89,  97,  101, 103, 107, 109, 113, 127,
    131, 137, 139, 149, 151, 157, 163, 167, 173, 179, 181, 191, 193, 197, 199,
    211, 223, 227, 229, 233, 239, 241, 251, 257, 263, 269, 271, 277, 281, 283,
    293, 307, 311, 313, 317, 331, 337, 347, 349, 353, 359, 367, 373, 587};

const char pStr[] =
    "0x65b48e8f740f89bffc8ab0d15e3e4c4ab42d083aedc88c425afbfcc69322c9cda7aac6c5"
    "67f35507516730cc1f0b4f25c2721bf457aca8351b81b90533c6c87b";

const int max_exp = 5;

int sign(int x) { return (x > 0) - (x < 0); }

void gen_private(privKey *K) {
  for (int i = 0; i < NUM_PRIMES; i++)
    K->e[i] = (rand() % (2 * max_exp + 1)) - max_exp;
}

int validate(const fp A) {
  mpz_t mul, r2p, d;
  mpz_init(mul);
  mpz_init(r2p);
  mpz_init(d);

  mpz_sqrt(r2p, global_params.p);
  mpz_mul_si(r2p, r2p, 4);

  mpz_set_si(d, 1);

  fp x;
  fp_init(&x);
  fp_random(&x);

  MP P, Q, MP_temp, A_1;
  MP_init(&P);
  MP_init(&Q);
  MP_init(&MP_temp);
  MP_init(&A_1);

  MP_set_x(&A_1, A);
  MP_set_x(&P, x);

  int flag = 0;
  for (int i = 0; i < NUM_PRIMES; i++) {
    mpz_add_ui(mul, global_params.p, 1);
    mpz_div_ui(mul, mul, primes[i]);

    xMUL(&MP_temp, P, mul, A_1);
    MP_copy(&Q, MP_temp);

    mpz_set_si(mul, primes[i]);
    xMUL(&MP_temp, Q, mul, A_1);
    if (!fp_iszero(MP_temp.Z)) {
      flag = 0;
      break;
    }

    if (!fp_iszero(Q.Z))
      mpz_mul_si(d, d, primes[i]);

    if (mpz_cmp(d, r2p) > 0) {
      flag = 1;
      break;
    }
  }

  mpz_clear(mul);
  mpz_clear(r2p);
  mpz_clear(d);

  fp_clear(&x);

  MP_clear(&P);
  MP_clear(&Q);
  MP_clear(&MP_temp);
  MP_clear(&A_1);
  return flag;
}

void action(fp *Aout, const fp A, const privKey Key) {
  fp x, rhs;
  fp_init(&x);
  fp_init(&rhs);

  mpz_t k, p_mul, q_mul;
  mpz_init(k);
  mpz_init(p_mul);
  mpz_init(q_mul);

  MP P, Q, R, MP_temp1, MP_temp2;
  MP_init(&P);
  MP_init(&Q);
  MP_init(&R);
  MP_init(&MP_temp1);
  MP_init(&MP_temp2);

  int S[NUM_PRIMES], s;

  MP MP_A;
  MP_init(&MP_A);
  MP_set_x(&MP_A, A);

  int e[NUM_PRIMES];
  for (int i = 0; i < NUM_PRIMES; i++) {
    e[i] = Key.e[i];
  }

  int flag = 0;
  for (int i = 0; i < NUM_PRIMES; i++) {
    if (e[i] != 0) {
      flag = 1;
      break;
    }
  }
  while (flag) {
    fp_random(&x);
    mont_rhs(&rhs, MP_A, x);

    s = mpz_kronecker(rhs.a, global_params.p);
    if (s == 0)
      continue;

    mpz_set_si(k, 1);
    memset(S, 0, sizeof S);
    for (int i = 0; i < NUM_PRIMES; i++) {
      if (sign(e[i]) == s) {
        S[i] = 1;
        mpz_mul_si(k, k, primes[i]);
      }
    }
    if (mpz_cmp_si(k, 1) == 0)
      continue;

    mpz_add_ui(p_mul, global_params.p, 1);
    mpz_div(p_mul, p_mul, k);

    MP_set_x(&P, x);
    xMUL(&Q, P, p_mul, MP_A);

    for (int i = 0; i < NUM_PRIMES; i++) {
      if (S[i] == 0)
        continue;

      mpz_div_ui(q_mul, k, primes[i]);
      xMUL(&R, Q, q_mul, MP_A);
      if (fp_iszero(R.Z))
        continue;

      xISOG(&MP_temp1, &MP_temp2, MP_A, Q, R, primes[i]);
      MP_norm(&MP_A, MP_temp1);
      MP_copy(&Q, MP_temp2);

      e[i] -= s;
    }

    flag = 0;
    for (int i = 0; i < NUM_PRIMES; i++) {
      if (e[i] != 0) {
        flag = 1;
        break;
      }
    }
  }

  fp_copy(Aout, MP_A.X);

  MP_clear(&MP_A);
  MP_clear(&P);
  MP_clear(&Q);
  MP_clear(&R);
  MP_clear(&MP_temp1);
  MP_clear(&MP_temp2);

  fp_clear(&x);
  fp_clear(&rhs);

  mpz_clear(k);
  mpz_clear(p_mul);
  mpz_clear(q_mul);

  return;
}
