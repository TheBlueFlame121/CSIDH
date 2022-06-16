#ifndef CSIDH_H_
#define CSIDH_H_

#include "mont_curve.h"
#include <string.h>

/* #define NUM_PRIMES 14 */
#define NUM_PRIMES 74

typedef struct {
  int e[NUM_PRIMES];
} privKey;

int sign(int x);

void gen_private(privKey *K);

int validate(const fp A);

void action(fp *Aout, const fp A, const privKey Key);

#endif // CSIDH_H_
