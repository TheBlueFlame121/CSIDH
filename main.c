#include "csidh.h"
#include "fp.h"
#include <stdio.h>

extern const char pStr[];
extern const int primes[NUM_PRIMES];

int main() {
  mpz_t prime;
  mpz_init(prime);
  mpz_set_str(prime, pStr, 0);
  fp_setup(prime);
  fp A_start;
  fp_init(&A_start);

  // Alice round 1
  privKey privA;
  gen_private(&privA);
  printf("Alice private key: ");
  for (int i = 0; i < NUM_PRIMES; i++)
    printf("%d ", privA.e[i]);
  printf("\n");

  fp A_alice;
  fp_init(&A_alice);

  action(&A_alice, A_start, privA);
  gmp_printf("Alice's public key is: %#Zx\n", A_alice.a);
  printf("Valid: %d\n\n\n", validate(A_alice));

  // Bob round 1
  privKey privB;
  gen_private(&privB);
  printf("Bob private key: ");
  for (int i = 0; i < NUM_PRIMES; i++)
    printf("%d ", privB.e[i]);
  printf("\n");

  fp A_bob;
  fp_init(&A_bob);

  action(&A_bob, A_start, privB);
  gmp_printf("Bob's public key is: %#Zx\n", A_bob.a);
  printf("Valid: %d\n\n\n", validate(A_bob));

  // Alice round 2
  fp sec_A;
  fp_init(&sec_A);
  action(&sec_A, A_bob, privA);
  gmp_printf("Shared secret Alice: %#Zx\n", sec_A.a);

  // Bob round 2
  fp sec_B;
  fp_init(&sec_B);
  action(&sec_B, A_alice, privB);
  gmp_printf("Shared secret Bob: %#Zx\n", sec_B.a);

  // check
  if (mpz_cmp(sec_A.a, sec_B.a) == 0)
    printf("\nSUCCESS\n");
  else
    printf("\nFAIL\n");

  return 0;
}
