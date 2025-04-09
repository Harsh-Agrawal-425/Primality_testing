#ifndef RSA_UTILS_H
#define RSA_UTILS_H

#include <gmp.h>

void rsa_keygen(mpz_t n, mpz_t e, mpz_t d, const mpz_t p, const mpz_t q);

#endif
