#ifndef MILLER_RABIN_H
#define MILLER_RABIN_H

#include <gmp.h>

bool miller_rabin_test_gmp(const mpz_t n, int k = 40);
void generate_large_prime(mpz_t prime, int bits, int &attempts, double &total_time_ms);

#endif
