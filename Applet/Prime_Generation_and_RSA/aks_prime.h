#ifndef AKS_PRIME_H
#define AKS_PRIME_H

#include <gmp.h>

void generate_large_prime_deterministically(mpz_t prime, int bits, int &attempts, double &total_time_ms);

#endif
