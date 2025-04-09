#ifndef COMMON_UTILS_H
#define COMMON_UTILS_H

#include <gmp.h>

void power_mod_gmp(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus);

#endif
