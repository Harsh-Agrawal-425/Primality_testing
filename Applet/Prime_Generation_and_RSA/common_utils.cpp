#include "common_utils.h"
#include <gmp.h>

void power_mod_gmp(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus) {
    mpz_powm(result, base, exponent, modulus);
}
