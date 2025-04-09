#include "rsa_utils.h"
#include <gmp.h>

void rsa_keygen(mpz_t n, mpz_t e, mpz_t d, const mpz_t p, const mpz_t q) {
    mpz_t phi, p1, q1, gcd;
    mpz_inits(phi, p1, q1, gcd, NULL);

    mpz_mul(n, p, q);
    mpz_sub_ui(p1, p, 1);
    mpz_sub_ui(q1, q, 1);
    mpz_mul(phi, p1, q1);

    mpz_set_ui(e, 65537);
    mpz_gcd(gcd, e, phi);
    while (mpz_cmp_ui(gcd, 1) != 0) {
        mpz_add_ui(e, e, 2);
        mpz_gcd(gcd, e, phi);
    }

    mpz_invert(d, e, phi);
    mpz_clears(phi, p1, q1, gcd, NULL);
}

void rsa_encrypt(mpz_t ciphertext, const mpz_t message, const mpz_t e, const mpz_t n) {
    mpz_powm(ciphertext, message, e, n);
}

void rsa_decrypt(mpz_t plaintext, const mpz_t ciphertext, const mpz_t d, const mpz_t n) {
    mpz_powm(plaintext, ciphertext, d, n);
}
