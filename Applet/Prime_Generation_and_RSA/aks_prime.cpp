#include "aks_prime.h"
#include "aks.h"
#include <chrono>
#include <gmp.h>

void generate_large_prime_deterministically(mpz_t prime, int bits, int &attempts, double &total_time_ms) {
    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, std::chrono::system_clock::now().time_since_epoch().count());

    mpz_t lower_bound;
    mpz_init(lower_bound);
    mpz_ui_pow_ui(lower_bound, 2, bits - 1);

    attempts = 0;
    auto start_time = std::chrono::high_resolution_clock::now();
    do {
        attempts++;
        mpz_urandomb(prime, state, bits);
        mpz_add(prime, prime, lower_bound);
    } while (!is_prime_deterministic(conv<ZZ>(mpz_get_str(NULL, 10, prime))));
    auto end_time = std::chrono::high_resolution_clock::now();

    total_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    gmp_randclear(state);
    mpz_clear(lower_bound);
}
