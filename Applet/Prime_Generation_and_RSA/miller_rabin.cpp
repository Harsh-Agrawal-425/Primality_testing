#include "miller_rabin.h"
#include "common_utils.h"
#include <chrono>
#include <gmp.h>

bool miller_rabin_test_gmp(const mpz_t n, int k) {
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0) return true;
    if (mpz_cmp_ui(n, 2) < 0 || mpz_even_p(n)) return false;

    mpz_t d, n_minus_1, a, x;
    mpz_inits(d, n_minus_1, a, x, NULL);
    mpz_sub_ui(n_minus_1, n, 1);
    mpz_set(d, n_minus_1);

    unsigned long r = 0;
    while (mpz_even_p(d)) {
        mpz_tdiv_q_2exp(d, d, 1);
        r++;
    }

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, std::chrono::system_clock::now().time_since_epoch().count());

    for (int i = 0; i < k; i++) {
        mpz_urandomm(a, state, n);
        mpz_add_ui(a, a, 2);
        power_mod_gmp(x, a, d, n);
        if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0) continue;

        bool witness = false;
        for (unsigned long j = 0; j < r - 1; j++) {
            mpz_powm_ui(x, x, 2, n);
            if (mpz_cmp_ui(x, 1) == 0) {
                mpz_clears(d, n_minus_1, a, x, NULL);
                gmp_randclear(state);
                return false;
            }
            if (mpz_cmp(x, n_minus_1) == 0) {
                witness = true;
                break;
            }
        }
        if (!witness) {
            mpz_clears(d, n_minus_1, a, x, NULL);
            gmp_randclear(state);
            return false;
        }
    }

    mpz_clears(d, n_minus_1, a, x, NULL);
    gmp_randclear(state);
    return true;
}

void generate_large_prime(mpz_t prime, int bits, int &attempts, double &total_time_ms) {
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
    } while (!miller_rabin_test_gmp(prime, 40));
    auto end_time = std::chrono::high_resolution_clock::now();

    total_time_ms = std::chrono::duration<double, std::milli>(end_time - start_time).count();

    gmp_randclear(state);
    mpz_clear(lower_bound);
}
