#include <iostream>
#include <fstream>
#include <chrono>
#include <iomanip>
#include <gmp.h>
#include <string>
#include <stdexcept>
#include <ctime>

class MillerRabinTester {
private:
    void power_mod_gmp(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus) {
        mpz_powm(result, base, exponent, modulus);
    }

    bool miller_rabin_test_gmp(const mpz_t n, int k) {
        if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0)
            return true;

        if (mpz_cmp_ui(n, 2) < 0 || mpz_even_p(n))
            return false;

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
            mpz_sub_ui(a, n, 3);
            mpz_urandomm(a, state, a);
            mpz_add_ui(a, a, 2);

            power_mod_gmp(x, a, d, n);
            if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0)
                continue;

            bool witness = false;
            for (unsigned long j = 0; j < r - 1; j++) {
                mpz_powm_ui(x, x, 2, n);
                if (mpz_cmp_ui(x, 1) == 0)
                    return false;
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

public:
    bool is_probably_prime(const mpz_t n, int k = 40) {
        return miller_rabin_test_gmp(n, k);
    }
};

std::string get_timestamp() {
    auto now = std::chrono::system_clock::now();
    auto t = std::chrono::system_clock::to_time_t(now);
    std::ostringstream oss;
    oss << std::put_time(std::localtime(&t), "%Y-%m-%d %H:%M:%S");
    return oss.str();
}

int main() {
    std::ofstream outfile("mr_runtime_results.csv");
    outfile << "Timestamp,BitLength,TimeTaken(ms),PrimeNumber\n";

    MillerRabinTester tester;
    gmp_randstate_t rstate;
    gmp_randinit_mt(rstate);
    gmp_randseed_ui(rstate, time(NULL));

    mpz_t rand_num, prime;
    mpz_inits(rand_num, prime, NULL);

    for (int bits = 2; bits <= 4096; ++bits) {
        std::cout << "Testing bit length: " << bits << std::endl;
        for (int i = 0; i < 20; ++i) {
            mpz_urandomb(rand_num, rstate, bits);
            mpz_nextprime(prime, rand_num);  // ensure it's prime

            auto start = std::chrono::high_resolution_clock::now();
            tester.is_probably_prime(prime, 40);
            auto end = std::chrono::high_resolution_clock::now();

            double elapsed = std::chrono::duration<double, std::milli>(end - start).count();
            std::string timestamp = get_timestamp();

            outfile << timestamp << "," << bits << "," << elapsed << "," << mpz_get_str(NULL, 10, prime) << "\n";
        }
        outfile.flush();
    }

    mpz_clears(rand_num, prime, NULL);
    gmp_randclear(rstate);
    outfile.close();

    std::cout << "Done! Results written to mr_runtime_results.csv" << std::endl;
    return 0;
}
