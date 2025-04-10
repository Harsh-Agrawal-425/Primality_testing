#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <unordered_set>
#include <chrono>
#include <iomanip>
#include <gmp.h>
#include <gmpxx.h>

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
            if (mpz_cmp_ui(n, 4) <= 0) {
                mpz_clears(d, n_minus_1, a, x, NULL);
                gmp_randclear(state);
                return true;
            }

            mpz_sub_ui(a, n, 3);
            mpz_urandomm(a, state, a);
            mpz_add_ui(a, a, 2);

            power_mod_gmp(x, a, d, n);
            if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0)
                continue;

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

public:
    bool is_probably_prime(unsigned long long n, int k = 40) {
        mpz_t mpz_n;
        mpz_init_set_ui(mpz_n, n);
        bool result = miller_rabin_test_gmp(mpz_n, k);
        mpz_clear(mpz_n);
        return result;
    }
};

void load_primes(const std::string& filename, std::unordered_set<unsigned long long>& prime_set) {
    std::ifstream infile(filename);
    if (!infile)
        throw std::runtime_error("Could not open prime file");

    std::string line;
    while (std::getline(infile, line)) {
        std::stringstream ss(line);
        unsigned long long num;
        while (ss >> num)
            prime_set.insert(num);
    }
}

int bit_length(unsigned long long n) {
    int length = 0;
    while (n) {
        n >>= 1;
        length++;
    }
    return length;
}

int main() {
    try {
        std::unordered_set<unsigned long long> prime_set;
        load_primes("./primes1_till_1e8.txt", prime_set);

        MillerRabinTester tester;
        std::ofstream outfile("sophie_germain_primes.csv");
        outfile << "Start_Timestamp,End_Timestamp,Prime,Bit_Length_Prime,2i+1,Bit_Length_2i+1\n";

        const unsigned long long LIMIT = 500; // i such that 2*i+1 <= 1e8
        for (unsigned long long i = 2; i <= LIMIT; ++i) {
            auto start = std::chrono::system_clock::now();

            if (!tester.is_probably_prime(i)) continue;
            if (prime_set.find(i) == prime_set.end()) continue;

            unsigned long long candidate = 2 * i + 1;
            if (!tester.is_probably_prime(candidate)) continue;
            if (prime_set.find(candidate) == prime_set.end()) continue;

            auto end = std::chrono::system_clock::now();

            std::time_t start_ts = std::chrono::system_clock::to_time_t(start);
            std::time_t end_ts = std::chrono::system_clock::to_time_t(end);

            outfile << std::put_time(std::localtime(&start_ts), "%Y-%m-%d %H:%M:%S") << ","
                    << std::put_time(std::localtime(&end_ts), "%Y-%m-%d %H:%M:%S") << ","
                    << i << "," << bit_length(i) << ","
                    << candidate << "," << bit_length(candidate) << "\n";

            if (i % 1000000 == 0)
                std::cout << "Checked up to i = " << i << std::endl;
        }

        outfile.close();
        std::cout << "Sophie Germain prime search complete.\n";
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
