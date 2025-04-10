#include <iostream>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <string>
#include <gmp.h>
#include <chrono>
#include <iomanip>
#include <mutex>
#include <thread>

using namespace std;

// mutex file_mutex; // Ensures thread-safe writing

// Load prime numbers into a set
unordered_set<unsigned long long> load_prime_numbers(const string &filename) {
    unordered_set<unsigned long long> prime_set;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(EXIT_FAILURE);
    }
    unsigned long long num;
    while (file >> num) {
        prime_set.insert(num);
    }
    file.close();
    return prime_set;
}

class MillerRabinTester {
private:
    void power_mod_gmp(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus) {
        mpz_powm(result, base, exponent, modulus);
    }

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
        gmp_randseed_ui(state, chrono::system_clock::now().time_since_epoch().count());

        for (int i = 0; i < k; i++) {
            mpz_sub_ui(a, n, 3);
            mpz_urandomm(a, state, a);
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

public:
    bool is_probably_prime(unsigned long long n, int k = 40) {
        mpz_t mpz_n;
        mpz_init_set_ui(mpz_n, n);
        bool result = miller_rabin_test_gmp(mpz_n, k);
        mpz_clear(mpz_n);
        return result;
    }

    void run_test_with_conflict_check(const string &primes_file, const string &log_file, int k) {
        unordered_set<unsigned long long> prime_set = load_prime_numbers(primes_file);
        
        // Open file once (append mode) and use mutex for thread safety
        ofstream outfile(log_file, ios::app);
        if (!outfile) {
            cerr << "Error: Could not open log file" << endl;
            return;
        }

        // Write CSV header if the file is empty
        if (outfile.tellp() == 0) {
            outfile << "Timestamp,Number,Time_Taken,Actual,Predicted\n";
        }

        int incorrect_count = 0;
        for (unsigned long long n = 2; n <= 100000000; ++n) {
            auto start_time = chrono::high_resolution_clock::now();
            bool is_mr_prime = is_probably_prime(n, k);
            auto end_time = chrono::high_resolution_clock::now();
            double time_taken = chrono::duration<double, milli>(end_time - start_time).count();

            bool is_actual_prime = (prime_set.find(n) != prime_set.end());

            if (is_mr_prime != is_actual_prime) {
                auto now = chrono::system_clock::now();
                auto timestamp = chrono::system_clock::to_time_t(now);

                // lock_guard<mutex> lock(file_mutex); // Prevent concurrent writes
                outfile << put_time(localtime(&timestamp), "%Y-%m-%d %H:%M:%S")
                        << "," << n
                        << "," << time_taken
                        << "," << (is_actual_prime ? "Prime" : "Composite")
                        << "," << (is_mr_prime ? "Prime" : "Composite")
                        << endl;
                incorrect_count++;
            }

            if (n % 1000000 == 0) {
                cout << "Processed up to " << n << ", conflicts found: " << incorrect_count << endl;
                outfile.flush();
            }
        }

        outfile.close();
        cout << "Total incorrect classifications for k = " << k << " : " << incorrect_count << endl;
    }
};

int main(int argc, char *argv[]) {
    if (argc != 4) {
        cerr << "Usage: " << argv[0] << " <primes_file> <log_file> <k>" << endl;
        return EXIT_FAILURE;
    }

    string primes_file = argv[1];
    string log_file = argv[2];
    int k = stoi(argv[3]);

    MillerRabinTester tester;
    tester.run_test_with_conflict_check(primes_file, log_file, k);

    return 0;
}
