#include <iostream>
#include <gmp.h>
#include <chrono>
#include <cstdlib>
#include "miller_rabin.h"
#include "aks_prime.h"
#include "rsa_utils.h"

using namespace std;

int main(int argc, char **argv)
{
    int bits = 40;
    if (argc > 1)
    {
        bits = atoi(argv[1]);
        if (bits <= 0)
        {
            cerr << "Invalid bit size provided." << endl;
            return 1;
        }
    }

    mpz_t prime, aksprime;
    mpz_inits(prime, aksprime, NULL);

    int attempts;
    double total_time_ms;

    // Miller-Rabin Prime
    generate_large_prime(prime, bits, attempts, total_time_ms);
    mpz_out_str(stdout, 10, prime); cout << endl;
    cout << attempts << endl;
    cout << total_time_ms << endl;
    cout << (total_time_ms / attempts) << endl;

    // AKS Prime
    int aksAttempts = 0;
    double total_time_ms_aks = 0;
    generate_large_prime_deterministically(aksprime, bits, aksAttempts, total_time_ms_aks);
    mpz_out_str(stdout, 10, aksprime); cout << endl;
    cout << aksAttempts << endl;
    cout << total_time_ms_aks << endl;
    cout << (total_time_ms_aks / aksAttempts) << endl;

    // RSA Setup (MR)
    mpz_t mr_p, mr_q, mr_n, mr_e, mr_d;
    mpz_inits(mr_p, mr_q, mr_n, mr_e, mr_d, NULL);
    auto start_keygen_mr = chrono::high_resolution_clock::now();
    generate_large_prime(mr_p, bits / 2, attempts, total_time_ms);
    generate_large_prime(mr_q, bits / 2, attempts, total_time_ms);
    rsa_keygen(mr_n, mr_e, mr_d, mr_p, mr_q);
    auto end_keygen_mr = chrono::high_resolution_clock::now();
    double keygen_time_mr = chrono::duration<double, milli>(end_keygen_mr - start_keygen_mr).count();

    cout << "\n[Miller-Rabin RSA]" << endl;
    cout << "Key Generation Time (ms): " << keygen_time_mr << endl;
    cout << endl;

    // RSA Setup (AKS)
    mpz_t aks_p, aks_q, aks_n, aks_e, aks_d;
    mpz_inits(aks_p, aks_q, aks_n, aks_e, aks_d, NULL);
    auto start_keygen_aks = chrono::high_resolution_clock::now();
    generate_large_prime_deterministically(aks_p, bits / 2, aksAttempts, total_time_ms_aks);
    generate_large_prime_deterministically(aks_q, bits / 2, aksAttempts, total_time_ms_aks);
    rsa_keygen(aks_n, aks_e, aks_d, aks_p, aks_q);
    auto end_keygen_aks = chrono::high_resolution_clock::now();
    double keygen_time_aks = chrono::duration<double, milli>(end_keygen_aks - start_keygen_aks).count();

    cout << "\n[AKS RSA]" << endl;
    cout << "Key Generation Time (ms): " << keygen_time_aks << endl;
    cout << endl;

    mpz_clears(prime, aksprime,
               mr_p, mr_q, mr_n, mr_e, mr_d,
               aks_p, aks_q, aks_n, aks_e, aks_d,
               NULL);

    cout.flush();
    return 0;
}
