#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <gmp.h>
#include <chrono>

int count_bits(const std::string &number_str)
{
    mpz_t num;
    mpz_init(num);
    mpz_set_str(num, number_str.c_str(), 10);

    int bit_count = mpz_sizeinbase(num, 2);

    mpz_clear(num);
    return bit_count;
}

unsigned long find_order(mpz_t n, unsigned long r)
{
    mpz_t mod_result, r_mpz;
    mpz_init(mod_result);
    mpz_init_set_ui(r_mpz, r);

    for (unsigned long k = 1; k <= r; k++)
    {
        mpz_powm_ui(mod_result, n, k, r_mpz);
        if (mpz_cmp_ui(mod_result, 1) == 0)
        {
            mpz_clear(mod_result);
            mpz_clear(r_mpz);
            return k;
        }
    }

    mpz_clear(mod_result);
    mpz_clear(r_mpz);
    return r + 1;
}

int is_perfect_power(mpz_t n)
{
    if (mpz_cmp_ui(n, 1) == 0)
        return 1;

    size_t bit_length = mpz_sizeinbase(n, 2);

    mpz_t root, result;
    mpz_init(root);
    mpz_init(result);

    for (unsigned long b = 2; b <= bit_length; b++)
    {
        mpz_root(root, n, b);

        mpz_pow_ui(result, root, b);
        if (mpz_cmp(result, n) == 0)
        {
            mpz_clear(root);
            mpz_clear(result);
            return 1;
        }
    }

    mpz_clear(root);
    mpz_clear(result);
    return 0;
}

int is_coprime(unsigned long a, unsigned long b)
{
    unsigned long temp;
    while (b != 0)
    {
        temp = b;
        b = a % b;
        a = temp;
    }
    return (a == 1);
}

int polynomial_congruence_test(mpz_t n, unsigned long r)
{
    mpz_t coefficient, term, poly_eval, x_pow;
    mpz_init(coefficient);
    mpz_init(term);
    mpz_init(poly_eval);
    mpz_init(x_pow);

    size_t log_n = mpz_sizeinbase(n, 2);

    unsigned long a_bound = (unsigned long)(sqrt(r) * log_n);
    if (a_bound < 1)
        a_bound = 1;

    for (unsigned long a = 1; a <= a_bound; a++)
    {
        mpz_set_ui(poly_eval, 0);

        for (unsigned long i = 0; i <= r - 1; i++)
        {
            mpz_bin_ui(coefficient, n, i);

            mpz_t n_minus_i;
            mpz_init(n_minus_i);
            mpz_sub_ui(n_minus_i, n, i);

            mpz_t a_mpz;
            mpz_init_set_ui(a_mpz, a);
            mpz_powm(x_pow, a_mpz, n_minus_i, n);

            mpz_mul(term, coefficient, x_pow);
            mpz_mod(term, term, n);

            mpz_add(poly_eval, poly_eval, term);
            mpz_mod(poly_eval, poly_eval, n);

            mpz_clear(n_minus_i);
            mpz_clear(a_mpz);
        }

        mpz_t expected;
        mpz_init(expected);
        mpz_set_ui(expected, a);

        if (mpz_cmp(poly_eval, expected) != 0)
        {
            mpz_clear(coefficient);
            mpz_clear(term);
            mpz_clear(poly_eval);
            mpz_clear(x_pow);
            mpz_clear(expected);
            return 0;
        }

        mpz_clear(expected);
    }

    mpz_clear(coefficient);
    mpz_clear(term);
    mpz_clear(poly_eval);
    mpz_clear(x_pow);
    return 1;
}

int aks_primality_test(mpz_t n)
{
    if (is_perfect_power(n))
    {
        return 0;
    }

    size_t log_n_squared = mpz_sizeinbase(n, 2);
    log_n_squared = log_n_squared * log_n_squared;

    unsigned long r = 2;
    while (r < mpz_get_ui(n))
    {
        if (is_coprime(r, mpz_get_ui(n)))
        {
            unsigned long ord = find_order(n, r);
            if (ord > log_n_squared)
            {
                break;
            }
        }
        r++;
    }

    mpz_t gcd_result;
    mpz_init(gcd_result);

    for (unsigned long a = 2; a <= r; a++)
    {
        mpz_gcd_ui(gcd_result, n, a);
        if (mpz_cmp_ui(gcd_result, 1) > 0 && mpz_cmp(gcd_result, n) < 0)
        {
            mpz_clear(gcd_result);
            return 0;
        }
    }
    mpz_clear(gcd_result);

    if (mpz_cmp_ui(n, r) <= 0)
    {
        return 1;
    }

    return polynomial_congruence_test(n, r);
}

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        std::cerr << "Usage: " << argv[0] << " <input_csv_file> <output_csv_file>" << std::endl;
        return 1;
    }

    std::string input_file = argv[1];
    std::string output_file = argv[2];

    std::ifstream infile(input_file);
    if (!infile.is_open())
    {
        std::cerr << "Failed to open input file: " << input_file << std::endl;
        return 1;
    }

    std::ofstream outfile(output_file);
    if (!outfile.is_open())
    {
        std::cerr << "Failed to open output file: " << output_file << std::endl;
        infile.close();
        return 1;
    }

    std::string line;

    auto total_start_time = std::chrono::high_resolution_clock::now();

    if (std::getline(infile, line))
    {
        outfile << line << ",Num_Bits,Determ_Output,AKS_Time" << std::endl;
    }

    long long total_numbers = 0;
    long long conflict_count = 0;
    long long prime_count_mr = 0;  
    long long prime_count_aks = 0; 

    while (std::getline(infile, line))
    {
        std::istringstream iss(line);
        std::string timestamp, number_str, time_taken, primality;

        if (std::getline(iss, timestamp, ',') &&
            std::getline(iss, number_str, ',') &&
            std::getline(iss, time_taken, ',') &&
            std::getline(iss, primality))
        {

            total_numbers++;

            int bits = count_bits(number_str);

            std::string determ_output;

            if (primality == "Composite")
            {
                determ_output = "Composite";
                outfile << line << "," << bits << "," << determ_output << ",NA" << std::endl;
            }
            else if (primality == "Prime")
            {
                prime_count_mr++;

                mpz_t num;
                mpz_init(num);
                mpz_set_str(num, number_str.c_str(), 10);

                auto start_time = std::chrono::high_resolution_clock::now();

                int is_prime = aks_primality_test(num);

                auto end_time = std::chrono::high_resolution_clock::now();
                std::chrono::duration<double> elapsed = end_time - start_time;
                double aks_time = elapsed.count();

                if (is_prime)
                {
                    determ_output = "Prime";
                    prime_count_aks++;
                }
                else
                {
                    determ_output = "Composite";
                    conflict_count++;
                    std::cout << "Conflict found for number: " << number_str << std::endl;
                }

                mpz_clear(num);

                outfile << line << "," << bits << "," << determ_output << "," << aks_time << std::endl;
            }
            else
            {
                determ_output = "Unknown";
                outfile << line << "," << bits << "," << determ_output << ",NA" << std::endl;
            }

            if (total_numbers % 1000 == 0)
            {
                std::cout << "Processed " << total_numbers << " numbers so far..." << std::endl;
            }
        }
    }

    auto total_end_time = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> total_elapsed = total_end_time - total_start_time;
    double total_runtime = total_elapsed.count();

    infile.close();
    outfile.close();

    std::cout << "Processing complete!" << std::endl;
    std::cout << "Total numbers processed: " << total_numbers << std::endl;
    std::cout << "Numbers declared Prime by Miller-Rabin: " << prime_count_mr << std::endl;
    std::cout << "Numbers confirmed Prime by AKS: " << prime_count_aks << std::endl;
    std::cout << "Conflicts found (MR=Prime, AKS=Composite): " << conflict_count
              << " out of " << total_numbers << " numbers ("
              << (double)conflict_count / total_numbers * 100 << "%)" << std::endl;
    std::cout << "Total execution time: " << total_runtime << " seconds" << std::endl;

    std::ofstream log_file(output_file, std::ios::app);
    if (log_file.is_open())
    {
        log_file << "Total Execution Time (s)," << total_runtime << std::endl;
        log_file.close();
    }

    return 0;
}
