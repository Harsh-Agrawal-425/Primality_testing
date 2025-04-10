#include <unordered_set>
#include <iomanip>  
#include <iostream>
#include <vector>
#include <chrono>
#include <gmp.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <cmath>
#include <ctime>
#include <fstream>
using namespace std;
using namespace NTL;

// Function declarations for AKS
int step0(ZZ n);
int step1(ZZ n);
long step2(ZZ n);
int step3(ZZ n, long r);
int step4(ZZ n, long r);
int step5(ZZ n, long r);
int step6();
ZZ gcd(ZZ m, long n);
ZZ order(ZZ a, long r);
ZZ phi(long x);

// Main AKS primality test function
int is_prime_deterministic(ZZ n) {
    int k;
    long r;

    k = step0(n);
    if (k) return k == 2;

    if (step1(n)) return 0;
    r = step2(n);
    if (step3(n, r)) return 0;
    if (step4(n, r)) return 1;
    if (step5(n, r)) return 0;
    return step6();
}

// Extra functions
ZZ gcd(ZZ m, long n) {
    ZZ z = n + ZZ(0);
    if (z < m) swap(z, m);
    while (m % z != 0) {
        ZZ k = m % z;
        m = z;
        z = k;
    }
    return z;
}

ZZ order(ZZ a, long r) {
    ZZ o = ZZ(1), k = a, z = ZZ(r);
    while (k != 1) {
        k *= a; k %= z; o++;
    }
    return o;
}

ZZ phi(long x) {
    ZZ n = ZZ(x), ph = ZZ(1), p;
    int k;
    for (PrimeSeq s; n > 1; ) {
        p = s.next();
        k = 0;
        while (n % p == 0) {
            k++; n /= p;
        }
        if (k > 0)
            ph *= power(p, k - 1) * (p - 1);
    }
    return ph;
}

// Steps of AKS
int step0(ZZ n) {
    long p;
    PrimeSeq s;
    while ((p = s.next()) && p < 2000) {
        if ((n % p) == 0)
            return (1 + (n == p));
    }
    return 0;
}

int step1(ZZ n) {
    RR A, B, N = to_RR(n), log_2000_n = log(N) / log(to_RR(2000));
    long b = 2;
    while (b < log_2000_n) {
        B = 1.0 / b;
        A = pow(N, B);
        if ((A - to_RR(FloorToZZ(A))) == 0.0) {
            if (power(FloorToZZ(A), b) == n)
                return 1;
        }
        b++;
    }
    return 0;
}

long step2(ZZ n) {
    RR log2n = to_RR(NumBits(n)), l = power(log2n, 2);
    ZZ k = FloorToZZ(l);
    for (long r = 2;; r++) {
        if (gcd(n, r) == 1) {
            if (order(n, r) > k) return r;
        }
    }
}

int step3(ZZ n, long r) {
    for (long a = 2003; a < r; a++) {
        if ((gcd(n, a) % n) > 1)
            return 1;
    }
    return 0;
}

int step4(ZZ n, long r) {
    return n <= r ? 1 : 0;
}

int step5(ZZ n, long r) {
    ZZ l = FloorToZZ(sqrt(to_RR(phi(r))) * NumBits(n));
    NTL::ZZ_p::init(n);
    ZZ_pX polymod(r, 1); polymod -= 1;
    ZZ_pXModulus mod(polymod);
    ZZ_pX RHS(1, 1); PowerMod(RHS, RHS, n, mod);

    for (long a = 1; a <= l; a++) {
        ZZ_pX LHS(1, 1); LHS += a;
        PowerMod(LHS, LHS, n, mod);
        LHS -= a;
        if (LHS != RHS) return 1;
    }
    return 0;
}

int step6() {
    return 1;
}

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

int main() {
    try {
        std::unordered_set<unsigned long long> prime_set;
        load_primes("../primes1_till_1e8.txt", prime_set);

        std::ofstream outfile("sophie_germain_primes_aks.csv");
        outfile << "Start_Timestamp,End_Timestamp,Prime,Bit_Length_Prime,2i+1,Bit_Length_2i+1\n";
	const unsigned long long TIME_LIMIT = 13549;
        const unsigned long long LIMIT = 500; // i such that 2*i+1 <= 1e8
	auto init_start = std::chrono::system_clock::now();
	for (unsigned long long i = 2; i <= LIMIT; ++i) {
	    auto now = std::chrono::system_clock::now();
	    auto elapsed = std::chrono::duration_cast<std::chrono::seconds>(now - init_start).count();
	    if (elapsed >= TIME_LIMIT) {
	        std::cout << "Time limit of " << TIME_LIMIT << " seconds reached. Stopping loop." << std::endl;
	        break;
	    }

	    auto start = std::chrono::system_clock::now();
	    ZZ one(i);
	    if (!is_prime_deterministic(one)) continue;
	    if (prime_set.find(i) == prime_set.end()) continue;

	    unsigned long long candidate = 2 * i + 1;
	    ZZ two(candidate);
	    if (!is_prime_deterministic(two)) continue;
	    if (prime_set.find(candidate) == prime_set.end()) continue;

	    auto end = std::chrono::system_clock::now();

	    std::time_t start_ts = std::chrono::system_clock::to_time_t(start);
	    std::time_t end_ts = std::chrono::system_clock::to_time_t(end);

	    outfile << std::put_time(std::localtime(&start_ts), "%Y-%m-%d %H:%M:%S") << ","
	            << std::put_time(std::localtime(&end_ts), "%Y-%m-%d %H:%M:%S") << ","
	            << i << "," << NumBits(ZZ(i)) << ","
	            << candidate << "," << NumBits(ZZ(candidate)) << "\n";

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
