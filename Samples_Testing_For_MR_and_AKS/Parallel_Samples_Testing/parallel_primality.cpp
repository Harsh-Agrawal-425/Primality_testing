#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <chrono>
#include <gmp.h>
#include <NTL/ZZ.h>
#include <sstream>
#include <string>
#include <cmath>
#include <stdio.h>
#include <time.h>
#include <mpi.h>

#include <math.h>
#include <NTL/ZZ.h>
#include <NTL/RR.h>
#include <NTL/ZZ_pE.h>
#include <NTL/ZZ_pEX.h>
#include <ctime>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>
#include <csignal>

using namespace std;
using namespace NTL;

ZZ generate_random_prime(int bits) {
    return GenPrime_ZZ(bits); 
}

void power_mod_gmp(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus)
{
    mpz_powm(result, base, exponent, modulus);
}

bool miller_rabin_test_gmp(const mpz_t n, int k)
{
    if (mpz_cmp_ui(n, 2) == 0 || mpz_cmp_ui(n, 3) == 0)
    {
        return true;
    }

    if (mpz_cmp_ui(n, 2) < 0 || mpz_even_p(n))
    {
        return false;
    }

    mpz_t d, n_minus_1, a, x;
    mpz_init(d);
    mpz_init(n_minus_1);
    mpz_init(a);
    mpz_init(x);

    mpz_sub_ui(n_minus_1, n, 1);

    mpz_set(d, n_minus_1);
    unsigned long r = 0;

    while (mpz_even_p(d))
    {
        mpz_tdiv_q_2exp(d, d, 1);
        r++;
    }

    gmp_randstate_t state;
    gmp_randinit_mt(state);
    gmp_randseed_ui(state, std::chrono::system_clock::now().time_since_epoch().count());

    for (int i = 0; i < k; i++)
    {
        if (mpz_cmp_ui(n, 4) <= 0)
        {
            mpz_clear(d);
            mpz_clear(n_minus_1);
            mpz_clear(a);
            mpz_clear(x);
            gmp_randclear(state);
            return true;
        }

        mpz_sub_ui(a, n, 3);
        mpz_urandomm(a, state, a);
        mpz_add_ui(a, a, 2);

        power_mod_gmp(x, a, d, n);
        if (mpz_cmp_ui(x, 1) == 0 || mpz_cmp(x, n_minus_1) == 0)
        {
            continue;
        }

        bool witness = false;
        for (unsigned long j = 0; j < r - 1; j++)
        {
            mpz_powm_ui(x, x, 2, n);

            if (mpz_cmp_ui(x, 1) == 0)
            {
                mpz_clear(d);
                mpz_clear(n_minus_1);
                mpz_clear(a);
                mpz_clear(x);
                gmp_randclear(state);
                return false;
            }

            if (mpz_cmp(x, n_minus_1) == 0)
            {
                witness = true;
                break;
            }
        }

        if (!witness)
        {
            mpz_clear(d);
            mpz_clear(n_minus_1);
            mpz_clear(a);
            mpz_clear(x);
            gmp_randclear(state);
            return false;
        }
    }

    mpz_clear(d);
    mpz_clear(n_minus_1);
    mpz_clear(a);
    mpz_clear(x);
    gmp_randclear(state);

    return true;
}

int step0(ZZ n);
int step1(ZZ n);
long step2(ZZ n);
int step3(ZZ n, long r);
int step4(ZZ n, long r);
int step5(ZZ n, long r);
int step6();
ZZ  gcd(ZZ m, ZZ n );
ZZ  order(ZZ  r, ZZ a);
ZZ phi(ZZ n);

int count_bits(const string &number_str)
{
    ZZ num(INIT_VAL, number_str.c_str());
    return NumBits(num);
}

int is_prime_deterministic(ZZ n){
    int k;  
    long r; 
    
    k = step0(n);
    if(k){
        if(k == 1){
            return 0;
        }
        return 1;
    }
    
    k = step1(n);
    if(k){
        return 0;
    }

    r = step2(n);

    k = step3(n, r);
    if(k){
        return 0;
    }

    k = step4(n, r);
    if(k){
        return 1;
    }

    k = step5(n, r);
    if(k){
        return 0;
    }

    k = step6();
    if(k){
        return 1;
    }
    
    return 0;
}

/////////////////////////////////////////////////////////
//
// Extra functions.
//
// These include the Euclidean algorithm, the order of
// an element, and Euler's totient function.
//
/////////////////////////////////////////////////////////

/*
 * Euclidean algorithm, should be self explanitory.
 */

ZZ gcd( ZZ m, long n ){
    ZZ k, z;
    z = 0;
    z = n+z; 
    if (z<m) {
        swap(z,m);
    }
    
    while (m % z != 0) {
        k = m % z;
        m = z;
        z = k;
    }
    return z;
}

/*
 * The order of a modulo r is the smallest integer
 * k such that a^k = 1 mod r.  It's denoted o_r(a).
 */

ZZ order(ZZ a, long r){

    ZZ k,o,z;
    o = 1;
    z = 1;
    k = a;
    z = r + z;

    while(k!=1){
        k*=a; k%=z; o++;
    }

    return o;
}

/*
 * Phi is Euler's totient function. This is not optimized
 * at all, and does not work if the input has a prime 
 * factor larger than the largest prime in the NTL
 * PrimeSeq library, which isn't usually a problem.
 */


ZZ phi(long x){

    ZZ n, ph,p; ph = 1;
    n = 0; n = n + x; 

    int k;

    for(PrimeSeq s; n>1; ){
        p = s.next();
        k = 0;
        while(n%p == 0){
            k++; n/=p;
        }
        if(k>0)
        ph *= power(p,k-1)*(p-1);
    }
    return ph;
}

/////////////////////////////////////////////////////////
//
// Steps.
//
// These are the steps of the algorithm as outlined
// in the paper.  A few of them could have been
// put in the main function due to their simplicity,
// but they were kept as steps to maintain consistency
// with the paper instead.
//
/////////////////////////////////////////////////////////

/*  Step 0:
 *  It is not difficult to test the first primes
 *  less than 2000 to see if they are divisors of
 *  n.  When we do this, step1 becomes much faster
 *  as we only need to check if n = a^b for a>2000,
 *  since if n = a^b for some a, then there exists
 *  a prime p <= a that divides n.
 *
 *  This code snippet was taken from Victor Shoup's
 *  website (the author of NTL).  It can be found
 *  here:
 *
 *  http://www.shoup.net/ntl/doc/tour-ex1.html
 *
 *  This function is not a part of the AKS algorithm.
 *
 *  Returns 1 if composite, 2 if prime, 0 if it's
 *  not sure.
 *
 */

int step0(ZZ n){

    long p;

    PrimeSeq s;
    
    p = s.next();  
    while (p && p < 2000) {
        if ((n % p) == 0){
            return (1 + (n == p)); 
        } 
        p = s.next();              
    }
    return 0; 
}

/* Step 1:
 * If n=a^b for a in N and b>1, output composite.
 * The idea here is to test different values of b.
 * We start with b=2, and take the b-th root of n.
 * We call this number A, which is in the class RR.
 * If, when we round A down, we get the same thing,
 * we can test if A^b == n.  (We do this second step
 * to deal with rounding errors.  Just because when we
 * round down we get 0, our precision may not be perfect
 * so we test the exponetiation for these few cases.)  
 *
 * If A^b=n, our test is over, else we increase b by one.
 *
 * We keep doing this until b reaches log_2000 (n).
 * Since we tested divisibility of all primes less
 * than 2000, we can stop here.
 *
 * Returns 1 if composite, 0 if we don't know.
 */

int step1(ZZ n){

    RR A, B, N, log_2000_n;
    long b = 2;

    log_2000_n = log(n)/log(2000);
    N = to_RR(n);

    while(b<log_2000_n){

        B = 1.0/b;    
        A = pow(N,B); 
        if((A - to_RR(FloorToZZ(A)))==0.0){ 
            if(power(FloorToZZ(A),b) == n)  
                return 1;
        }
        b++;
    }

    return 0; 

}

/* Step 2:
 * This is a function that searches for the 
 * smallest r such that ord_r(n) > log_2(n)^2.
 */

long step2(ZZ n){

    RR log2n, l;
    ZZ k;
    long r;

    log2n = NumBits(n);
    l = power(log2n,2);
    k = FloorToZZ(l);
    
    for(r=2;;r++){
        if(gcd(n,r)==1){
            if(order(n,r) > k){
                return r; 
            }
        }
    }
}

/* Step 3:
 * This computes gcd(a,n) for values less than
 * r, where r is the value computed in step 3.
 * Notice that since we have already determined that
 * for any value k<2000, k does not divide n, we
 * do not need to check values of gcd(a,n) for values
 * a<2000.
 *
 * Returns 1 if composite, 0 if we don't know.
 */

int step3(ZZ n, long r){

    long a;

    for(a=2002;++a<r;){
        if ( (gcd(n,a)%n)>1 ) 
            return 1; 
    }
    
    return 0;
}

/* Step 4:
 * This is an extremely simple function and
 * is only laid out here to be consistent with the 
 * paper.  If n<=r, we know that n is prime.
 */

int step4(ZZ n, long r){
    if(n <= r)
        return 1;
    return 0;
}

/* Step 5:
 * This is the heart of the algorithm.  It tests
 * for different values of a, whether or not
 * (x+a)^n = x^n + a, over the polynomial ring
 * F = Z_n[x] / x^r - 1.  If for any value of a the
 * equality doesn't hold, we output composite.  
 * Otherwise we continue.
 *
 * Unfortunately the input of polynomials using NTL is
 * not straight forward. Comments should provide insight
 * as to what is being defined and why.
 *
 * Inspiration for step 5 is taken from George Poulose's
 * AKS implementation.  Nothing was copied directly but
 * certain methods were immitated.  In particular, the 
 * use of polynomials, and their input.
 *
 * His program may be found here:
 *
 * http://www.gpoulose.com/gc/AKS_cpp.txt
 *
 * Function returns 1 if composite, 0 if we don't know.
 * Although, for all intents and purposes, this is the 
 * last step, so seeing a 0 means we have a prime.  But,
 * in keeping with the outline of the algorithm, we 
 * pretend that we don't know what happens when we see
 * a 0.
 */

int step5(ZZ n, long r){

    ZZ l;
    long a;

    NTL::ZZ_p::init(n); 
                       

    ZZ_pX polymod(r,1); 
    polymod -=1 ;       

    ZZ_pXModulus mod(polymod); 

    ZZ_pX RHS(1, 1);
    PowerMod(RHS, RHS, n, mod);
    
    l = FloorToZZ(sqrt(to_RR(phi(r)))*NumBits(n));

    for(a=1; a<=l; a++){
        ZZ_pX LHS(1,1); 
        LHS += a;       
        PowerMod(LHS,LHS,n,mod); 
        LHS -= a;       
        if(LHS != RHS)  
          return 1;     
    }

    return 0;
}

/* Step 6:
 * If we've got this far, we know that n is prime.
 * Again, this function is merely a placeholder
 * to be consistent with the paper.
 */

int step6(){
    return 1;
}

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    
    int bit_length = 34 + rank; 
    
    if (bit_length > 60) {
        MPI_Finalize();
        return 0;
    }

    string filename = "primality_" + to_string(bit_length) + "bits.csv";
    ofstream csv_file(filename);
    csv_file << "Bits,Number,AKS_Time(ms),MR_Time(ms),AKS_Result,MR_Result\n";

    int iterations = 10;

    for (int i = 0; i < iterations; i++) {
        ZZ number = generate_random_prime(bit_length);
        ostringstream oss;
        oss << number;
        string num_str = oss.str();
        cout << "Processor " << rank << " | Bit-Length: " << bit_length << " | Test " << i+1 << " | Number: " << num_str << endl;

        mpz_t n_gmp;
        mpz_init(n_gmp);
        mpz_set_str(n_gmp, num_str.c_str(), 10);

        auto start_aks = chrono::high_resolution_clock::now();
        bool is_prime_aks = is_prime_deterministic(number);
        auto end_aks = chrono::high_resolution_clock::now();
        double aks_time = chrono::duration<double, milli>(end_aks - start_aks).count();

        auto start_rm = chrono::high_resolution_clock::now();
        bool is_prime_rm = miller_rabin_test_gmp(n_gmp, 40);
        auto end_rm = chrono::high_resolution_clock::now();
        double rm_time = chrono::duration<double, milli>(end_rm - start_rm).count();

        csv_file << bit_length << "," << num_str << "," << aks_time << "," << rm_time << ","
                 << (is_prime_aks ? "Prime" : "Composite") << ","
                 << (is_prime_rm ? "Prime" : "Composite") << "\n";

        mpz_clear(n_gmp);
    }

    csv_file.close();
    cout << "Processor " << rank << " finished. Results saved to " << filename << endl;

    MPI_Finalize();
    return 0;
}
