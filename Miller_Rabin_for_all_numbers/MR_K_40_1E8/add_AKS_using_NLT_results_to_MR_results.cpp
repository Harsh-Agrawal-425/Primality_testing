#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <cmath>
#include <NTL/ZZ.h>
#include <chrono>
#include <stdio.h>
#include <time.h>
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

int main(int argc, char *argv[])
{
    if (argc != 3)
    {
        cerr << "Usage: " << argv[0] << " <input_csv_file> <output_csv_file>" << endl;
        return 1;
    }

    string input_file = argv[1];
    string output_file = argv[2];

    ifstream infile(input_file);
    if (!infile.is_open())
    {
        cerr << "Failed to open input file: " << input_file << endl;
        return 1;
    }

    ofstream outfile(output_file);
    if (!outfile.is_open())
    {
        cerr << "Failed to open output file: " << output_file << endl;
        infile.close();
        return 1;
    }

    string line;
    auto total_start_time = chrono::high_resolution_clock::now();

    if (getline(infile, line))
    {
        outfile << line << ",Num_Bits,Determ_Output,AKS_Time" << endl;
    }

    long long total_numbers = 0, conflict_count = 0, prime_count_mr = 0, prime_count_aks = 0;
    
    while (getline(infile, line))
    {
        istringstream iss(line);
        string timestamp, number_str, time_taken, primality;

        if (getline(iss, timestamp, ',') && getline(iss, number_str, ',') &&
            getline(iss, time_taken, ',') && getline(iss, primality))
        {
            total_numbers++;
            int bits = count_bits(number_str);
            string determ_output;

            if (primality == "Composite")
            {
                determ_output = "Composite";
                outfile << line << "," << bits << "," << determ_output << ",NA" << endl;
            }
            else if (primality == "Prime")
            {
                prime_count_mr++;
                ZZ num(INIT_VAL, number_str.c_str());

                auto start_time = chrono::high_resolution_clock::now();
                int is_prime = is_prime_deterministic(num);
                auto end_time = chrono::high_resolution_clock::now();
                chrono::duration<double> elapsed = end_time - start_time;
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
                    cout << "Conflict found for number: " << number_str << endl;
                }

                outfile << line << "," << bits << "," << determ_output << "," << aks_time << endl;
            }
            else
            {
                determ_output = "Unknown";
                outfile << line << "," << bits << "," << determ_output << ",NA" << endl;
            }

            if (total_numbers % 1000 == 0)
            {
                cout << "Processed " << total_numbers << " numbers so far..." << endl;
            }
        }
    }

    auto total_end_time = chrono::high_resolution_clock::now();
    chrono::duration<double> total_elapsed = total_end_time - total_start_time;
    double total_runtime = total_elapsed.count();

    infile.close();
    outfile.close();

    cout << "Processing complete!" << endl;
    cout << "Total numbers processed: " << total_numbers << endl;
    cout << "Numbers declared Prime by Miller-Rabin: " << prime_count_mr << endl;
    cout << "Numbers confirmed Prime by AKS: " << prime_count_aks << endl;
    cout << "Conflicts found (MR=Prime, AKS=Composite): " << conflict_count
         << " out of " << total_numbers << " numbers ("
         << (double)conflict_count / total_numbers * 100 << "%)" << endl;
    cout << "Total execution time: " << total_runtime << " seconds" << endl;

    ofstream log_file(output_file, ios::app);
    if (log_file.is_open())
    {
        log_file << "Total Execution Time (s)," << total_runtime << endl;
        log_file.close();
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
