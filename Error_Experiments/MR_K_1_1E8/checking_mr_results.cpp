#include <iostream>
#include <fstream>
#include <unordered_set>
#include <sstream>
#include <string>

using namespace std;

// Function to load prime numbers from the text file into an unordered set
unordered_set<int> load_prime_numbers(const string &filename) {
    unordered_set<int> prime_set;
    ifstream file(filename);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << filename << endl;
        exit(EXIT_FAILURE);
    }
    int num;
    while (file >> num) {
        prime_set.insert(num);
    }
    file.close();
    return prime_set;
}

int main() {
    string prime_file = "primes1_till_1e8.txt";
    string csv_file = "primality_test_results_k1.csv";
    unordered_set<int> prime_set = load_prime_numbers(prime_file);
    
    ifstream file(csv_file);
    if (!file.is_open()) {
        cerr << "Error: Could not open " << csv_file << endl;
        return EXIT_FAILURE;
    }
    
    string line;
    getline(file, line); // Skip header line
    
    int incorrect_count = 0;
    int line_count = 0;
    
    while (getline(file, line)) {
        line_count++;
        if (line_count % 100000 == 0) {
            cout << "Processed " << line_count << " lines. Incorrect classifications so far: " << incorrect_count << endl;
        }
        
        stringstream ss(line);
        string timestamp, num_str, time_taken, primality;
        getline(ss, timestamp, ',');
        getline(ss, num_str, ',');
        getline(ss, time_taken, ',');
        getline(ss, primality, ',');
        
        int number = stoi(num_str);
        
        if (primality == "Prime" && prime_set.find(number) == prime_set.end()) {
            incorrect_count++;
        }
    }
    
    file.close();
    cout << "Final incorrect prime classifications: " << incorrect_count << endl;
    return 0;
}
