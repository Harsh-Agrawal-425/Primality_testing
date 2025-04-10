#!/bin/bash

# Number of parallel jobs (adjust based on your system)
NUM_JOBS=20

# Paths to source files
MR_SOURCE="miller_rabin_with_check.cpp"
# CHECKER_SOURCE="checking_mr_results.cpp"

# Compile Miller-Rabin test program if not compiled
#if [ ! -f "miller_rabin_executable" ]; then
#    echo "Compiling Miller-Rabin test program..."
g++ -o executable "$MR_SOURCE" -lgmp -lgmpxx -O2 -pthread
#fi

## Compile the checking script if not compiled
#if [ ! -f "checking_executable" ]; then
#    echo "Compiling checking script..."
#    g++ -o checking_executable "$CHECKER_SOURCE" -O2
#fi

rm -rf data
mkdir data
cd data
mkdir logs
mkdir csvs
cd ../

# Function to run Miller-Rabin test
run_mr_test() {
    K=$1
#    DIR_NAME="MR_K_${K}_1E8"
    LOG_FILE="data/logs/miller_rabin_k${K}.log"
    OUTPUT_FILE="data/csvs/primality_conflicts_k${K}.csv"
#    CHECK_LOG="$DIR_NAME/checking_log.log"
    PRIME_FILE="primes1_till_1e8.txt"
    
    # Create directory
#    mkdir -p "$DIR_NAME"
    
    # Run the test and save the output
    echo "Running Miller-Rabin test for K=$K..." | tee "$LOG_FILE"
    ./executable "$PRIME_FILE" "$OUTPUT_FILE" "$K" &>> "$LOG_FILE"
    
    # Check results
#    echo "Checking results for K=$K..." | tee -a "$LOG_FILE"
#    ./checking_executable "$PRIME_FILE" "$OUTPUT_FILE" &>> "$CHECK_LOG"
    
    echo "Completed K=$K" | tee -a "$LOG_FILE"
}

export -f run_mr_test

# Run the tests in parallel for K=1 to 20
parallel -j $NUM_JOBS run_mr_test ::: {1..20}

echo "All tests completed."
