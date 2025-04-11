#!/bin/bash

NUM_JOBS=20

MR_SOURCE="miller_rabin_with_check.cpp"

g++ -o executable "$MR_SOURCE" -lgmp -lgmpxx -O2 -pthread

rm -rf data
mkdir data
cd data
mkdir logs
mkdir csvs
cd ../

run_mr_test() {
    K=$1
    LOG_FILE="data/logs/miller_rabin_k${K}.log"
    OUTPUT_FILE="data/csvs/primality_conflicts_k${K}.csv"
    PRIME_FILE="primes1_till_1e8.txt"
    echo "Running Miller-Rabin test for K=$K..." | tee "$LOG_FILE"
    ./executable "$PRIME_FILE" "$OUTPUT_FILE" "$K" &>> "$LOG_FILE"
    echo "Completed K=$K" | tee -a "$LOG_FILE"
}

export -f run_mr_test

parallel -j $NUM_JOBS run_mr_test ::: {1..20}

echo "All tests completed."
