#!/bin/bash

# Step 1: Go into Prime_Generation_and_RSA, run make clean and then make
echo "🔧 Cleaning and building C++ backend..."
cd Prime_Generation_and_RSA || { echo "Directory Prime_Generation_and_RSA not found!"; exit 1; }

make clean
make || { echo "❌ Make failed in Prime_Generation_and_RSA"; exit 1; }

# Step 2: Compile check_prime_with_aks.cpp from Prime_Validate
echo "🧪 Compiling AKS checker..."
cd ../Prime_Validate || { echo "Directory Prime_Validate not found!"; exit 1; }

g++ check_prime_with_aks.cpp -o check_prime_with_aks -lntl -lgmp || { echo "❌ Compilation of check_prime_with_aks.cpp failed!"; exit 1; }

# Step 3: Compile the Java file
echo "☕ Compiling Java applet..."
cd ..  # Back to base directory

javac PrimeGeneratorApplet.java || { echo "❌ Java compilation failed!"; exit 1; }

# Step 4: Launch the applet using appletviewer
echo "📟 Running applet using appletviewer..."
appletviewer PrimeGeneratorApplet.html || { echo "❌ Failed to launch appletviewer"; exit 1; }

echo "✅ All steps completed!"

