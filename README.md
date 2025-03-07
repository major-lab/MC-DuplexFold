# MC-DuplexFold
# MC-DuplexFold (mcdf) v 1.1
# Copyright 2025 Simon Chasles <simon.chasles@umontreal.ca>
# For basic compilation:   gcc -o mcdf -fopenmp mcduplexfold.c -lm
# For basic execution:     ./mcdf -s {first RNA sequence} -t {second RNA sequence}
# Sirt-1:miR34a example:   ./mcdf -s ACACCCAGCUAGGACCAUUACUGCCA -t UGGCAGUGUCUUAGCUGGUUGU -p 1
# Sirt-1:miR34a example:   ./mcdf -s ACACCCAGCUAGGACCAUUACUGCCA -t UGGCAGUGUCUUAGCUGGUUGU -f 5
# HNF4a:mir34a example:    ./mcdf -s aacauggccuaagggccacaucccacugcca -t UGGCAGUGUCUUAGCUGGUUGU -p 1
# For help:                ./mcdf -h
# Note: mcff must be executable in the current directory.
