#!/bin/bash

echo "Starting Forward FEP (A→B)"

# Run forward FEP windows sequentially (0 to 14)
for i in $(seq 0 14); do
    echo "========================================"
    echo "Running forward window $i"
    echo "========================================"
    namd2 +p32 fep_forward_${i}.conf | tee fep_forward_${i}.out
    echo "Window $i completed"
    echo ""
done

echo "Forward FEP completed!"
echo "FEP output files: fep_forward_*.fepout"
