#!/bin/bash

echo "Starting backward FEP (A→B)"

# Run forward FEP windows sequentially
for i in {15..1}; do
    echo "========================================"
    echo "Running forward window $i"
    echo "========================================"
    namd2 +p32 fep_backward_${i}.conf | tee fep_backward_${i}.out
    echo "Window $i completed"
    echo ""
done

echo "Backward FEP completed!"
echo "FEP output files: fep_backward_*.fepout"
