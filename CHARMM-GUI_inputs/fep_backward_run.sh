#!/bin/bash

echo "Starting Backward FEP (B→A)"

# Run backward FEP windows sequentially (15 down to 0)
for i in $(seq 15 -1 0); do
    echo "========================================"
    echo "Running backward window $i"
    echo "========================================"
    namd2 +p32 fep_backward_${i}.conf | tee fep_backward_${i}.out
    echo "Window $i completed"
    echo ""
done

echo "Backward FEP completed!"
echo "FEP output files: fep_backward_*.fepout"
