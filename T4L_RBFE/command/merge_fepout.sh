#!/bin/bash
#Script to merge output files of .fepout
#including complex_forward.fepout, complex_backward.fepout, solvent_forward.fepout, and solvent_backward.fepout, total 4 files.

snum=1
enum=16

# 1. Merge complex forward
echo "Merging complex_forward.fepout..."
> complex_forward.fepout   # clear or create empty file
while [[ $snum -le $enum ]]; do
    if [ -f "md_forward_${snum}.fepout" ]; then
        cat md_forward_${snum}.fepout >> complex_forward.fepout
    else
        echo "  WARNING: md_forward_${snum}.fepout not found, skipping."
    fi
    ((snum++))
done
echo "  Done. Lines: $(wc -l < complex_forward.fepout)"

# 2. Merge complex backward
snum=1
echo "Merging complex_backward.fepout..."
> complex_backward.fepout
while [[ $snum -le $enum ]]; do
    if [ -f "md_backward_${snum}.fepout" ]; then
        cat md_backward_${snum}.fepout >> complex_backward.fepout
    else
        echo "  WARNING: md_backward_${snum}.fepout not found, skipping."
    fi
    ((snum++))
done
echo "  Done. Lines: $(wc -l < complex_backward.fepout)"

# 3. Merge solvent forward
snum=1
echo "Merging solvent_forward.fepout..."
> solvent_forward.fepout
while [[ $snum -le $enum ]]; do
    if [ -f "md_forward_${snum}.fepout" ]; then
        cat md_forward_${snum}.fepout >> solvent_forward.fepout
    else
        echo "  WARNING: md_forward_${snum}.fepout not found, skipping."
    fi
    ((snum++))
done
echo "  Done. Lines: $(wc -l < solvent_forward.fepout)"

# 4. Merge solvent backward
snum=1
echo "Merging solvent_backward.fepout..."
> solvent_backward.fepout
while [[ $snum -le $enum ]]; do
    if [ -f "md_backward_${snum}.fepout" ]; then
        cat md_backward_${snum}.fepout >> solvent_backward.fepout
    else
        echo "  WARNING: md_backward_${snum}.fepout not found, skipping."
    fi
    ((snum++))
done
echo "  Done. Lines: $(wc -l < solvent_backward.fepout)"

echo ""
echo "All 4 .fepout files merged:"
echo "  complex_forward.fepout"
echo "  complex_backward.fepout"
echo "  solvent_forward.fepout"
echo "  solvent_backward.fepout"
echo ""
echo "Now run ParseFEP in VMD on each file separately, then compute:"
echo "  ΔΔG = ΔG_complex - ΔG_solvent"
