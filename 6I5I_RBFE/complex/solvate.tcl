package require solvate
package require autoionize
mol delete all
mol load psf complex.psf pdb complex.pdb
solvate complex.psf complex.pdb -t 15 -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -neutralize -o ionized
set all [atomselect top all]
puts "BOX_MINMAX [measure minmax $all]"
puts "BOX_CENTER [measure center $all]"
exit
