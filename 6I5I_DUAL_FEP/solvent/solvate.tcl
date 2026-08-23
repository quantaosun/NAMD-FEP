package require solvate
package require autoionize
mol delete all
mol load psf ligand.psf pdb ligand.pdb
solvate ligand.psf ligand.pdb -t 15 -o solvated
autoionize -psf solvated.psf -pdb solvated.pdb -neutralize -o ionized
mol delete all
mol load psf ionized.psf pdb ionized.pdb
set all [atomselect top all]
puts "CELL [molinfo top get {a b c}]"
puts "ORIGIN [molinfo top get {center}]"
exit
