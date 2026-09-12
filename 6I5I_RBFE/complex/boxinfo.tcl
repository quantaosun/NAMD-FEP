mol delete all
mol load psf ionized.psf pdb ionized.pdb
set all [atomselect top all]
set mm [measure minmax $all]
set c [measure center $all]
puts "BOXMIN $mm"
puts "BOXCEN $c"
exit
