mol delete all
mol load psf ionized.psf pdb ionized.pdb
puts "CELL [molinfo top get {a b c}]"
puts "ORIGIN [molinfo top get {center}]"
exit
