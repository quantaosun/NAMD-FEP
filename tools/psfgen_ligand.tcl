# Minimal local preparation template.
# Supply a validated ligand topology/parameter stream before using this file.
if {$argc != 4} {
    puts stderr "usage: vmd -dispdev text -e psfgen_ligand.tcl -- input.pdb output.psf output.pdb ligand.str"
    exit 2
}

set input_pdb [lindex $argv 0]
set output_psf [lindex $argv 1]
set output_pdb [lindex $argv 2]
set ligand_stream [lindex $argv 3]

package require psfgen
topology $ligand_stream
segment LIG { pdb $input_pdb }
coordpdb $input_pdb LIG
guesscoord
writepsf $output_psf
writepdb $output_pdb
