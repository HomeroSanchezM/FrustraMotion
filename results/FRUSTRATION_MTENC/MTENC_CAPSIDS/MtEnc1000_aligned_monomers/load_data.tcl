# Load 60 monomers MtEnc1000_monomer1.pdb ... MtEnc1000_monomer60.pdb

set num_monomers 60
set base_filename "MtEnc1000_monomer"

# Colors
set available_colors [list 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]

for {set i 1} {$i <= $num_monomers} {incr i} {
    set filename "${base_filename}${i}.pdb"
    puts "Loading $filename..."

    # load PDB
    mol new $filename type pdb waitfor all

    # Apply colors
    set molID [molinfo top]
    mol modstyle 0 $molID NewCartoon
    mol modcolor 0 $molID ColorID [expr {($i - 1) % [llength $available_colors]}]

    # Transparence (0.0 = opaque, 1.0 = invisible)
    mol material Transparent
}

