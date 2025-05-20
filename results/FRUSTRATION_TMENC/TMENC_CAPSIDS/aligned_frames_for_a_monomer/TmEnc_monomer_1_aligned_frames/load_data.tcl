# Load 491 monomers TmEnc0_monomer1.pdb ... TmEnc4900_monomer1.pdb

set num_monomers 491
set base_filename "../results/FRUSTRATION_TMENC/TMENC_CAPSIDS/aligned_frames_for_a_monomer/TmEnc_monomer_1_aligned_frames/TmEnc"
set monomer_suffix "_monomer1.pdb"

# Colors - using VMD's default color IDs
set available_colors [list 0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]

# Calculate frame numbers (0, 10, 20, ..., 4900)
set frame_numbers {}
for {set i 0} {$i < $num_monomers} {incr i} {
    lappend frame_numbers [expr {$i * 10}]
}

# Load each structure with proper naming
set color_index 0
foreach frame $frame_numbers {
    set filename "${base_filename}${frame}${monomer_suffix}"
    puts "Loading $filename..."
    
    # Load PDB file
    mol new $filename type pdb waitfor all
    
    # Apply color (cycle through available colors)
    mol modcolor 0 [molinfo top] "ColorID" [lindex $available_colors $color_index]
    set color_index [expr {($color_index + 1) % [llength $available_colors]}]
    
    # Apply representation (optional)
    mol modstyle 0 [molinfo top] "NewCartoon"
}

puts "Successfully loaded $num_monomers structures."

