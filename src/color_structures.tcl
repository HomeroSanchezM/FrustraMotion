
# Auto-generated VMD Tcl script

mol new "/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results/FRUSTRATION_MTENC/MTENC_CAPSIDS/FRUSTRATION_monomer_for_a_frame/MtEnc0_monomers/MtEnc0_monomer5.pdb" type pdb
mol new "/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results/FRUSTRATION_TMENC/TMENC_CAPSIDS/FRUSTRATION_monomer_for_a_frame/TmEnc0_monomers/TmEnc0_monomer26.pdb" type pdb

# Liste de couleurs à appliquer
set colors {grey grey yellow green yellow grey grey tan grey yellow green grey grey grey grey grey grey grey green grey tan grey grey grey grey yellow grey yellow grey grey green grey grey grey grey green green grey green grey grey red grey grey grey tan grey grey grey tan tan grey grey yellow yellow tan yellow green tan grey grey grey grey tan yellow tan yellow grey green yellow grey grey yellow red green green grey green grey tan grey green grey green tan yellow tan red green tan grey tan yellow tan grey grey tan grey yellow grey yellow yellow yellow green grey grey grey yellow yellow grey green grey yellow tan grey red tan yellow green green grey grey yellow grey tan grey tan green yellow grey green yellow grey grey grey grey grey grey yellow green yellow tan yellow grey tan grey tan yellow green yellow grey green green grey grey yellow tan yellow grey tan grey tan grey yellow tan tan grey grey grey yellow green green tan grey grey grey tan yellow grey grey green grey tan grey grey grey tan tan yellow yellow tan grey yellow yellow tan yellow yellow green green tan grey grey green green grey grey grey grey green yellow grey grey green green green grey grey green grey grey yellow green yellow green yellow green grey grey red green grey green grey grey grey grey grey grey grey tan tan green grey green grey green grey red grey green grey green tan green grey grey grey tan grey yellow green yellow green grey}

# Fonction pour appliquer la couleur à une molécule
proc apply_colors {molid colors} {
    set sel_all [atomselect $molid "name CA"]
    set resid_list [$sel_all get resid]
    set resid_unique [lsort -unique $resid_list]

    for {set i 0} {$i < [llength $colors]} {incr i} {
        set color [lindex $colors $i]
        set resid [lindex $resid_unique $i]
        if {$resid eq ""} { continue }
        set color_id [colorinfo colorid rgb $color]
        set sel [atomselect $molid "resid $resid"]
        $sel set colorID $color_id
        $sel delete
    }

    mol modcolor 0 $molid ColorID
    mol modstyle 0 $molid NewCartoon
}

# Appliquer aux deux molécules
apply_colors 0 $colors
apply_colors 1 $colors
