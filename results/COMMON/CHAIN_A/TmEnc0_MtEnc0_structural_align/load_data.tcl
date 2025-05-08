# Chargement des fichiers PDB
mol new MtEnc0_monomer2.pdb
mol new TmEnc0_monomer1.pdb

# Configuration des représentations graphiques

# Pour MtEnc_monomer2 (mol 0)
mol modstyle 0 0 NewCartoon
mol modcolor 0 0 ColorID 1 ; # Rouge
mol selupdate 0 0 on
mol smoothrep 0 0 0

# Pour TmEnc0_monomer1 (mol 1)
mol modstyle 0 1 NewCartoon
mol modcolor 0 1 ColorID 4 ; # Bleu
mol selupdate 0 1 on
mol smoothrep 0 1 0

# Ajustement de la vue
display resetview
scale by 1.2
rotate x by -30
rotate y by 15

