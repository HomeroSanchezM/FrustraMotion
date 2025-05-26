
    # FrustratometeR calculation
    library(frustratometeR)

    # Lire le fichier PDB
    pdb_file <- "/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results/FRUSTRATION_MTENC/MTENC_CAPSIDS/FRUSTRATION_monomer_for_a_frame/MtEnc0_monomers/MtEnc0_monomer26.pdb"
    output_dir <- "/home/homero/Documentos/M1/S2/Stage/FrustraMotion/results/FRUSTRATION_MTENC/MTENC_CAPSIDS/FRUSTRATION_monomer_for_a_frame/MtEnc0_frustration_seqdist_12_isolate"

    # Calculer la frustration
    results <- calculate_frustration(PdbFile = pdb_file, Mode = "singleresidue", SeqDist = 12, ResultsDir = output_dir, Graphics = FALSE)

    