# FrustraMotion

## Use of scrips

### parser_pdb.py :

This script take one or two PDB file and realise a global alignment to found the similar domains
Usage: 

get the sequence of one PDB 
```bash
python3 parser_pdb.py path/to/file.pdb
```
Get the sequence for two PBD and compare it with a global alignment
```bash
python3 parser_pdb.py path/to/file1.pdb path/to/file2.pdb
```

### structural_align.py :

This scripts take one PDB File containing multiples monomers (chains) and make a structural alignment of the monomers and mesure the mean RMSF for all the monomers per residue.
A plot is added to the plots directory, with name of type `rmsf_with_std_per_res_<TmEnc|MtEnc>_monomer_<(t)>.png`
The aligned PDB files are added to the result directory.
```bash
python3 structural_align.py path/to/file.pdb
```
If 2 files given as parameter, the script will align the 2 files structures entirely (not considering monmers)

Usage:
```bash
python3 structural_align.py path/to/file1.pdb path/to/file2.pdb
```

### structural_align_multiples_files.py

This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a structural alignment of this monomer and mesure the  mean RMSF for this monomers in each file (frames).
A plot is added to the plot directory, with name of type `rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png`
The aligned PDB files are added to the result directory.
```bash
python3 structural_align_multiples_files.py path/to/directory/ number_of_the_wanted_monomer
```
### frustration_plots.py

This scripts take one PDB File containing multiples monomers (chains) and make a frustration study of the monomers using FrustratometeR 
the results files are added to `results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_monomer_for_a_frame`. 
A plot is added to the `plots/frustration directory`, with name of type
`frustration_per_res_<TmEnc|MmEnc>_monomer_<(i)>.png`

```bash
python3 frustration_plots.py path/to/file.pdb
```

If 2 file given, frustration will be calculated for the monomers of the 2 files and a scatter-plot of the frustration means for each residue will be made
```bash
python3 frustration_plots.py path/to/file1.pdb path/to/file2.pdb 
```
### frustration_plots_multiples_files.py

This scripts take all the PDB Files of a given directory containing multiples monomers (chains)
extract a wanted monomer, and make a frustration study of the monomers in each file (frames) using FrustratometeR 
the results files are added to `results/FRUSTRATION_<MTENC|TMENC>/<MTENC|TMENC>_CAPSIDS/FRUSTRATION_frames_for_a_monomer`. 
A plot is added to the `plots/frustration directory`, with name of `type rmsf_with_std_per_frame_<TmEnc|MmEnc>_monomer_<(i)>.png`.

```bash
python3 frustration_plots_multiples_files.py path/to/directory/ number_of_the_wanted_monomer
```

If 2 file given, frustration will be calculated for the monomers of the 2 files and a scatter-plot of the frustration means for each residue will be made
```bash
python3 frustration_plots_multiples_files.py path/to/file1.pdb path/to/file2.pdb 
```