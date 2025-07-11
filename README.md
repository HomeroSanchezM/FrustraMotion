# FrustraMotion
FrustraMotion is a tool for the study of local energetic frustration of the dynamics states of proteins.
Is composed of 3 main scripts: 
-Frustration.py
-Parser.py / Parser_contact.py
-visualisation.py / visualisation_contact.py 
 
## Use of scrips

### Parser.py :
This script takes as a parameter the directory containing the frustration data and creates various dataframes with the structure: frame × residues
for the single residue frustration data. 
It take as parameter a directory containing the frustration data obtain by Frustration.py or directly by using FrsutratometeR 
Usage: 

Parse the frustration data for protein complex, were directories have as code name `ProteinName_FrameNumber.done` 
```bash
python3 Parser.py path/to/directory/
```
This generates a directory `single_residue_dataframes/Not_isolated/ProteinName/` and store one CSV file for each chain of 
the protein complex with the frustration index (FI) of each residue for all the frames.

Parse the frustration data for a protein complex were the chain have been frustrated separetely, and directories have as code name `ProteinName_FrameNumber_ChainName.done`
```bash
python3 Parser.py path/to/directory/ --isolate
```
This generate a directory `sigle_residue_dataframe/Isolated/ProteinName/` and store one CSV file for each chain of the protein complex.

Parse the frustration data for a mono-chain protein, were directories have as code name `ProteinName_FrameNumber.done`
```bash
python3 Parser.py path/to/directory/ --true_isolate
```
This generate a directory `single_residue_dataframe/True_isolated/ProteinName` and store a single CSV file with the frustration index (FI) of each residue for all the frames.

### Visualisation.py :
This script takes as a parameter the directory containing the dataframes 



