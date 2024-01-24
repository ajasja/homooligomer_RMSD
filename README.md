# homooligomer_RMSD
Script to calculate RMSD between symmetric homo-oligmoers.
This is useful if you want to align a C10 homooligomer to C10 that has a different chain order. 
The scripts first aligns chain A in the reference structure to chain A in the mobile structure. 

# Instalation instructions
Biopython is needed for alignment and PDB parsing.
It can be installed using

`pip install biopython` or `conda install biopython`

# Example run
align_oligomers.py data/refc10.pdb data/mobilec10.pdb --out mobile_aligned.pdb

# Outputs
The aligned PDB is saved to a file and the RMSD is printed to the output (and saved into the PDB as a remark).
