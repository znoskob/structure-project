# structure-project
### Introduction
This project contains scripts to run an automated protocol for searching the [CoSSMos](http://cossmos.slu.edu) database for 3D structural patterns among secondary structures using the Python programming language.

### Files
The following scripts are included, in the `scripts` folder:
* **clip.py**: This script uses the RNA CoSSMos search results file as input and downloads the corresponding .pdb files.  Then each result is clipped from the .pdb file and renumbered.  H atom coordinates are removed from the clipped .pdb files.
* **check_files.py**: This script performs a quality check on the clipped .pdb files. If clipped .pdb files are empty, missing atom coordinates, or contain multiple coordinate sets for the same atom, the files are flagged and moved to the "bad_clips" folder.
* **Avg_Structure.py**: This script calculates an average structure for each unique sequence (including closing base pairs) and then calculates an all-atom RMSD of each structure to the average structure.  The structure with the smallest RMSD is selected as the sequence-representative structure.
* **RMSD.py**: This script calculates all-against-all RMSD values for all of the sequence-representative structures. The RMSD is calculated using all atoms in the sugar-phosphate backbone and three atoms from each base (N9, C8, and C4 for purines or N1, C2, and C6 for pyrimidines).
* **RMSD_Tree_Folders.py**: This script calculates an UPGMA tree from the RMSD values of the sequence-representative structures.  The degenerate sequence is determined at each inner node on the tree.  The user is prompted to supply an RMSD cut-off and the script creates folders containing sequence-representative structures that have branch-length values within the RMSD cut-off.
* **Avg_Structure-cluster.py**: This script calculates an average structure for the sequence-representative structures in the selected folder and then calculates an RMSD of each structure to the average structure using all atoms from the backbone and three atoms from each base.  The structure with the smallest RMSD is selected as the cluster-representative structure.

Within the `scripts` folder, there are options to characterize on Linux or Windows operating systems. The follosing scripts are included in the `Characterization` folders:
* **dssr.sh**: This script runs DSSR on all PDB files in a given directory and exports all data to a merged .JSON file.
* **json_stacks_10_31_17.py**: This script parses the JSON file for stacking information and exports the data to a text file.
* **json_conf_pucker.py**: This script parses the JSON file for base conformation and sugar pucker information and exports the data to a text file.
* **json_hbond.py**: This script parses the JSON file for hydrogen bond information (non pairing) and exports the data to a text file.
* **json_bp.py**: This script parses the JSON file for base pairing information and exports the data to a text file.
* **json_detailedbp_5_26.py**: This script parses the JSON file for pairing hydrogen bond information and exports the data to a text file.
* **create.py**: This script converts all text files to .CSV files for easier analysis.
* **count_updated_10_31_17.py**: This script calculates the interactions common to a majority of representative structures in a cluster, and outputs the data as percentages in the terminal and exports the data to respective text files.

### Installation
1. Ensure you have [Python v2.7](https://www.python.org/downloads/), (Python 3 is not supported), [Biopython](https://www.github.com/biopython/biopython/), and [DSSR](http://forum.x3dna.org/site-announcements/download-instructions/).
2. Clone or download this repository.
3. Copy `matrix.py` from the this repository's `biopython_mods` folder and replace `Bio/motifs/matrix.py` in your Biopython working directory.  This modification adds a 'U' to the degenerate dictionary, which is required for RMSD_Tree_Folders.py to run, but does not affect any other functionality of the original Biopython file.

### Running the Scripts
Refer to instructions.pdf (coming soon!) for instructions on how to run the scripts.

### Collaboration
Pull requests are always welcome!

### License
See the [LICENSE](LICENSE.md) file for licensing information. 
