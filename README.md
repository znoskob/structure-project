# structure-project
### Introduction
This project contains scripts to run an automated protocol for searching the [CoSSMos](http://cossmos.slu.edu) database for 3D structural patterns among secondary structures using the Python programming language.

### Files
The following scripts are included, in the `scripts` folder:
* **clip.py**: description coming soon!
* **check_files.py**: description coming soon!
* **Avg_Structure.py**: description coming soon!
* **RMSD.py**: description coming soon!
* **RMSD_Tree_Folders.py**: description coming soon!
* **Avg_Structure-cluster.py**: description coming soon!

### Installation
1. Ensure you have [Python v2.7](https://www.python.org/downloads/) (Python 3 is not supported) and [Biopython](https://www.github.com/biopython/biopython).
2. Clone or download this repository.
3. Copy `matrix.py` from the this repository's `biopython_mods` folder and replace `Bio/motifs/matrix.py` in your Biopython working directory.  This modification adds a 'U' to the degenerate dictionary, which is required for RMSD_Tree_Folders.py to run, but does not affect any other functionality of the original Biopython file.

### Running the Scripts
Refer to instructions.pdf (coming soon!) for instructions on how to run the scripts.

### Collaboration
Pull requests are always welcome!

### License
See the [LICENSE](LICENSE.md) file for licensing information. 
