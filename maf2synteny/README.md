This directory contains the scripts for computing synteny blocks for a subset of the VGP bird species:
* `workflow.py` is a gwf (https://gwf.app) workflow, which will:
    * Extract the selected species from the HAL alignment into MAF format using `hal2maf`. 
    * For each pair of species containing a focal species and flamingo (used as a reference):
        * Extract the pairwise MAF alignment using `maffilter`.
        * Run `maf2synteny` on the pairwise alignment to compute the synteny blocks.
        * Extract the information from the output of `maf2synteny` using a custom python script (`get_info.py`)
* `./filter/control_file_options_new_CICMAG` is the control file for running `maffilter`. 
* `./CIGMAG` contains the output from the synteny computations. Each sub-directory contains the raw output from `maf2synteny` (`blocks_coords.txt`, `coverage_report.txt` and `genomes_permutations.txt`), and the post-processed output from `get_info.py` (the synteny block information in `block_info.csv`, and the sequence information in `seq_info.csv`).