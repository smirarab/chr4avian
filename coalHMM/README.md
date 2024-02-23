This directory contains all necessary scripts for running an automated CoalHMM pipeline to compute the proportion of incomplete lineage sorting (ILS):

* `workflow.py` is a gwf (https://gwf.app) workflow, which will:
    * Run `autocoalhmm` for all possible trees for *Mesitornis unicolor*, *Caloenas nicobarica*, *Crotophaga sulcirostris* and *Phoenicopterus ruber*, with *Gallus gallus* as an outgroup, for all chromosomes. 
    	* If another user would like to run the gwf workflow, the `autocoalhmm` GitHub repo should be downloaded in the current directory. More information on `autocoalhmm` can be found here: Iker Rivas-González (2022), rivasiker/autocoalhmm: v1.0.0 (v1.0.0), Zenodo, https://doi.org/10.5281/zenodo.7277715.
    * Compute the proportion of sites that belong to each of the hidden states of CoalHMM in non-overlapping 100 kb windows, using the `python` script in `./scripts/get_tables_asym.py`. 

* `./results/11_chr4_column_bin` contains the tables outputted by `./scripts/get_tables_asym.py`. 

* `./analyses/ILS_analysis.Rmd` is an R Markdown used for plotting the horizon charts in Fig. 2. It makes use of the tables in `./results/11_chr4_column_bin`.

* `./data` contains metadata about the bird dataset, which is necessary to run the workflow. If another user would like to run the gwf workflow, the alignments should be in this directory, separated by chromosome and gzipped. 

Note that because of the large size of the raw `autocoalhmm` output files (the posterior decoding information per site for the HMM), the files are not included in this repo. 