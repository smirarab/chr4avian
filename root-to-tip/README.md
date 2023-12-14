Before these analyses are performed, all gene trees are rooted.
To root gene trees, we used the following command to root each fully resolved gene tree at the LCA of Palaeognathae.

~~~ bash
 grep -v treefile 63430_named_random_region.gene.trees.filter_mRNA| nw_reroot - `grep Palaeognathae annotations.txt|cut -f1|tr '\n' ' '`

21952 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
32668 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
42529 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
43658 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
45524 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
54762 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
60281 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
61138 ERROR: Outgroup's LCA is tree's root - cannot reroot. Try -l.
~~~ 
Here, 

* `63430_named_random_region.gene.trees.filter_mRNA` is the set of gene trees without contracting low support branches. This file is available on an ftp site linked from the main paper repository: https://github.com/smirarab/avian-data/blob/master/ASTRAL/genetrees.txt
* `annotations.txt` is provided in the file [annotations.txt](annotations.txt)
* The output of the script above is saved in a file called `63430_random_region.gene.trees.filter_mRNA_rooted.tre.gz` (too big for github). 

### Main figure

The script [root-to-tip.r](root-to-tip.r) draws the root to tip distance figures in the supplement. 

It uses the following data files:

* outliers.txt: Indicates which genes are in the outlier regions. Gene IDs are for all genes not the rooted ones. 
* mono-mutime2.txt.xz: For each locus that finds Columbea as monophyletic, it gives us the following (in that order). For others, it gives NA. Gene ID is for the rooted ones. 
	* The locus id
	* The branch length above Columbea, 
	* The number of leaves from Columbea present. 

