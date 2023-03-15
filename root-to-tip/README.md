The script [root-to-tip.r](root-to-tip.r) draws the root to tip distance figures in the supplement. 

It uses the following data files:

* mono-mutime2.txt.xz: For each locus that finds Columbea as monophyletic, it gives us (in that order):
	* The locus id
	* The branch length above Columbea, 
	* The number of leaves from Columbea present. 
* mono-time.txt.xz: For every locus, it includes:
	*  The locus id
	*  An extra column with some (ignored) number if the Columbea is monophyletic. This ignored column is our unsuccesful attempt at dating each gene tree, an attempt that was not used in final results. 
* outliers-rooted.txt
	* The ID of outlier genes

Before these analyses are performed, all gene trees are rooted. The IDs above refer to the the IDs of rooted gene trees because 8 gene trees that could not be rooted are excluded from analyses (See below for the list). 
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
