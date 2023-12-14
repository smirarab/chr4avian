
* [analyze.r](analyze.r) draws the supplementary figures used in the paper. It uses the following files. 
* The three trees used in the analyses are given in [t1t2t3.tre](t1t2t3.tre).
* [chr4-summary-ELB.tsv](chr4-summary-ELB.tsv): This file includes the results of the omega selection test. For each of the 992 genes, it gives:
	*  coordinates of the gene on the chromosome (`mRNA_begin` and `mRNA_end`). 
	* Ignore column `label`. It is incorrect. The script will add the correct label using the correct coordinates. 
	* Columns `Tree1_two_ratio_lnL`, `Tree2_two_ratio_lnL` ,`Tree3_two_ratio_lnL` given the log likelihood with two ratio model. 
	* The actual foreground omega for each case is given in `Tree1_two_ratio_forground`, `Tree2_two_ratio_forground`, `Tree3_two_ratio_forground` columns.
	* The p values for each of the three topologies: `lrt_1_1_p_value`, `lrt_2_2_p_value`, `lrt_3_3_p_value`
	* The script draws foreground omega and p-values associated with the tree with the highest of the `ratio_lnL` values, according to the columns above. 
	*  Other columns can be safely ignored. 
