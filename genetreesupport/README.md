This directory includes the scripts for the main analyses showing quartet score metrics (BQS and QQS) and monophyly, in addition to taxon sampling results

### Quartet scores (BQS and QQS)

* Check out [clade-analysis](./clade-analysis/) folder for scripts that perform the actual quartet computing using tqDist and other scripts we developed. 

* [draw-movingaverage.r](draw-movingaverage.r): The main R code for drawing moving average based on the text files described below. 
	* This file also does some of the calculations necessary to compute QQS and other measures (see below).  
	* It also calculates the **p-values**. (look for `P-value` and function `pvalues` in the file)
* `63K_trees.names_header.txt.xz`: includes only basic information about 63,430 loci, such as their start and end position. 
* `all-clade.stat.xz` and `'new-all-clade.stat.xz'`: info for branch quartet score (BQS) for each group
	* The first file includes some extra clades not included as part of the final analyses.  The second file includes some clades missing from our original analyses. The R script selects and renames the appropriate clades. 
	* The third column is the gene number.
	* The fourth column is the number of quartets relevant to that clade (`ref`)
	* The fifth column is the number of quartets not matching ot that gene tree (`difference`)
	* The sixth column (called `ratio` in the R script) gives the actual BQS and is `1-difference/ref`. 
* `clade-rec.stat.xz`: The QQS results.
	* First two columns give the names of the two children of the clade being tested. Let's call these `C1` and `C2`. Let the sister of `C1` and `C2` be `S` and everything outside 	the parent of this clade be `O`.
	* Third column: have values 1, 2, 3, or 4:
		* 1 means the dominant topology:    `C1,C2|S,O`
		* 2 means the alternative topology: `C1,S|C2,O`
		* 3 means the alternative topology: `C2,S|C1,O`
		* 4 means the unresolved topology:  `C1,C2,S,O`
	* Fourth column is gene numbers.
	* Fifth column is the number of quartets that are not found in the quadripartition corresponding to column 3. Let these numbers be $x_1$, $x_2$, $x_3$, or $x_4$, depending on the third column. These are computed using tqDist tool. 
	* To compute QQS, 
		* Let $d_1, d_2,$ and $d_3$ be the number of *resolved* quartets conflicting with each topology (those not found in the unresolved quadripartition (4) are not relevant to this branch). Note that $$d_1=x_1-x_4 \\ d_2=x_2-x_4 \\ d_3=x_3-x_4.$$  The R script computes these values (columns T1, T2, and T3). 
		* Let $a_1, a_2,$ and $a_3$ be the number of *resolved* quartets supporting each topology. Note that $$d_1=a_2+a_3 \\ d_2=a_2+a_1 \\ d_3=a_1+a_2$$ because each quartet conflicting with any topology has to support one of the other two topologies.  
		* Simple manipulation shows that the QQS, defined as $$\frac{a_1}{a_1+a_2+a_3}$$ (for first topology) is given by $$\frac{d_2+d_3-d_1}{d_1+d_2+d_3}$$ and ditto for the other two topologies. These are computed in the R code. 
		* These are values shown in moving averages and are used to compute the p-values. Note that to compute monophyly, we simply test to see if $a_i=1$.

### Taxon sampling

* `removed-count.tsv`: maps the name of files to how many of each group is removed
* `all.stat.xz`: For each subsampled tree (first column), it shows:
	* Gene number 
	* main: score of S2023 topology
	* alt: score of J2014 topology
	* diff: the difference between the two. 
* The diff column (which needs to be flipped in sign) is what is shown in the figure. 