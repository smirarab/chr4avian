* `ry.R`: This is the R script that draws the results. 

* `rycoding-qqs.txt`: Shows the QQS results for ry coding. Settings are similar to `clade-rec.stat.xz` inside [genetreesupport](../genetreesupport/README.md) directory. See that README file. 
	* Note that the **first and second columns here have wrong indices** and should be ignored. 
	* Instead, we use `outlier-indices.txt`, `non-outlier-indices.txt`, and `control-indices.txt` to map indices of outliers region genes, non-outlier genes that happen to be on chr4, and control genes outside chr4 to indices of the 64K gene tree file. There are 1000 "true" control gene trees. These files are created by `figureoutproblems.sh`
	* The names of outlier genes are in `outliers.txt` 

* The actual gene trees are given in the following files:
	* With contraction:

	~~~
	control.original.alrt0.95.gene.trees.gz    
	control.ry_coding.alrt0.95.gene.trees.gz                 
	non-outlier_chr4.original.alrt0.95.gene.trees.gz                
	non-outlier_chr4.ry_coding.gene.alrt0.95.gene.trees.gz 
	outlier_chr4.ry_coding.alrt0.95.gene.trees.gz
	outlier_chr4.original.alrt0.95.gene.trees.gz 
	~~~
	
	* Without contraction:
	
	~~~
	control.original.gene.trees.gz 
	control.ry_coding.gene.trees.gz 
	non-outlier_chr4.original.gene.trees.gz 
	non-outlier_chr4.ry_coding.gene.trees.gz 
	outlier_chr4.original.gene.trees.gz 
	outlier_chr4.ry_coding.gene.trees.gz
	~~~
	
* Ignore: `outliers.chr4.txt`