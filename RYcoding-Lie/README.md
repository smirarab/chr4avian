* `ry.R`: This is the R script that draws the results. 

* `rycoding-qqs.txt`: Shows the QQS results for ry coding. Settings are similar to `clade-rec.stat.xz` inside [genetreesupport](../genetreesupport/README.md) directory. See that README file. 
	* Note that the **first and second columns here have wrong indices** and should be ignored. 
	* Instead, we use `outlier-indices.txt`, `non-outlier-indices.txt`, and `control-indices.txt` to map indices of outliers region genes, non-outlier genes that happen to be on chr4, and control genes outside chr4 to indices of the 64K gene tree file. There are 1000 "true" control gene trees. 
	* The names of outlier genes are in `outliers.txt` 

* `400loci-order.txt`: the coordinates of the 400 genes. 
* `monophyletics.txt`: the file showing which loci have Columbea as monophyletic
