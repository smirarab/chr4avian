### INSTALL: 

~~~bash
tar xvfz Triplet_rooting.tar.gz
gunzip resolved-genetrees.tre.gz
pip install TripletVoting
unzip examined-clades.zip
~~~

Note: to test and debug, you may want to use `head.tre` instead of the full `resolved-genetrees.tre` used by default. Each script has a way to give it a tree file; alternatively tick it with a symlink


### QQS

Text files:

* `newannotation.txt`: mapping species names to clade names
* `rev-annotation-nopar`: opposite mapping of the previous file
* `allgroups.txt`: The name of all clades one may encounter
* Files inside `examined-clades.zip` give definition of various clades. These are inputs to the next scripts. 


Scripts:

* `make-rec-group.sh`: Given three clade files (text files with one taxon per row) named `C1` `C2` and `S`, it generates the three resolved quadripartitions and an unresolved one in the following order. Each input clade file is given as a text file with one taxon per row (from the previous zip file). Note that the fourth part of the quadripartition (`O`) is implicit and doesn't need to be given. It's all leaves not in the first three. 
	* `(C1,C2)|(S,O)`
	* `(C1,S)|(C2,O)`
	* `(C2,S)|(C1,O)`
	* `(C1,2,S,O)`
* `score-rec-clade.sh`: uses tqDist through a [python script](Triplet_rooting/tq_distance_to_refs.py) to compute quartet distance to a quadripartition defined by three clades given as input: C1 C2 S, same as previous script
	* Output gives the number of quartets not matching each of the quadripartition topologies (in the given order above)
* `score-rec-clade-runs.sh`: Runs of the previous script; 
	* it also includes a one-liner to produce the final stat files used by the R scripts. 
	

### BQS 
To compute the BQS score for a group defined in a file called `examplegroup.txt`:

~~~bash
./score-clade.sh examplegroup.txt
~~~

To run on our predefined groups:

~~~bash
for x in `cat groups.txt`; do bash ./score-clade.sh $x; done
~~~

Script `summarize.sh` combines various files and creates the final `all-clade.stat` and other files used by the R script.
