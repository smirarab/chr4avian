INSTALL: 

~~~bash
tar xvfz Triplet_rooting.tar.gz
gunzip resolved-genetrees.tre.gz
pip install TripletVoting
unzip examined-clades.zip
~~~

To compute the quartet score for a group define in a file called `examplegroup.txt`:
~~~bash
./score-clade.sh examplegroup.txt
~~~

To run on our predefined groups:
~~~bash
for x in `cat groups.txt`; do bash ./score-clade.sh $x; done
~~~

To do monophyletic analyses, run:
~~~bash
./count-mopnophyletic.sh resolved-genetrees.tre `cat groups.txt|tr '\n' ' '`
~~~

Examples of clade definition that can be used later to define new groups: `clade-defenitions.zip`

Script `summarize.sh` combines various files and creates the final `all-clade.stat` and other files used by the R script.
