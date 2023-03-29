INSTALL: 

~~~bash
tar xvfz Triplet_rooting.tar.gz
gunzip resolved-genetrees.tre.gz
pip install TripletVoting
~~~

Run for a particular clade (e.g., `Accipitriformes.txt`):
~~~bash
./score-clade.sh clade-defs/Accipitriformes.txt
~~~

Examples of clade definition files are available under: `clade-defenitions.zip`

Script `summarize.sh` combines various files and creates the final `all-clade.stat` file.
