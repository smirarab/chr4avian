#!/bin/bash

if [ $# -lt 2 ]; then
	echo "USAGE: $0 [treefile] [tab-deliminated-mapping-file] [outgropus (optional; provided the clade names; best to privde more than one)]"
fi

tree=$1
shift
mapping=$1
shift

if [ $# -lt 1 ]; then
	outgroup='Palaeognathae'
else
	outgroup=$*
fi

nw_rename $tree $mapping |python3 -W ignore -c '
import treeswift as ts;
import sys; 
tl=ts.read_tree_newick(sys.stdin.read());
tl = tl if isinstance(tl, list) else [tl]
for t in tl:
	#sys.stderr.write(str(t))
	t.reroot(t.mrca("'$outgroup'".split(",")),branch_support=True);t.is_rooted=True;
	print(t.newick()[5:].strip());' |nw_condense -
