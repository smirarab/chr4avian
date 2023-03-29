#!/bin/bash

# arg 1 : gene tree file
# arg 2, 3, ...: name of files including clades to be tested  

g=$1
shift

while (( "$#" )); do
	echo processing $1 ....
	c=`cat $1|tr '\n' ' '`
	while read l; do
		a=$( echo $l | nw_clade -m - $c 2>/dev/null|nw_stats -|grep leaves )
		echo $a
	done < $g > monphyletic-$1
	grep -n -o  -f  $1 $g |sed -e "s/:/\t/g"|cut -f1| uniq -c > present-$1
	echo writtten in  present-$1 and  monphyletic-$1
	shift 
done
