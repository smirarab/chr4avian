#!/bin/bash
#SBATCH --job-name="clades"
#SBATCH --output="astral.%j.%N.out"
#SBATCH --partition=shared
#SBATCH --export=ALL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem=2G
#SBATCH -t 24:00:00
#SBATCH -A uot138
#SBATCH --export=ALL

# module load cpu/0.15.4  gcc/10.2.0 python/3.8.5

# Inputs: [a clade file like dove.txt] [input gene trees (default: resolved-genetrees.tre)]
# requires newick utilities

# Note: cat main.unresolved.tre alternative.unresolved.tre > references.tre
# coco.txt dove.txt flamingo.txt grebe.txt mesite.txt sandgrouse.txt

SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#SCRIPT_DIR=/home/smirarab/projects/avian/avian-data/ASTRAL/scripts/clade-analysis/

set -x
set -e

t=${1/.txt/}

gt=resolved-genetrees.tre
p=""
if [ $# == 2 ]; then
	gt=$2
	p="-"$gt
	nw_stats $gt|grep leaves|cut -f2 > leaves$p.txt
fi
cat $1|$SCRIPT_DIR/make-group.sh > ref-$t.tre

mkdir -p `dirname clade-$t$p.stat`
mkdir -p `dirname maxclade-$t$p.txt`

paste <( nw_prune -v $gt $( cat $1)|nw_stats - |grep leaves|cut -f 2 ) leaves$p.txt |awk '{print $1*($1-1)*($2-$1)*($2-$1-1)/4}' > maxclade-$t$p.txt 

python $SCRIPT_DIR//Triplet_rooting/tq_distance_to_refs.py -i $gt -r ref-$t.tre -m quartet_raw -o clade-$t$p.txt
paste maxclade-$t$p.txt clade-$t$p.txt |awk '{print $2,$3,$1,$4,($1-$4)/($1==0?1:$1)}' > clade-$t$p.stat
#f=clade-$t.txt; l=$(( $( cat $f |wc -l ) / 2 )); paste <( head -n $l $f ) <( tail -n $l $f )|awk '{print $2, $3, $6, $3-$6}' > clade-$t.stat

#zip -m textfiles.zip ref-$t.tre maxclade-$t$p.txt clade-$t$p.txt
