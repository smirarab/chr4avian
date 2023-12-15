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

# Inputs: [a clade file A like dove.txt] [reciprocal clade of A] [sister of A and and its reciprocal] [input gene trees (default: resolved-genetrees.tre)]
# requires newick utilities




SCRIPT_DIR=$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )
#SCRIPT_DIR=/home/smirarab/projects/avian/avian-data/ASTRAL/scripts/clade-analysis/

set -x
set -e

t=${1/.txt/}_${2/.txt/}

gt=resolved-genetrees.tre
#gt=head.tre
p=""
if [ $# == 4 ]; then
	gt=$4
	p="-"$gt
	nw_stats $gt|grep leaves|cut -f2 > leaves$p.txt
fi
$SCRIPT_DIR/make-rec-group.sh $1 $2 $3 > ref-rec-$t.tre

mkdir -p `dirname recclade-$t$p.stat`
mkdir -p `dirname maxrecclade-$t$p.txt`


python $SCRIPT_DIR//Triplet_rooting/tq_distance_to_refs.py -i $gt -r ref-rec-$t.tre -m quartet_raw -o clade-rec-$t$p.txt
