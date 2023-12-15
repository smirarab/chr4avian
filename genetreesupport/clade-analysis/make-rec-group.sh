#!/bin/bash
set -e
set -x

# Given three clade files (text files with one taxon per row) it generates the the three resolved quadripartitions and an unresolved one

g1=$(cat $1 |xargs -I@ grep @ newannotation.txt |cut -f2|sort|uniq|tr '\n' ',')
g2=$(cat $2 |xargs -I@ grep @ newannotation.txt |cut -f2|sort|uniq|tr '\n' ',')
gs=$(cat $3 |xargs -I@ grep @ newannotation.txt |cut -f2|sort|uniq|tr '\n' ',')
gall=$(cat $1 $2 $3 |xargs -I@ grep @ newannotation.txt |cut -f2|sort|uniq)

echo "((("$g1"),("$g2")),("$gs"),("$( diff <( echo $gall|tr ' ' '\n' ) <( cat allgroups.txt )|grep ">"|sed -e "s/> //g"|tr '\n' ',')"));"|sed -e "s/,)/)/g"|nw_rename - rev-annotation-nopar.txt 
echo "((("$g1"),("$gs")),("$g2"),("$( diff <( echo $gall|tr ' ' '\n' ) <( cat allgroups.txt )|grep ">"|sed -e "s/> //g"|tr '\n' ',')"));"|sed -e "s/,)/)/g"|nw_rename - rev-annotation-nopar.txt 
echo "((("$g2"),("$gs")),("$g1"),("$( diff <( echo $gall|tr ' ' '\n' ) <( cat allgroups.txt )|grep ">"|sed -e "s/> //g"|tr '\n' ',')"));"|sed -e "s/,)/)/g"|nw_rename - rev-annotation-nopar.txt 
echo "(("$g2"),("$gs"),("$g1"),("$( diff <( echo $gall|tr ' ' '\n' ) <( cat allgroups.txt )|grep ">"|sed -e "s/> //g"|tr '\n' ',')"));"|sed -e "s/,)/)/g"|nw_rename - rev-annotation-nopar.txt 
