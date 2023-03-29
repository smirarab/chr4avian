#!/bin/bash
set -e
set -x

g=$(cat - |xargs -I@ grep @ newannotation.txt |cut -f2|sort|uniq)
echo "(("$(echo $g|tr ' ' ',')"),("$( diff <( echo $g|tr ' ' '\n' ) <( cat allgroups.txt )|grep ">"|sed -e "s/> //g"|tr '\n' ',')"));"|sed -e "s/,)/)/g"|nw_rename - rev-annotation-nopar.txt 
