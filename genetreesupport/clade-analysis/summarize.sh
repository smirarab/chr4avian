#!/bin/bash
bash rename.sh
rm all-*stat.gz
#grep " " minus-*stat|sed -e "s/:/ /" -e "s/.stat//" > all-minus.stat
#grep " " clade-*stat|sed -e "s/:/ /" -e "s/.stat//" > all-clade.stat
grep " " clade-*stat|sed -e "s/:/ /" -e "s/.stat//" -e "s/clade-//g" > all-clade.stat
gzip all*stat
grep "" monphyletic-*.txt|sed -e "s/:/\t/" -e "s/://g" -e "s/#leaves//"| sed -e "s/monphyletic-//g" -e "s/.txt//" > all-monphyletic.txt
grep "" present-*.txt|sed -e "s/[:]/\t/g" -e "s/  */\t/g" -e "s/present-//g" -e "s/.txt//" > all-present.txt
