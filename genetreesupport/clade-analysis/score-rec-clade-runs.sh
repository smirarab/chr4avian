sbatch score-rec-clade.sh Accipitriformes.txt Strigiformes.txt CPBTL-Coli.txt
sbatch score-rec-clade.sh Aequornithes.txt Phaethontimorphae.txt Strisores.txt
sbatch score-rec-clade.sh CPBTL-Coli.txt StrigiformesAccipitriformes.txt Australaves.txt
sbatch score-rec-clade.sh CPBTL.txt Coliiformes.txt StrigiformesAccipitriformes.txt
sbatch score-rec-clade.sh Casuariiformes.txt Apterygiformes.txt TinamiformesRheiformes.txt
sbatch score-rec-clade.sh Charadriiformes.txt Gruiformes.txt Opisthocomiformes.txt
sbatch score-rec-clade.sh Columbiformes.txt OtherColumbimorphae.txt Otidimorphae.txt
sbatch score-rec-clade.sh Columbimorphae.txt OtherColumbiformes.txt Otidimorphae.txt
sbatch score-rec-clade.sh Columbimorphae.txt Phoenicopteriformes.txt Passera.txt
sbatch score-rec-clade.sh Cuculiformes.txt Otidiformes.txt Musophagiformes.txt
sbatch score-rec-clade.sh ElementavesTelluraves.txt Columbaves.txt  Phoenicopteriformes.txt
sbatch score-rec-clade.sh ElementavesTelluraves.txt Columbaves.txt Phoenicopteriformes.txt
sbatch score-rec-clade.sh GruiformesCharadriiformes.txt Opisthocomiformes.txt StrisoresAequornithesPhaethontimorphae.txt
sbatch score-rec-clade.sh Otidimorphae.txt Columbimorphae.txt ElementavesTelluraves.txt
sbatch score-rec-clade.sh PhaethontimorphaeAequornithes.txt Strisores.txt GruiformesCharadriiformesOpisthocomiformes.txt
sbatch score-rec-clade.sh StrisoresAequornithesPhaethontimorphae.txt GruiformesCharadriiformesOpisthocomiformes.txt Telluraves.txt
sbatch score-rec-clade.sh Telluraves.txt Elementaves.txt Columbaves.txt
sbatch score-rec-clade.sh Tinamiformes.txt Rheiformes.txt ApterygiformesCasuariiformes.txt
sbatch score-rec-clade.sh Tinamiformes.txt Rheiformes.txt OtherPalaeognathae.txt

grep  "" clade-rec-*txt|sed -e "s/clade-rec-//" -e "s/.txt:/ /" -e "s/_/ /"  > clade-rec.stat
