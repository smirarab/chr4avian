1. Add "Birds" prefix to ancestor names
~~~bash
wget https://raw.githubusercontent.com/ComparativeGenomicsToolkit/cactus/v2.0.4/src/cactus/cactus_progressive_config.xml
sed -i cactus_progressive_config.xml -e "s/default_internal_node_prefix=\"Anc\"/default_internal_node_prefix=\"BirdsAnc\"/g"
~~~

2. Make the WDL, which is uploaded then run directly on Terra.
~~~bash
cactus-prepare --wdl ./ucsc_birds_scaled_consensus.txt --noLocalInputs --preproc
essBatchSize 10 --alignDisk 3000G --halAppendDisk 3000G --preprocessDisk 3000G --defaultDisk 1000G --defaultCores 64 --gpu --gpuCount 8 --defaultMemory 385G --configFile ./cactus_progressive_config.xml
~~~

3. To make the assembly hub (note we're still using tools as part of Cactus v2.0.4):
Chicken was renamed galGal6 to be compatible with UCSC Genome Browser
~~~bash
printf "GCF_000002315.5_GalGal6\tgalGal6" > rename_chicken.tsv
halRenameGenomes birds.hal rename_chicken.tsv
~~~

4. Maf was constructed with:
~~~bash
for seq in $(halStats birds.hal --sequenceStats galGal6 | awk '{print $1}') ; do echo ${seq::-1}; done | paralell -j 16 "hal2maf birds.hal birds_{}.maf --refGenome galGal6 --refSequence {} --onlyOrthologs --inMemory --noAncestors"
~~~

5. And the hub with:
~~~
hal2assemblyHub.py ./js birds.hal  vg_birds_hub
~~~

6. And the MAF was added manually as a track to galGal6 in the hub afterwards.

Included:
* 57-alignment.xlsx: The list of included species
* ucsc_birds_scaled_consensus.txt: The input guide tree
* ucsc_birds_scaled_consensus.wdl: output .wdl for "cactus-prepare".
