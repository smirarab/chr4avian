from gwf import Workflow
import os
import pandas as pd
import sys

gwf = Workflow()

names = pd.read_csv('./data/count_per_species.txt', sep='\t')
dist_root = pd.read_csv('./data/distance_to_root.csv')

full_chr_lst = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2', '20', '21', '22', '23', '24', '25', '26', '27', '28', '3', '4', '5', '6', '7', '8', '9', 'LGE64', 'M', 'W', 'Z']

sp = ['CALNIC', 'MESUNI', 'PHORUB', 'CROSUL']

comb = []
for i in range(len(sp)):
    for j in range(i+1, len(sp)):
        for k in range(j+1, len(sp)):
            print([i, j, k])
            st_2 = [sp[i], sp[j], sp[k]]
            comb.append('_'.join(st_2))
            st_2 = [sp[i], sp[k], sp[j]]
            comb.append('_'.join(st_2))
            st_2 = [sp[j], sp[k], sp[i]]
            comb.append('_'.join(st_2))

for i in range(0, len(comb), 3):
    print("'{}_GALGAL', '{}_GALGAL', '{}_GALGAL',".format(comb[i], comb[i+1], comb[i+2]))

# If the directory does not exist
if not os.path.isdir('./results/11_chr4_columb/'):
    # Create directory
    os.mkdir('./results/11_chr4_columb/')

for i in comb:
    spe = ''
    for sp in i.split('_'):
        spe += ' '+names[names['code'] == sp]['taxon'].item()
    spe += ' Gallus_gallus'
    spe = spe[1:]
    sp_sep = spe.split(' ')[0:3]
    three_sp = dist_root[dist_root['species'].isin(sp_sep)]
    three_sp = three_sp.sort_values('tab', ascending=False)
    species = (spe, ' '.join(list(three_sp['species'])[0:2]))
    # If the directory does not exist
    if not os.path.isdir('./results/11_chr4_columb/{}_GALGAL'.format(i)):
        # Create directory
        os.mkdir('./results/11_chr4_columb/{}_GALGAL'.format(i))
    for chrom in full_chr_lst:
        chrom_name = chrom
        # Define target
        gwf.target('{}_chr{}'.format(i+'_GALGAL', chrom_name),
                    inputs=['./data/maf_files_new/chr{}.maf.gz'.format(chrom), '../autocoalhmm/autocoalhmm.py'], 
                    outputs=['./results/11_chr4_columb/{}/chr_{}/final_table.HDF'.format(i+'_GALGAL', chrom)],
                    cores=1,
                    memory='2g',
                    walltime= '00:20:00',
                    account='Primategenomes') << """
            mkdir ./results/11_chr4_columb/{}/chr_{}
            cd ./results/11_chr4_columb/{}/chr_{}
            python ../../../../autocoalhmm/autocoalhmm.py {} Gallus_gallus.chr{} ../../../../data/maf_files_new/chr{}.maf.gz {}
            """.format(i+'_GALGAL', chrom, i+'_GALGAL', chrom, species[0], chrom_name, chrom, species[1])

# If the directory does not exist
if not os.path.isdir('./results/11_chr4_columb_bin/'):
    # Create directory
    os.mkdir('./results/11_chr4_columb_bin/')
# If the directory does not exist
if not os.path.isdir('./results/11_chr4_columb_bin/info_asym/'):
    # Create directory
    os.mkdir('./results/11_chr4_columb_bin/info_asym/')

# For each species
for sp_short in comb:
    sp_short = sp_short+'_GALGAL'
    for chrom in full_chr_lst:
        chrom_name = chrom
        if not os.path.isfile('./results/11_chr4_columb/{}/chr_{}/final_table.HDF'.format(sp_short, chrom)):
            continue
        # Define target
        gwf.target('ILS_{}_chr{}'.format(sp_short, chrom_name),
                    inputs=['./results/11_chr4_columb/{}/chr_{}/final_table.HDF'.format(sp_short, chrom)], 
                    outputs=['./results/11_chr4_columb_bin/info_asym/{}_chr{}.csv'.format(sp_short, chrom)],
                    cores=4,
                    memory='64g',
                    walltime= '04:00:00',
                    account='Primategenomes') << """
            python ./scripts/get_tables_asym.py {} {}
            """.format(sp_short, chrom)
