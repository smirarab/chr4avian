from gwf import Workflow
import os

gwf = Workflow()

sp_lst = ['GCF_000002315', 'GCA_009769525', 'GCA_901699155', 'GCA_009769465', 'GCA_017976375', 'GCF_003957565', 'GCA_009819775']

ref = 'GCA_017639555'

for chr in ['CM030196.1']:
    gwf.target('hal2maf_CICMAG_{}'.format(chr),
                inputs=['vgp_birds_60way.hal'],
                outputs=['vgp_birds_cicmag_{}.maf'.format(chr)],
                cores=4,
                memory='16g',
                walltime= '6:00:00',
                account='Primategenomes') << """
        hal2maf vgp_birds_60way.hal vgp_birds_cicmag_{}.maf --refGenome GCA_017639555.1_bCicMag1.pri --refSequence {} --targetGenomes GCA_009769525.1_bPteGut1.pri,GCA_901699155.2_bStrTur1.1,GCA_009769465.1_bTauEry1.pri,GCA_017976375.1_bCucCan1.pri,GCF_000002315.5_GalGal6,GCF_003957565.2_bTaeGut1.4.pri,GCA_009819775.1_bPhoRub2.pri --noDupes
        """.format(chr, chr)

if not os.path.exists('CICMAG/'):
    os.makedirs('CICMAG/')

for chr in ['CM030196.1']:
    if not os.path.exists('CICMAG/{}/'.format(chr)):
        os.makedirs('CICMAG/{}/'.format(chr))
    for sp in sp_lst:
        if not os.path.exists('CICMAG/{}/{}'.format(chr, sp)):
            os.makedirs('CICMAG/{}/{}'.format(chr, sp))
        gwf.target('Maffilter_new_CICMAG_{}_{}'.format(chr, sp),
                inputs=['vgp_birds_cicmag_{}.maf'.format(chr)],
                outputs=['CICMAG/{}/{}/filtered_new.maf'.format(chr, sp)],
                cores=8,
                memory='16g',
                walltime= '1:00:00',
                account='Primategenomes') << """
        cd ./filter
        ./maffilter param=./control_file_options_new_CICMAG CHR={} SP1={} SP2={} > /dev/null
        """.format(chr, ref, sp)
        gwf.target('maf2synteny_CICMAG_{}_{}'.format(chr, sp),
                        inputs=['CICMAG/{}/{}/filtered_new.maf'.format(chr, sp)],
                        outputs=['CICMAG/{}/{}/5000'.format(chr, sp)],
                        cores=4,
                        memory='16g',
                        walltime= '01:00:00',
                        account='Primategenomes') << """
                cd ./CICMAG/{}/{}/
                maf2synteny -o ./ filtered_new.maf -b 5000
                """.format(chr, sp)
        gwf.target('convert_to_csv_CICMAG_{}_{}'.format(chr, sp),
                        inputs=['CICMAG/{}/{}/5000'.format(chr, sp)],
                        outputs=['CICMAG/{}/{}/5000/block_info.csv'.format(chr, sp)],
                        cores=1,
                        memory='2g',
                        walltime= '01:00:00',
                        account='Primategenomes') << """
                python get_info.py {} {} {}
                """.format(chr, 'CICMAG', sp)
