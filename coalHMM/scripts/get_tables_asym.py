# This script can be used to get the ILS estimates for a specified chromosome, 
# in 100 kb windows. It outputs two different
# tables a table containing the ILS estimates from binary tags. 

import pandas as pd
import os.path
import sys
import numpy as np

chrom_name = sys.argv[2]
species = sys.argv[1]


def count_V(state, V_type):
    return int(list(state).count(V_type))
count_V0 = lambda x: count_V(x, 'V0')
count_V1 = lambda x: count_V(x, 'V1')
count_V2 = lambda x: count_V(x, 'V2')
count_V3 = lambda x: count_V(x, 'V3')


with pd.HDFStore('../../results/11_chr4_columb/{}/chr_{}/final_table.HDF'.format(species, chrom_name),  mode='r') as store:
    df = pd.DataFrame()
    for i, reread in enumerate(store.select(
        'Gallus_gallus_chr'+chrom_name, 
        chunksize = 1000000, 
        columns = ['Gallus_gallus','V0', 'V1', 'V2', 'V3']
        )):
        # Bin the human coordinates of the region in 100 kb windows
        reread['bin'] = (reread['Gallus_gallus']//100000)*100000
        # Get rid of NaN human coordinates
        run_df = reread[reread['Gallus_gallus'] != -1]
        # Back in the unbinned table, calculate maximum state per position
        run_df['state_V'] = run_df[['V0','V1','V2','V3']].idxmax(axis=1)
        run_df = run_df[['bin', 'state_V']]
        # Calculate ILS and count per bin
        new = run_df.groupby('bin').agg(
            {
                'state_V':[count_V0, count_V1, count_V2, count_V3],
            }
        ).reset_index()
        new.columns = ["_".join(x) for x in new.columns.ravel()]
        df = df.append(new)
    
df = df.reset_index()

wm = lambda x: np.average(x, weights=df.loc[x.index, "count"])

df = df.reset_index().groupby('bin_').agg(
    {
        'state_V_<lambda_0>':'sum',
        'state_V_<lambda_1>':'sum',
        'state_V_<lambda_2>':'sum',
        'state_V_<lambda_3>':'sum'}
    ).reset_index().rename(columns={
        'bin_':'position',
        'state_V_<lambda_0>':'V0',
        'state_V_<lambda_1>':'V1',
        'state_V_<lambda_2>':'V2',
        'state_V_<lambda_3>':'V3'
    })

df.to_csv('../../results/11_chr4_columb_bin/info_asym/{}_chr{}.csv'.format(species, chrom_name), index=False)
