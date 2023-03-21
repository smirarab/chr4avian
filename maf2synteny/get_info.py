#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 22 12:38:12 2022

@author: iker
"""

import pandas as pd
import sys

chr = sys.argv[1]
ref = sys.argv[2]
sp = sys.argv[3]

for res in [5000]:
    with open('./{}/{}/{}/{}/blocks_coords.txt'.format(ref, chr, sp, res), 'r') as file:
        arr = []
        for line in file.readlines():
            if 'Seq_id' in line:
                names = line.strip().split('\t')
            elif '-----' in line:
                df = pd.DataFrame(data = arr, columns = names)
                df.to_csv('./{}/{}/{}/{}/seq_info.csv'.format(ref, chr,sp, res), index = False)
                break
            else:
                arr.append(line.strip().split('\t'))


for res in [5000]:
    with open('./{}/{}/{}/{}/blocks_coords.txt'.format(ref, chr, sp, res), 'r') as file:
        arr = []
        acc = False
        for line in file.readlines():
            if '-----' in line:
                acc = True
                continue
            if not acc:
                continue
            elif 'Block' in line:
                block = line.strip().split(' ')[1].replace('#', '')
            elif 'Seq_id' in line:
                continue
            else:
                arr.append(line.strip().split('\t')+[block])
    names = ['Seq_id', 'Strand', 'Start',   'End',     'Length', 'Block']
    df = pd.DataFrame(data = arr, columns = names)
    df.to_csv('./{}/{}/{}/{}/block_info.csv'.format(ref, chr, sp, res), index = False)
