#!/usr/bin/env python3


import warnings
warnings.simplefilter(action='ignore')
import pandas as pd
import numpy as np
import os
import glob2
import re
import progressbar
import time
import matplotlib.pyplot as plt
import seaborn as sns

path= "/content/drive/MyDrive/GBC/s_filtered/" #Importing all csv files of samples

## path of annotated files ----IMP---
files= glob2.glob(os.path.join(path, "*.csv"))
ce_files = []
abs_files =[]

## removes exome and germline samples ---IMP---
for i in files:
  if '-SE8-' in i or '-SSE-' in i or 'SSE' in i or 'SE8' in i or '-B-' in i:
    abs_files.append(i)
  else:
    ce_files.append(i)
    
    
## filtration
for i in progressbar.progressbar(range(len(ce_files))):
  df = pd.read_csv(ce_files[i]).drop_duplicates('IGV_link')
  sp1=re.split(r'/|\\',ce_files[i])[-1]
  col_n=sp1.split('_')[0]
  df[f'{col_n}:DP'] = df[f'{col_n}:DP'].replace('.', '0').astype('float')
  df = df[df['FILTER']=='PASS']

  ## Pop frequency
  columns_to_process = ['esp6500siv2_all', 'ExAC_ALL', 'ExAC_SAS', 'AF', 'AF_sas', '1000g2015aug_all', '1000g2015aug_all', '1000g2015aug_SAS', 'Ref_Depth']

  for column in columns_to_process:
    df[column] = df[column].replace('.', '0').astype('float')

  df = df[df['esp6500siv2_all']<= 0.01]
  df = df[df['ExAC_ALL']<= 0.01]
  df = df[df['ExAC_SAS']<= 0.01]
  df = df[df['AF']<=0.01]
  df = df[df['AF_sas']<= 0.01]
  df = df[df['1000g2015aug_all']<= 0.01]
  df = df[df['1000g2015aug_all']<= 0.01]
  df = df[df['1000g2015aug_SAS']<= 0.01]
  df = df[df['Ref_Depth']>=2]
  df = df[df['ExonicFunc.ensGene']!='synonymous SNV']
  df = df[~df['InterVar_automated'].str.contains('benign', case=False)]
  df.to_csv(f'/content/drive/MyDrive/GBC/s_filtered/{col_n}_filtered.csv', index=False)