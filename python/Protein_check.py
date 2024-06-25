#! /usr/bin/python
### Protein check input function
#Packages
import numpy as np
import pandas as pd

def Protein_check(df):
    seqs=df['Sequence'].values
    ind=0
    std = list("ACDEFGHIKLMNPQRSTVWY")
    for i in range(len(seqs[0])):
        if not (seqs[0][i] in std):
            ind = i
            break
    if ind!=0:
        string='Error : Unexpected character in protein sequence in position: '+str(ind)+' '+str(seq[0][ind])
        return string
    else:
        string='Correct format of the selected protein for the analysis. You can proceed to Parameters or Execute Tab panel'
        return string

        
