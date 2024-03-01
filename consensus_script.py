import pandas as pd
import numpy as np
import os

def consensus_prediction_fnc(path,cutoff):
    data=pd.read_csv(path,sep="\t",header=0)

    data_interest=data[data['consensus_percentile_rank']<cutoff]
    if data_interest.empty:
      return data_interest
    all_alleles=np.unique(data_interest['start'].values)
    new_dict={'allele':[],'number of alleles':[],'peptide':[],'start':[],'end':[],'consensus_percentile_rank':[],'ann_ic50':[],'ann_rank':[],'smm_ic50':[],'smm_rank':[],'minimum consensus score':[]}
    names=['allele','number of alleles','peptide','start','end','consensus_percentile_rank','ann_ic50','ann_rank','smm_ic50','smm_rank','minimum consensus score']
    for unique_start in all_alleles:
        #columns all - 'allele','seq_num','start','end','length','peptide', 'consensus_percentile_rank', 'ann_ic50', 'ann_rank', 'smm_ic50','smm_rank', 'comblib_sidney2008_score', 'comblib_sidney2008_rank'
        data_allele=data_interest[data_interest['start']==unique_start]
        
        alleles=data_allele['allele'].values.tolist()
        number_alleles=len(alleles)
        peptide=data_allele['peptide'].values.tolist()[0]
        start=data_allele['start'].values.tolist()[0]
        end=data_allele['end'].values.tolist()[0]
        cons_score=data_allele['consensus_percentile_rank'].values.tolist()
        min_cons=min(data_allele['consensus_percentile_rank'])
        ann_ic50=data_allele['ann_ic50'].values.tolist()
        ann_rank=data_allele['ann_rank'].values.tolist()
        smm_ic50=data_allele['smm_ic50'].values.tolist()
        smm_rank=data_allele['smm_rank'].values.tolist()
        vals=[alleles,number_alleles,peptide,start,end,cons_score,ann_ic50,ann_rank,smm_ic50,smm_rank,min_cons]

        for index,name in enumerate(names):
            new_dict[name].append(vals[index])

    df=pd.DataFrame.from_dict(new_dict)
    df=df.sort_values(by=['number of alleles','minimum consensus score'],ascending=[False,True])
    return df


