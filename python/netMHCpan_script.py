
###Packages
import os
import pandas as pd
import numpy as np


#### Functions
def netMHCpan_fnc(path,allel_map):
    ###########################
    dat1=netMHCpan_analysis(path,allel_map)
    return dat1


def netMHCpan_analysis(res_path,allel_map):
    
    with open(res_path,'r') as f:
        lin=f.readlines()
    df_a=pd.DataFrame()
    for j in range(len(allel_map)):
        inds=[i if allel_map[j] in lin[i] else -1 for i in range(len(lin))]
        inds=np.array(inds)
        b=inds>-1
        final_inds=inds[b]
        dat=[lin[i].split() for i in final_inds[:-1]]
        df=pd.DataFrame(dat)
        if df.shape[1]>17:
            df_1=df[df[17]=='SB']
            df_2=df[df[17]=='WB']
        else:
            continue
        if df_a.empty:
            df_a=pd.concat([df_1,df_2],ignore_index=True)  
        else:
            df_b=pd.concat([df_1,df_2],ignore_index=True)

            df_a=pd.concat([df_a,df_b],ignore_index=True)
    if df_a.empty:
        return df_a

    dict_dat={'Position':df_a[0].values.tolist(),'Alleles':df_a[1].values.tolist(),'Sequences':df_a[2].values.tolist(),'Protein':df_a[10].values.tolist(),'netMHCpan_Rank_Score_EL':df_a[12].values.tolist(),'netMHCpan_Rank_Score_BA':df_a[14].values.tolist(),'Affinity(nM)':df_a[15].values.tolist(),'Bind_level':df_a[17].values.tolist()}
    df_f=pd.DataFrame(dict_dat)
    df_f=df_f.sort_values(by='netMHCpan_Rank_Score_EL',ascending=False)

    proteins=np.unique(df_f['Protein'].values.tolist())
    for i in range(len(proteins)):
        df_i=df_f[df_f['Protein']==proteins[i]]
        pos_i=np.unique(df_i['Position'].values.tolist())
        for j in range(len(pos_i)):
            df_j=df_i[df_i['Position']==pos_i[j]]
            position_j=pos_i[j]
            alleles_j=df_j['Alleles'].values.tolist()
            sequence_j=np.unique(df_j['Sequences'].values.tolist())[0]
            protein_j=proteins[i]
            scores_j=df_j['netMHCpan_Rank_Score_EL'].values.tolist()
            scores_j=[np.float32(k) for k in scores_j]
            scores_jj=df_j['netMHCpan_Rank_Score_BA'].values.tolist()
            scores_jj=[np.float32(k) for k in scores_jj]
            affin_j=df_j['Affinity(nM)'].values.tolist()
            affin_j=[np.float32(k) for k in affin_j]
            bind_level_j=df_j['Bind_level'].values.tolist()
            sb_counts=bind_level_j.count('SB')
            wb_counts=bind_level_j.count('WB')
            temp_df=pd.DataFrame({'Position':position_j,'Alleles':[alleles_j],'Sequences':sequence_j,'Protein':protein_j,'netMHCpan_Rank_Scores_EL':[scores_j],'netMHCpan_Rank_Scores_BA':[scores_jj],'Affinities(nM)':[affin_j],'Bind-levels':[bind_level_j],'Strong_Binders_Counts':sb_counts,'Weak_Binders_Counts':wb_counts})
            if (i==0) and (j==0):
                unique_df=temp_df
            else:
                unique_df=pd.concat([unique_df,temp_df],ignore_index=True)
    unique_df=unique_df.sort_values(by=['Strong_Binders_Counts','Weak_Binders_Counts'],ascending=[False,False])
    #return [df_f,unique_df]
    return unique_df
