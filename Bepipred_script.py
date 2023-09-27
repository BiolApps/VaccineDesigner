
#Packages
import os
import subprocess
import numpy as np
import pandas as pd
import time
import re


##Functions
def Bepipred3(dataframe,path_bepipred,results_bepipred_dir='',thres_bepi=0.1512,base_neg=1,amino_region=10,secondary_score=0,bepipred_opt=1):
    path=os.getcwd()
    seq=dataframe['Sequence'][0]
    with open('temp.fasta','w+') as f:
        f.write('>seq1\n')
        f.write(seq)
    file_path=path+'/temp.fasta'
    #Bepipred
    os.chdir(path_bepipred)
    command = 'python3 bepipred3_CLI.py -i '+file_path+' -o '+ results_bepipred_dir +' -pred vt_pred'
    os.system(command)
    df = Bepipred3_analysis(results_bepipred_dir,'seq1',thres_bepi,base_neg,amino_region,secondary_score,bepipred_opt)
    df['Bepipred mean score']=np.float64(df['Bepipred mean score'].values)
    return df

def Bepipred3_fnc(fasta_file,code_path,result_path):
    os.chdir(code_path)
    command = 'python3 bepipred3_CLI.py -i '+fasta_file+' -o '+ result_path +' -pred vt_pred'
    process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,shell=True)
    stdout, stderr =process.communicate()
    return 0


def Bepipred3_analysis(results_bepipred_dir,name,thres_bepi,base_neg=1,amino_region=10,secondary_score=0,bepipred_opt=1):
    os.chdir(results_bepipred_dir)
    data=pd.read_csv('raw_output.csv')
    data=data.reset_index(drop=True)
    scores=data['BepiPred-3.0 score']
    mask=data['BepiPred-3.0 score']>thres_bepi
    mask_values=mask.values
    inds=[]
    groups=[]
    if bepipred_opt==1:
      for k in range(len(mask_values)):
          if mask_values[k] == True:
              inds.append(k)
      for i in range(len(inds)):
          if i==0:
              groups.append([inds[i]])
              counter=0
          else:
              if secondary_score == 0:
                  cond = inds[i] <= groups[counter][-1]+1+base_neg
              else:
                  cond = (inds[i] <= groups[counter][-1]+1+base_neg)& (scores[inds[i]] > secondary_score)
              
              if cond :
                  if len(groups[counter])==1:
                      groups[counter].append(inds[i])
                  else:
                      groups[counter][-1]=inds[i]
              else:
                  groups.append([inds[i]])
                  counter=counter+1
    else:
      pattern=[True]*amino_region
      for i in range(len(mask_values) - len(pattern) + 1):
        if all(mask_values[i:i+len(pattern)] == pattern):
            groups.append([i,i+len(pattern)-1])
              
    dd={'start':[],'end':[],'length':[]}
    for k in range(len(groups)):
        if len(groups[k])>1:
            start=groups[k][0]
            end=groups[k][1]
            length=abs(end-start)
        else:
            start=groups[k][0]
            end=groups[k][0]
            length=1
        dd['start'].append(start)
        dd['end'].append(end)
        dd['length'].append(length)
    df=pd.DataFrame.from_dict(dd)

    df_filtered=df[df['length']>=amino_region -2]
    start=df_filtered['start'].values
    end=df_filtered['end'].values
    final={'Sequences':[],'Bepipred Scores':[],'Length Region':[],'Start':[],'End':[]}
    for j in range(len(start)):
        st=start[j]
        en=end[j]
        dd=data.loc[st:en,:]
        chars=dd['Residue'].values.tolist()
        seq=''.join(chars)
        scores=dd['BepiPred-3.0 score'].values.tolist()
        final['Sequences'].append(seq)
        final['Bepipred Scores'].append(scores)
        final['Length Region'].append(len(seq))
        final['Start'].append(st)
        final['End'].append(en)
    final=pd.DataFrame.from_dict(final)
    #Mean scores
    sc=final['Bepipred Scores'].values
    mean=[np.mean(i) for i in sc]
    final['Bepipred mean score']=mean
    final=final.drop('Bepipred Scores',axis=1)
        
    
    return final
