
###Script for the generation of multiepitope sequences
#Packages
import itertools
import os
import pandas as pd
import numpy as np

order_dict={1:[0,1,2],2:[0,2,1],3:[2,0,1],4:[2,1,0],5:[1,2,0],6:[1,0,2]}
##Functions
def generate_combinations(strings,num):
    combinations = []
    combinations.extend(itertools.permutations(strings, num))
    return combinations

def final_components(seqs,linker):
    final_comp=[]
    for j in range(len(seqs)):
        temp=seqs[j]
        str_temp=''
        for k in temp:
            str_temp=str_temp+linker+k
        final_comp.append(str_temp)
    return final_comp

def epitope_combinations(sequence,list_epitopes):
    new_seqs=[]
    for i in list_epitopes:
        temp=sequence+i
        new_seqs.append(temp)
    return new_seqs

def multiepitope_vaccine(components,order,start=''):
    vaccines=[]
    subunit=[]
    for j in range(len(order)):
        comp_j=components[order[j]]
        if len(comp_j)==0:
            continue
        if len(subunit)==0:
            for i in comp_j:
                subunit.append(i)
        else:
            new_ls=[]
            for i in subunit:
                temp_seqs=epitope_combinations(i,comp_j)
                new_ls.extend(temp_seqs)
            subunit=new_ls
    subunit=[start+i for i in subunit]
    return subunit

def multiepitope_fnc(dfs,num_epitopes,linkers,order_opt,N_term):
    
    #Dataframes
    b_df=dfs[0]
    c_df=dfs[1]
    h_df=dfs[2]
    if not b_df.empty:
        b_linker=linkers[0]
        b_seqs=b_df['Sequences'].values
        if len(b_seqs)>5:
            b_seqs=b_seqs[:5]
        else:
            b_seqs=b_seqs
        b_num=num_epitopes[0]
        if b_num > len(b_seqs):
          b_num=len(b_seqs)
        comb_b_seqs=generate_combinations(b_seqs,b_num)
        comp_b=final_components(comb_b_seqs,b_linker)
    else:
        comp_b=''
    if not c_df.empty:
        c_linker=linkers[1]
        c_num=num_epitopes[1]
        c_seqs=c_df['Sequences'].values
        if len(c_seqs)>5:
            c_seqs=c_seqs[:5]
        else:
            c_seqs=c_seqs
        if c_num>len(c_seqs):
          c_num=len(c_seqs)
        comb_c_seqs=generate_combinations(c_seqs,c_num)
        comp_c=final_components(comb_c_seqs,c_linker)
    else:
        comp_c=''
    if not h_df.empty:
        h_linker=linkers[2]
        h_num=num_epitopes[2]
        h_seqs=h_df['Sequences'].values
        if len(h_seqs)>5:
            h_seqs=h_seqs[:5]
        else:
            h_seqs=h_seqs
        if h_num > len(h_seqs):
          h_num=len(h_seqs)
        comb_h_seqs=generate_combinations(h_seqs,h_num)
        comp_h=final_components(comb_h_seqs,h_linker)
    else:
        comp_h=''
    #Order
    order=order_dict[int(order_opt)]
    #Multiepitope Sequences
    vacc=multiepitope_vaccine([comp_b,comp_c,comp_h],order,N_term)
    #Final DataFrame
    df=pd.DataFrame({'Sequences':vacc})
    return df
