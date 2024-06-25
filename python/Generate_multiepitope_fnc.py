
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

def combination_function(combination,epitopes):
    new_list_with_all_epitopes=[]
    
    for new_epitope in epitopes:
        temp_list=[]
        if type(combination) is tuple:
            temp_list.append(combination)
        elif type(combination) is list:
            temp_list.extend(combination)
        
        temp_list.append(new_epitope)
        
        new_list_with_all_epitopes.append(temp_list)
    return new_list_with_all_epitopes

def epitopes_on_vaccine(epitopes_all,order,nums):
    b_cell=[]
    ctl=[]
    htl=[]
    all_epitopes_combined=[]
    if len(epitopes_all[0])==0:
      nums[0]=0
    elif len(epitopes_all[1])==0:
      nums[1]=0
    elif len(epitopes_all[2])==0:
      nums[2]=0
      
    for j in range(len(order)):
        epitopes=epitopes_all[order[j]]
        if len(epitopes)==0:
            continue
        if len(all_epitopes_combined)==0:
            for epitope in epitopes:
                all_epitopes_combined.append(epitope)
        else:
            new_list=[]
            for i in all_epitopes_combined:
                temp_epitopes=combination_function(i,epitopes)
                
                new_list.extend(temp_epitopes)
            all_epitopes_combined=new_list
    
    #b1=np.array(all_epitopes_combined)
    #b_cells=np.array([[i] for i in b1[:,0]]).reshape([b1[:,0].shape[0],len(b1[0,0])])
    #ctls=np.array([[i] for i in b1[:,1]]).reshape([b1[:,0].shape[0],len(b1[0,1])])
    #htls=np.array([[i] for i in b1[:,2]]).reshape([b1[:,0].shape[0],len(b1[0,2])])
    names=['B-cell','CTL','HTL']
    if nums[order[0]]==0:
          df={names[order[0]]:['' for row in all_epitopes_combined],names[order[1]]:[list(row[0]) for row in all_epitopes_combined],names[order[2]]:[list(row[1]) for row in all_epitopes_combined]}
    
    elif nums[order[1]]==0:
          df={names[order[0]]:[list(row[0]) for row in all_epitopes_combined],names[order[1]]:['' for row in all_epitopes_combined],names[order[2]]:[list(row[1]) for row in all_epitopes_combined]}

    elif nums[order[2]]==0:
          df={names[order[0]]:[list(row[0]) for row in all_epitopes_combined],names[order[1]]:[list(row[1]) for row in all_epitopes_combined],names[order[2]]:['' for row in all_epitopes_combined]}

    else:
      df={names[order[0]]:[list(row[0]) for row in all_epitopes_combined],names[order[1]]:[list(row[1]) for row in all_epitopes_combined],names[order[2]]:[list(row[2]) for row in all_epitopes_combined]}
    
    return df


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
        if len(b_seqs)>6:
            b_seqs=b_seqs[:6]
        else:
            b_seqs=b_seqs
        b_num=num_epitopes[0]
        if b_num > len(b_seqs):
          b_num=len(b_seqs)
        comb_b_seqs=generate_combinations(b_seqs,b_num)
        ####comb_b_seqs --- b cell epitope combinations
        comp_b=final_components(comb_b_seqs,b_linker)
    else:
        comp_b=''
        comb_b_seqs=[]

    if not c_df.empty:
        c_linker=linkers[1]
        c_num=num_epitopes[1]
        c_seqs=c_df['Sequences'].values
        if len(c_seqs)>6:
            c_seqs=c_seqs[:6]
        else:
            c_seqs=c_seqs
        if c_num>len(c_seqs):
          c_num=len(c_seqs)
        comb_c_seqs=generate_combinations(c_seqs,c_num)
        ####comb_c_seqs --- ctl epitope combinations
        comp_c=final_components(comb_c_seqs,c_linker)
        
    else:
        comp_c=''
        comb_c_seqs=[]

    if not h_df.empty:
        h_linker=linkers[2]
        h_num=num_epitopes[2]
        h_seqs=h_df['Sequences'].values
        if len(h_seqs)>6:
            h_seqs=h_seqs[:6]
        else:
            h_seqs=h_seqs
        if h_num > len(h_seqs):
          h_num=len(h_seqs)
        comb_h_seqs=generate_combinations(h_seqs,h_num)
        comp_h=final_components(comb_h_seqs,h_linker)
    else:
        comp_h=''
        comb_h_seqs=[]
    #Order
    order=order_dict[int(order_opt)]
    #Multiepitope Sequences
    vacc=multiepitope_vaccine([comp_b,comp_c,comp_h],order,N_term)
    epit=epitopes_on_vaccine([comb_b_seqs,comb_c_seqs,comb_h_seqs],order,num_epitopes)
    
    #Final DataFrame
    df=pd.DataFrame({'Sequences':vacc})
    df['B_cell']=epit['B-cell']
    df['CTL']=epit['CTL']
    df['HTL']=epit['HTL']
    return df
