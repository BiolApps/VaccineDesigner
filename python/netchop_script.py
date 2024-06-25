###Packages
import os
import pandas as pd
import numpy as np
def netchop_fnc(path,vaccine_sequence):
    data=netchop_analysis(path,vaccine_sequence)
    return data
def filtering_fnc(string_interest):
    string_as=string_interest.split(' ')
    string_f=[i for i in string_as if i!='']
    return string_f
def filter_lines(lines):
    final_lines=[i for i in lines if len(i)==5]
    return final_lines
def transform_line(row):
    row_0=int(row[0])
    row_3=float(row[3])
    row[0]=row_0
    row[3]=row_3
    return row

def netchop_analysis(results_path,vaccine_sequence):
    with open(results_path,'r') as f:
        lin=f.readlines()
    dict_new={}

    lines_filtered=[filtering_fnc(i) for i in lin]
    final_lines=filter_lines(lines_filtered)
    tabs=final_lines[0]
    tabs[-1]=tabs[-1][:-1]
    data=final_lines[1:]
    data_numpy =np.array(data)
    rows=[]
    for row in data_numpy:
        try:
            new_row=transform_line(row)
            rows.append(new_row)
        except:
            continue
    data_numpy=np.array(rows)
    data_numpy[:,-1]=[i[:-1] for i in data_numpy[:,-1]]
    dict_new={}
    for index,tab in enumerate(tabs):
        dict_new[tab]=data_numpy[:,index]

    df=pd.DataFrame.from_dict(dict_new)
    final=df[['pos','AA','C','score']]
    return final
