import os
import textwrap
from bs4 import BeautifulSoup
from pprint import pprint
from urllib.parse import urljoin
import webbrowser
import sys
from requests_html import HTMLSession
import re
import pandas as pd
import time
from requests.exceptions import ReadTimeout


def Prot_param_fnc(protein_seq):
    data={'prot_id': '', 'mandatory': '', None: '', 'sequence': protein_seq}
    success=False
    #Protparam URL
    url='https://web.expasy.org/cgi-bin/protparam/protparam'
    while not success:
        try:
        # pprint(data)
            session=HTMLSession()
            res = session.post(url, data=data,timeout=20)
            ##Protparam text
            text=res.html.full_text
            res.close()
            success=True
        except TimeoutError:
            time.sleep(0.1)
            session.close()
        except ReadTimeout:
            time.sleep(0.1)
            session.close()
        #Parsing to text of html page
    #Molecular weight
    molecular_weight_string=re.findall('Molecular weight: -?[0-9]+\.[0-9]+',text)
    mol=float(re.findall('-?[0-9]+\.[0-9]+',molecular_weight_string[0])[0])
    #Instability index
    ind_inst=text.index('The instability index (II) is computed to be')
    tt=text[ind_inst:ind_inst+100]
    inst_str=re.findall('computed to be -?[0-9]+\.[0-9]+',tt)
    inst_score=float(re.findall('-?[0-9]+\.[0-9]+',inst_str[0])[0])
    #Instability label
    ind_lab=text.index('This classifies the protein as')
    tt=text[ind_lab:ind_lab+50]
    lab_str=re.findall('unstable',tt)
    
    if len(lab_str)==0:
        inst='Stable'
    else:
        inst='unstable'
    #Aliphatic index
    aliph_str=re.findall('Aliphatic index: -?[0-9]+\.[0-9]+',text)
    alip=float(re.findall('-?[0-9]+\.[0-9]+',aliph_str[0])[0])
    #GRAVY score
    gravy_str=re.findall('\(GRAVY\): -?[0-9]\.[0-9]+',text)
    gravy_score=float(re.findall('-?[0-9]+\.[0-9]+',gravy_str[0])[0])

    a={'Molecular weight':mol,'Instability index':inst_score,'Instability':inst,'Aliphatic index':alip,'GRAVY score':gravy_score}
       
    return a

def Protparam(seqs):
    mol=[]
    inst_ind=[]
    inst_lab=[]
    alip_ind=[]
    gravy=[]
    for k in seqs:
        res=Prot_param_fnc(k)
        mol.append(res['Molecular weight'])
        inst_ind.append(res['Instability index'])
        inst_lab.append(res['Instability'])
        alip_ind.append(res['Aliphatic index'])
        gravy.append(res['GRAVY score'])
    dataframe=pd.DataFrame({'Sequences':seqs,'Molecular weight':mol,'Instability index':inst_ind,'Instability':inst_lab,'Aliphatic index':alip_ind,'GRAVY score':gravy})
    return dataframe

# def Predictions_per_30(df,thres):
#     seqs=df.Sequences.values
#     final_preds=[]
#     final_pred=[]
#     for k in seqs:
#         seqs_parts=textwrap.wrap(k,30)
#         preds=predictions(seqs_parts)
#         pred=predictions([k])[0]
#         final_preds.append(preds)
#         final_pred.append(pred)
#     result=[True if i>thres else False for i in final_pred]
#     a=pd.DataFrame({'Sequences':seqs,'Prediction_per_30':final_preds,'Prediction_Score_DeepVac':final_pred,'Prediction_deepvacpred':result})
#     return a
