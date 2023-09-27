#Packages
import os
import subprocess
import numpy as np
import pandas as pd
import time
import re
from bs4 import BeautifulSoup
from pprint import pprint
from urllib.parse import urljoin
import webbrowser
import sys
from requests_html import HTMLSession
from requests.exceptions import ReadTimeout


def Vaxijen_fnc(protein_seq,target='bacteria',threshold=0.4):
    success=False
    #Vaxijen URLs
    url='http://www.ddg-pharmfac.net/vaxijen/VaxiJen/VaxiJen.html'
    url_script='../scripts/VaxiJen_scripts/VaxiJen3.pl'
    #Parameters
    seq=protein_seq
    target=target
    url=urljoin(url,url_script)
    our_data={'uploaded_file': '', 'Verbose': '', 'SequenceOnOff': '', 'SummaryMode': '', 'threshold': str(threshold), 'reset': '', 'Target': target, 'seq':seq}
    while not success:
        try:
        # pprint(data)
            session=HTMLSession()
            res = session.post(url, data=our_data,timeout=20)
            ##Vaxijen results text
            text=res.html.full_text
            res.close()
            success=True
        except TimeoutError:
            time.sleep(5)
            session.close()
        except ReadTimeout:
            time.sleep(5)
            session.close()
    #get the numbers with decimals from the text
    numbers=re.findall('-?[0-9]+\.[0-9]+',text)
    #prediction score
    pred_score=float(numbers[-1])
    #probable antigen or not
    indication=len(re.findall('NON-ANTIGEN',text))
    if indication==0:
        antigen=1
    else:
        antigen=0
    a={'prediction_score':pred_score,'ANTIGEN':antigen}
        
    return a


def Vaxijen(df,target,threshold):
    seqs=df['Sequences']
    pred_score=[]
    antigen=[]
    for k in seqs:
        temp=Vaxijen_fnc(k,target,threshold)
        pred_score.append(temp['prediction_score'])
        antigen.append(bool(temp['ANTIGEN']))
        time.sleep(1)

    data={'Sequences':seqs,'Prediction_Score':pred_score,'Antigen':antigen}
    df=pd.DataFrame.from_dict(data)
    
    return df
