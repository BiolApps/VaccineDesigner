#Packages
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import tkinter as tk
from tkinter import Tk,filedialog
import pandas as pd

a='C:/Users/30697/Desktop/Rotation/Automation_of_web_tools/Local_Tools/toxinpred2/toxinpred2'
dir='C:/Users/30697/Desktop/Rotation/Automation_of_web_tools/webserver/webserver/temp_data/Toxinpred'
df=pd.DataFrame({'Sequences':['NNQNQNQNQNQNQ','QNQNQNQNQNQNQNQQ']})
def Toxinpred_fnc(df,path_toxinpred,results_toxinpred_dir,model,thresh_toxin):
    path=os.getcwd()
    seq=df['Sequences']
    with open('temp.fasta','w+') as f:
        for k in range(len(seq)):
            f.write('>seq'+str(k)+'\n')
            f.write(seq[k]+'\n')
    file_path=path+'\\temp.fasta'
    #Toxinpred
    os.chdir(path_toxinpred)
    command = 'python toxinpred2.py -i '+file_path+' '+' -m '+str(model)+' -t '+str(thresh_toxin)+ ' -o '+ results_toxinpred_dir+'/results_toxinpred2.csv'
    os.system(command)
    os.chdir(results_toxinpred_dir)
    res=pd.read_csv('results_toxinpred2.csv')
    os.chdir(path)
    return res

