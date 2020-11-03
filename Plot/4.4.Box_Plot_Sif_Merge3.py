# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 14:58:57 2019

@author: Administrator
"""
import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font',family='Times New Roman')
import os 
import warnings
warnings.filterwarnings('ignore')

in_file = r"C:\Users\Administrator\Desktop\Paper\statistic.csv"
plot_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step4.Table\Polygon\4.4.Box_Plot_Merge"

Table = pd.read_csv(in_file, index_col = 0 )    
Errors = []
"""
计算MAE，Mean，Median值
"""
def Statistical_Index(arr):
    arr = arr[~np.isnan(arr)]
    MAE = np.nanmean(abs(arr))
    MEAN = np.nanmean(arr)
    MEDIAN = np.nanmedian(arr)
    return (round(MAE, 2),round(MEAN, 2),round(MEDIAN, 2))

#sos       
Tro_SC_Sta_sos = np.array(Table["SC-GPP_SOS"]).reshape(-1,1)
Tro_SR_Sta_sos = np.array(Table["SR-GPP_SOS"]).reshape(-1,1)                
Tro_Sta_sos = np.array(Table["Obs-GPP_SOS"]).reshape(-1,1)
NIR_Sta_sos = np.array(Table["NIRv-GPP_SOS"]).reshape(-1,1)
Mod_Sta_sos = np.array(Table["EVI-GPP_SOS"]).reshape(-1,1)

Tro_SC_Sta_sos = np.delete(Tro_SC_Sta_sos, Errors, 0) 
Tro_SR_Sta_sos = np.delete(Tro_SR_Sta_sos, Errors, 0)
Tro_Sta_sos = np.delete(Tro_Sta_sos, Errors, 0)     
NIR_Sta_sos = np.delete(NIR_Sta_sos, Errors, 0)
Mod_Sta_sos = np.delete(Mod_Sta_sos, Errors, 0)                     
print("Tro_SC_sos: {0}".format(Statistical_Index(Tro_SC_Sta_sos)))
print("Tro_SR_sos: {0}".format(Statistical_Index(Tro_SR_Sta_sos)))          
print("Tro_sos:    {0}".format(Statistical_Index(Tro_Sta_sos))) 
print("NIR_sos:    {0}".format(Statistical_Index(NIR_Sta_sos))) 
print("EVI_sos:    {0}".format(Statistical_Index(Mod_Sta_sos))) 
print("\n") 

data_sos = [Mod_Sta_sos[~np.isnan(Mod_Sta_sos)],\
            NIR_Sta_sos[~np.isnan(NIR_Sta_sos)],\
            Tro_Sta_sos[~np.isnan(Tro_Sta_sos)],\
            Tro_SC_Sta_sos[~np.isnan(Tro_SC_Sta_sos)],\
            Tro_SR_Sta_sos[~np.isnan(Tro_SR_Sta_sos)]]

#eos
Tro_SC_Sta_eos = np.array(Table["SC-GPP_EOS"]).reshape(-1,1)
Tro_SR_Sta_eos = np.array(Table["SR-GPP_EOS"]).reshape(-1,1)                
Tro_Sta_eos = np.array(Table["Obs-GPP_EOS"]).reshape(-1,1)
NIR_Sta_eos = np.array(Table["NIRv-GPP_EOS"]).reshape(-1,1)
Mod_Sta_eos = np.array(Table["EVI-GPP_EOS"]).reshape(-1,1)

Tro_SC_Sta_eos = np.delete(Tro_SC_Sta_eos, Errors, 0)
Tro_SR_Sta_eos = np.delete(Tro_SR_Sta_eos, Errors, 0)
Tro_Sta_eos = np.delete(Tro_Sta_eos, Errors, 0)     
NIR_Sta_eos = np.delete(NIR_Sta_eos, Errors, 0)
Mod_Sta_eos = np.delete(Mod_Sta_eos, Errors, 0)      
         
print("Tro_SC_eos: {0}".format(Statistical_Index(Tro_SC_Sta_eos)))
print("Tro_SR_eos: {0}".format(Statistical_Index(Tro_SR_Sta_eos)))          
print("Tro_eos:    {0}".format(Statistical_Index(Tro_Sta_eos))) 
print("NIR_eos:    {0}".format(Statistical_Index(NIR_Sta_eos))) 
print("EVI_eos:    {0}".format(Statistical_Index(Mod_Sta_eos))) 
print("\n") 

data_eos = [Mod_Sta_eos[~np.isnan(Mod_Sta_eos)],\
            NIR_Sta_eos[~np.isnan(NIR_Sta_eos)],\
            Tro_Sta_eos[~np.isnan(Tro_Sta_eos)],\
            Tro_SC_Sta_eos[~np.isnan(Tro_SC_Sta_eos)],\
            Tro_SR_Sta_eos[~np.isnan(Tro_SR_Sta_eos)]]

#eos
Tro_SC_Sta_los = np.array(Table["SC-GPP_GSL"]).reshape(-1,1)
Tro_SR_Sta_los = np.array(Table["SR-GPP_GSL"]).reshape(-1,1)                
Tro_Sta_los = np.array(Table["Obs-GPP_GSL"]).reshape(-1,1)
NIR_Sta_los = np.array(Table["NIRv-GPP_GSL"]).reshape(-1,1)
Mod_Sta_los = np.array(Table["EVI-GPP_GSL"]).reshape(-1,1)

Tro_SC_Sta_los = np.delete(Tro_SC_Sta_los, Errors, 0)
Tro_SR_Sta_los = np.delete(Tro_SR_Sta_los, Errors, 0)
Tro_Sta_los = np.delete(Tro_Sta_los, Errors, 0)                      
NIR_Sta_los = np.delete(NIR_Sta_los, Errors, 0)
Mod_Sta_los = np.delete(Mod_Sta_los, Errors, 0)   

print("Tro_SC_los: {0}".format(Statistical_Index(Tro_SC_Sta_los)))
print("Tro_SR_los: {0}".format(Statistical_Index(Tro_SR_Sta_los)))          
print("Tro_los:    {0}".format(Statistical_Index(Tro_Sta_los))) 
print("NIR_los:    {0}".format(Statistical_Index(NIR_Sta_los))) 
print("EVI_los:    {0}".format(Statistical_Index(Mod_Sta_los))) 
print("\n") 


data_los = [Mod_Sta_los[~np.isnan(Mod_Sta_los)],\
            NIR_Sta_los[~np.isnan(NIR_Sta_los)],\
            Tro_Sta_los[~np.isnan(Tro_Sta_los)],\
            Tro_SC_Sta_los[~np.isnan(Tro_SC_Sta_los)],\
            Tro_SR_Sta_los[~np.isnan(Tro_SR_Sta_los)]]
matplotlib.rc('font',family='Times New Roman')
fig, axes = plt.subplots(nrows=1, ncols=3, figsize=(12, 4))#, sharey=True)
fs = 16 # fontsize
labels = ['EVI','NIR$_\mathdefault{V}$',\
          'SIF$_\mathdefault{Obs}$',\
          'SIF$_\mathdefault{total-SC}$',\
          'SIF$_\mathdefault{total-SR}$']  
lim_y = 55
inter_y = 25
gsl_offset = 50
line_width = 1  
rotation_angle = 45

axes[0].boxplot(data_sos, patch_artist=False, boxprops = {'color':'blue'}, flierprops = {'markeredgecolor':'b'}, \
                meanline =True, showmeans=True, capprops = {'color':'blue'}, whiskerprops={'color':'blue'}, showfliers=False)
axes[0].set_title('(a) Start of season', fontsize=fs, family="Times New Roman")
axes[0].set_xticklabels(labels, fontsize = fs/1.2, family="Times New Roman", rotation = rotation_angle, ha= "center")
axes[0].set_ylabel('Difference (Days)',fontsize = fs, family="Times New Roman")
axes[0].set_ylim(-lim_y, lim_y)
axes[0].set_yticks(np.arange(-lim_y + 5, lim_y , inter_y))   
axes[0].set_yticklabels(np.arange(-lim_y+ 5, lim_y, inter_y), fontsize = fs, family="Times New Roman")                 
axes[0].axhline(linewidth=line_width, color='black', linestyle = "dashed", fillstyle = "bottom")        
#axes[0].tick_params(axis='both', labelsize = 12) 

axes[1].boxplot(data_eos,patch_artist=False, boxprops = {'color':'blue'}, flierprops = {'markeredgecolor':'b'}, \
                meanline =True, showmeans=True, capprops = {'color':'blue'}, whiskerprops={'color':'blue'}, showfliers=False)
axes[1].set_title('(b) End of season', fontsize=fs, family="Times New Roman")
axes[1].set_xticklabels(labels, fontsize = fs/1.2, family="Times New Roman", rotation = rotation_angle, ha= "center") 
axes[1].set_ylim(-lim_y, lim_y)
axes[1].set_yticks(np.arange(-lim_y + 5, lim_y , inter_y))   
axes[1].set_yticklabels(np.arange(-lim_y+ 5, lim_y, inter_y), fontsize = fs, family="Times New Roman")    
"""
axes[1].set_ylim(-lim_y - 10, lim_y + 10)
axes[1].set_yticks(np.arange(-lim_y - 5, lim_y + 10, inter_y + 5))   
axes[1].set_yticklabels(np.arange(-lim_y - 5, lim_y + 10, inter_y + 5), fontsize = fs, family="Times New Roman") 
"""                
axes[1].axhline(linewidth=line_width, color='black', linestyle = "dashed", fillstyle = "bottom")        
#axes[1].tick_params(axis='both', labelsize = 12) 
            
axes[2].boxplot(data_los, patch_artist=False, boxprops = {'color':'blue'},flierprops = {'markeredgecolor':'b'}, \
                meanline =True, showmeans=True, capprops = {'color':'blue'}, whiskerprops={'color':'blue'}, showfliers=False)
axes[2].set_title('(c) Growing Season Length', fontsize=fs, family="Times New Roman")
axes[2].set_xticklabels(labels, fontsize = fs/1.2, family="Times New Roman", rotation = rotation_angle, ha= "center")
axes[2].set_ylim(-lim_y - gsl_offset, lim_y + gsl_offset)
axes[2].set_yticks(np.arange(-lim_y - gsl_offset + 5, lim_y + gsl_offset + 5, inter_y*2))   
axes[2].set_yticklabels(np.arange(-lim_y - gsl_offset + 5, lim_y + gsl_offset + 5, inter_y*2), fontsize = fs, family="Times New Roman")                 
axes[2].axhline(linewidth=line_width, color='black', linestyle = "dashed", fillstyle = "bottom")        
#axes[2].tick_params(axis='both', labelsize = 12)     

fig.subplots_adjust()
fig.tight_layout() 
Plot_path = os.path.join(plot_path, "DoubleLogBeck_Trs_Boreal.jpg")
plt.show()  
fig.savefig(Plot_path,dpi=600,quality=100)        
