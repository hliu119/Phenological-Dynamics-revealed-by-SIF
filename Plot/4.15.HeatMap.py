# -*- coding: utf-8 -*-
"""
Created on Sat Aug 22 02:52:36 2020

@author: Administrator
"""

import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
ftsize = 5
fig,ax = plt.subplots(figsize=(3, 3)) 
#Generate Color Map
sns.set(rc={'font.family':"Times New Roman"})

colormap = sns.diverging_palette(220, 10, as_cmap=True)
uniform_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\compare.csv", index_col = 0)
uniform_pd.columns = ['\u0394SOS_SC', '\u0394EOS_SC','\u0394GSL_SC',\
                      '\u0394SOS_SR', '\u0394EOS_SR','\u0394GSL_SR']                           
                                          
ax = sns.heatmap(uniform_pd,cmap=colormap,annot = True, annot_kws = {"size": ftsize}, cbar=False)
"""
cbar = ax.collections[0].colorbar
# here set the labelsize by 20
cbar.ax.tick_params(length = 0, labelsize= ftsize)
"""
#add the column names as labels
ax.set_yticklabels(uniform_pd.index, rotation = 0, fontsize = ftsize/1.2, family="Times New Roman")
ax.set_xticklabels(uniform_pd.columns, rotation = 0, fontsize = ftsize/1.2, family="Times New Roman")
ax.set_xlabel('Phenological Metrics', fontsize = ftsize*1.5, family="Times New Roman")
ax.set_ylabel('SiteName', fontsize = ftsize*1.5, family="Times New Roman")
Plot_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step4.Table\Polygon\4.14.HeatMap\Appendix1.jpg"
fig.tight_layout() 
plt.show()
fig.savefig(Plot_path, dpi=600, quality=100)