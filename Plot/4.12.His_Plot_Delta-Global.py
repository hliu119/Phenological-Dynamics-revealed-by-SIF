# -*- coding: utf-8 -*-
"""
Created on Thu Dec 26 14:58:57 2019

@author: Administrator

MCD12 LandCover Types：
DBF == 4：DBF,5：MF
EBF == 2：EBF
NF == 1：ENF,3：DNF
CRO == 12: CRO, 14: CRO&NV
GRA == 10: GRA
SHR == 6：CSH, 7：OSH
SAV == 8：WSA, 9：SAV
"""

import numpy as np
import matplotlib.pyplot as plt
import os 
from osgeo import gdal_array
import math
import warnings
warnings.filterwarnings('ignore')

in_file_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step14.Boreal_Map\14.2.Delta"
landcover_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step14.Boreal_Map\MOD12_LandUse_N.tif"
plot_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step4.Table\Polygon\4.10.Box_Plot_Delta"
SIF_Products = ["SC"]
VI_Products = ["EVI","NIRv"]
SIF_Name = ['SIF$_\mathdefault{total-SC}$']
VI_Name = ["EVI",'NIR$_\mathdefault{V}$']    
Metrics = [1,3,5] #1,3,5分别代表SOS,EOS,LOS

nrows = 2
ncols = 1
figure_scale_row = nrows * 1.5      
figure_scale_col = ncols * 4  
fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(figure_scale_col, figure_scale_row))#, sharey=True)
fs = 5 # fontsize
x_labels =  ["DBF", "NF", "SAV",  "GRA", "SHR", "CRO"]   
y_lables = [-100, -50, 0, 50, 100]
min_y = -120
max_y = 120
bar_width = 0.2 
"""
bar_width,capsize,capthick 误差线格式
"""
capsize = 1
capthick = 0.75
elinewidth = 0.75
linewidth = 0.8 #边框宽度
ftsize = 8 #字体大小
for i in range(len(SIF_Products)):
    for j in range(len(VI_Products)):
        mean = []
        err = []
        for k in range(len(Metrics)):
            lc_array = gdal_array.LoadFile(landcover_file)
            in_file = os.path.join(in_file_path, "{0}-{1}_{2}.tif".format(SIF_Products[i], VI_Products[j], Metrics[k]))    
            in_array = gdal_array.LoadFile(in_file)                 
            DBF = [] 
            NF = []
            CRO = [] 
            GRA = []
            SHR = []            
            SAV = []

            for m in range(len(lc_array)):
                for n in range(len(lc_array[0])):
                   if abs(in_array[m,n]) < 999: 
                      if lc_array[m,n] == 4 or lc_array[m,n] == 5:
                          DBF.append(in_array[m,n])
                      elif  lc_array[m,n] == 1 or lc_array[m,n] == 3: 
                          NF.append(in_array[m,n])                        
                      elif lc_array[m,n] == 12 or lc_array[m,n] == 14:  
                          CRO.append(in_array[m,n])                           
                      elif lc_array[m,n] == 10: 
                          GRA.append(in_array[m,n])                           
                      elif lc_array[m,n] == 6 or lc_array[m,n] == 7: 
                          SHR.append(in_array[m,n])
                      elif lc_array[m,n] == 8 or lc_array[m,n] == 9: 
                          SAV.append(in_array[m,n])

            mean.append([np.nanmean(DBF), np.nanmean(NF), np.nanmean(SAV),\
                         np.nanmean(GRA), np.nanmean(SHR), np.nanmean(CRO)])
               
            err.append([np.nanstd(DBF), np.nanstd(NF), np.nanstd(SAV),\
                        np.nanstd(GRA), np.nanstd(SHR), np.nanstd(CRO)])
            """
            #计算标准误
            err.append([np.nanstd(DBF)/math.sqrt(len(DBF)), np.nanstd(NF)/math.sqrt(len(NF)), np.nanstd(SAV)/math.sqrt(len(SAV)),\
                        np.nanstd(GRA)/math.sqrt(len(GRA)), np.nanstd(SHR)/math.sqrt(len(SHR)), np.nanstd(CRO)/math.sqrt(len(CRO))])                  
            """ 
        #画直方图        
        x = np.arange(int(len(x_labels)))
            
        SOS = axes[j].bar(x + 0 * bar_width, mean[0], bar_width, yerr = err[0], error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="green", label = "SOS", align="center", alpha=0.75)
        EOS = axes[j].bar(x + 1 * bar_width, mean[1], bar_width, yerr = err[1], error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="blue",  label = "EOS", align="center", alpha=0.75)
        LOS = axes[j].bar(x + 2 * bar_width, mean[2], bar_width, yerr = err[2], error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="red",   label = "LOS", align="center", alpha=0.75)
        
        handles = [SOS, EOS, LOS]
        labels = ['SOS','EOS', 'GSL'] 
        
        axes[j].set_title("({0}) {1} - {2}".format(chr(97+i*2+j),SIF_Name[i],VI_Name[j]), fontsize = ftsize, family="Times New Roman")  
        axes[j].set_xticks(x + 1 * bar_width)
        axes[j].set_yticks(y_lables)        
        axes[j].set_ylim(min_y, max_y)
        axes[j].set_xticklabels(x_labels, fontsize = fs, family="Times New Roman")
        axes[j].set_yticklabels(y_lables, fontsize = fs, family="Times New Roman")
        axes[j].spines['left'].set_linewidth(linewidth)
        axes[j].spines['right'].set_linewidth(linewidth)
        axes[j].spines['top'].set_linewidth(linewidth)
        axes[j].spines['bottom'].set_linewidth(linewidth)
        axes[j].tick_params(axis='both', length = 1.5, width = 0.8, labelsize = ftsize/1.5)
    
        axes[j].set_ylabel('Difference Values (days)', fontsize = ftsize/1.5, family="Times New Roman")
        axes[1].set_xlabel('Land Cover', fontsize = ftsize/1.2, family="Times New Roman")
        if i == 0 and j == 0:
            axes[j].legend(handles, labels, loc ='upper left', fancybox = False, shadow = False,frameon = False, 
                       columnspacing = 0.5, ncol = 3, prop={'family':"Times New Roman", 'size':ftsize/2})  
    
fig.tight_layout() 
fig.subplots_adjust(left = None, right = None, bottom=0.15)   

Plot_path = os.path.join(plot_path, "VIs-TROPOMI_Global.jpg")
plt.show()
fig.savefig(Plot_path, dpi=600, quality=100,bbox_inches='tight')
        