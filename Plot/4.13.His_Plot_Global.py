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
import matplotlib.gridspec as gridspec
import warnings
warnings.filterwarnings('ignore')

in_file_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step14.Boreal_Map"
landcover_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step14.Boreal_Map\MOD12_LandUse_N.tif"
plot_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step4.Table\Polygon\4.11.Box_Plot_His"

Products = ["TROPOMI_SC","TROPOMI_SR","NIRv","EVI"]
LandCover = ["DBF", "NF", "SAV", "GRA", "SHR", "CRO"]
LC_CODE = [[4,5],[1,3],[8,9],[10],[6,7],[12,14]]
Metrics = [1,3,5] #1,3,5分别代表SOS,EOS,LOS
Name = ["DBF", "NF", "SAV", "GRA", "SHR", "CRO"]
nrows = 2
ncols = 3
figure_scale_row = nrows * 2.0      
figure_scale_col = ncols * 2.0 
fig = plt.figure(figsize=(figure_scale_col, figure_scale_row))#, sharey=True)
gs = gridspec.GridSpec(nrows, ncols, wspace=0.1, hspace=0.6)  

fs = 5 # fontsize
x_labels = ['SOS', 'EOS', 'GSL'] #'Tropomi_SR', 'MODIS', 'OCO-2']   
y_lables = [0, 100, 200, 300]
min_y = 0
max_y = 365
ny = 4 #y轴刻度个数
bar_width = 0.15 #柱状体宽度
capsize = 1.2 #柱状体标准差参数1
capthick = 0.8  #柱状体标准差参数2
elinewidth = 0.8 #柱状体标准差参数3
linewidth = 1.0 #边框线宽度
ftsize = 10 #字体大小
ftfamily = "Times New Roman"
axlength = 2.0 #轴刻度长度
axwidth = 1.2 #轴刻度宽度
legendcols = 5 #图例一行的个数

for i in range(len(LandCover)):
    mean = []
    err = [] 
    for j in range(len(Metrics)):  
        lc_array = gdal_array.LoadFile(landcover_file)
        SIF_SC_array = gdal_array.LoadFile(os.path.join(in_file_path, "14.1.SC_GSL_Mask/SC_{0}.tif".format(Metrics[j])))
        SIF_SR_array = gdal_array.LoadFile(os.path.join(in_file_path, "14.1.SR_GSL_Mask/SR_{0}.tif".format(Metrics[j]))) 
        dcSIF_array = gdal_array.LoadFile(os.path.join(in_file_path, "14.1.dcSIF_GSL_Mask/dcSIF_{0}.tif".format(Metrics[j])))         
        EVI_array = gdal_array.LoadFile(os.path.join(in_file_path, "14.1.EVI_GSL_Mask/EVI_{0}.tif".format(Metrics[j])))
        NIRv_array = gdal_array.LoadFile(os.path.join(in_file_path, "14.1.NIRv_GSL_Mask/NIRv_{0}.tif".format(Metrics[j])))

        SIF_SC = []
        SIF_SR = [] 
        dcSIF = []         
        NIRv = []        
        EVI = []        

        for m in range(len(lc_array)):
            for n in range(len(lc_array[0])):
               if lc_array[m,n] in LC_CODE[i]:
                  if abs(SIF_SC_array[m,n]) < 999: 
                      SIF_SC.append(SIF_SC_array[m,n])
                  if abs(SIF_SR_array[m,n]) < 999: 
                      SIF_SR.append(SIF_SR_array[m,n])
                  if abs(dcSIF_array[m,n]) < 999: 
                      dcSIF.append(dcSIF_array[m,n])                      
                  if abs(NIRv_array[m,n]) < 999: 
                      NIRv.append(NIRv_array[m,n])                                       
                  if abs(EVI_array[m,n]) < 999:                   
                      EVI.append(EVI_array[m,n])
                  
        mean.append([np.nanmean(SIF_SC), np.nanmean(SIF_SR), np.nanmean(dcSIF), np.nanmean(NIRv), np.nanmean(EVI)])
        err.append([np.nanstd(SIF_SC), np.nanstd(SIF_SR), np.nanstd(dcSIF), np.nanstd(NIRv), np.nanstd(EVI)]) 
    mean = np.array(mean)
    err = np.array(err)
    
    #画直方图        
    x = np.arange(int(len(x_labels)))
    col = int(i % ncols)
    row = int(i / ncols)
    axes = fig.add_subplot(gs[row, col])
        
    SIF_SC = axes.bar(x + 0 * bar_width, mean[:,0], bar_width, yerr = err[:,0],  error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="olivedrab", label = "SIF$_T$$_R$$_O$$_P$$_O$$_M$$_I$$_\_$$_t$$_o$$_t$$_a$$_l$$_-$$_S$$_C$", align="center", alpha=1)
    SIF_SR = axes.bar(x + 1 * bar_width, mean[:,1], bar_width, yerr = err[:,1],  error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="yellowgreen",  label = "SIF$_T$$_R$$_O$$_P$$_O$$_M$$_I$$_\_$$_t$$_o$$_t$$_a$$_l$$_-$$_S$$_R$", align="center", alpha=1)
    dcSIF = axes.bar(x + 2 * bar_width, mean[:,2], bar_width, yerr = err[:,2],  error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="forestgreen",   label = "SIF$_T$$_R$$_O$$_P$$_O$$_M$$_I$$_\_$$_O$$_b$$_s$", align="center", alpha=1)
    NIRv = axes.bar(x + 3 * bar_width, mean[:,3], bar_width, yerr = err[:,3],  error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="darkgoldenrod",   label = "NIRv", align="center", alpha=1)
    EVI = axes.bar(x + 4 * bar_width, mean[:,4], bar_width, yerr = err[:,4],  error_kw = {'ecolor' : '0.2', 'elinewidth':elinewidth, 'capsize' :capsize, 'capthick' :capthick}, color="gold",   label = "EVI", align="center", alpha=1)

    axes.set_title("({0}) {1}".format(chr(97+i),Name[i]), fontsize = ftsize/1.2, family = ftfamily)  
    axes.set_xticks(x + 2 * bar_width)
    axes.set_xticklabels(x_labels, fontsize = fs, family = ftfamily)
    axes.set_ylim(min_y, max_y)
    axes.spines['left'].set_linewidth(linewidth)
    axes.spines['right'].set_linewidth(linewidth)
    axes.spines['top'].set_linewidth(linewidth)
    axes.spines['bottom'].set_linewidth(linewidth)
    axes.tick_params(axis='both', length = axlength, width = axwidth, labelsize = ftsize/1.5)
    if col == 0:
        axes.set_ylabel('Day of Year (days)', fontsize = ftsize/1.6, family=ftfamily)
        axes.set_yticks(np.linspace(min_y, max_y - 65, ny))
        axes.set_yticklabels(y_lables, fontsize = fs + 2, family=ftfamily)
    else:
        axes.yaxis.set_visible(False)
    axes.set_xlabel('Phenological Metrics', fontsize = ftsize/1.5, family=ftfamily)

    handles = [SIF_SC, SIF_SR, dcSIF, NIRv, EVI]
   
    labels = ['SIF$_\mathdefault{total-SC}$',\
              'SIF$_\mathdefault{total-SR}$',\
              'SIF$_\mathdefault{Obs}$',\
              'NIR$_\mathdefault{V}$','EVI']        
    """    
    if i == 0:
        axes.legend(handles, labels, loc ='upper left', fancybox = False, shadow = False,frameon = False, 
                   ncol = legendcols, prop={'family':ftfamily, 'size':ftsize/3})
    """


fig.legend(handles, labels, loc ='lower center', fancybox = False, shadow = False,frameon = False, 
          ncol = legendcols, handletextpad = 0.2, columnspacing = 1.5, prop={'family':"Times New Roman", 'size':ftsize/1.3})  
fig.tight_layout() 
fig.subplots_adjust(left = None, right = None, bottom = 0.15)
Plot_path = os.path.join(plot_path, "Rs2-Global.jpg")
plt.show()
fig.savefig(Plot_path, dpi=600, quality=100,bbox_inches='tight')
        