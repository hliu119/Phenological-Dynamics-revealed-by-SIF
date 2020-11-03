# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 04:59:19 2020

@author: Administrator
"""

import pandas as pd
import os
import numpy as np

def Merge_Ori(array, day_length):
    Year_Period = 365
    arr_merge = np.full((int(Year_Period/day_length), 2), -9999.0)    
    for m in range(0, int(Year_Period/day_length)):
        count = 0
        value = 0
        for n in range(0, len(array[:,0])):       
            if array[n,0] >= m * day_length and array[n,0] < (m + 1) * day_length:  
                if not np.isnan(array[n,1]):
                    if int(array[n,1]) != -9999:
                       count += 1
                       value += array[n,1]
        arr_merge[m,0] = (m + 1) * day_length                       
        if count != 0:
            arr_merge[m,1] = value/count 
    return  arr_merge 

Datasets = ["GPP"]#["SC","SR","Obs"]
site_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\Merge_Info.csv")
inroot = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version0\1Day"
outroot = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version0\8Day"
day_length = 8

for i in range(len(site_pd)):
    site_name = site_pd.iloc[i,0] #获取站点名称
    if site_name.split("-")[0] == "US":
        region = "NorthAmerica"
    else:
        region = "Europe"
    for dataset in Datasets:    
        file = os.path.join(inroot, "Output_{0}\{1}\SitePixel\Polygon\{2}.csv".format(dataset,region,site_name))
        arr = np.loadtxt(file, delimiter = ",")    
        dist_arr = Merge_Ori(arr, day_length)
        outfile = os.path.join(outroot, "Output_{0}\{1}\SitePixel\Polygon\{2}.csv".format(dataset,region,site_name))
        np.savetxt(outfile, dist_arr, delimiter = ',')
    
