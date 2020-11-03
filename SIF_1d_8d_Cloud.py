f# -*- coding: utf-8 -*-
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

Datasets = ["SC","SR","Obs"]
Years = [2018,2019]
site_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\Merge_Info.csv")
inroot = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version4\1Day"
outroot = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version4\8Day"
day_length = 8
cloud_threshold = 0.1

for i in range(len(site_pd)):
    site_name = site_pd.iloc[i,0] #获取站点名称
    if site_name.split("-")[0] == "US":
        region = "NorthAmerica"
    else:
        region = "Europe"
    for dataset in Datasets: 
        _arr = np.full((365,2),-9999.0)
        for year in Years:
            file = os.path.join(inroot, "Output_{0}_{1}\{2}\SitePixel\Polygon\{3}.csv".format(dataset, year, region, site_name))
            arr = np.loadtxt(file, delimiter = ",") 
            cloud_file = os.path.join(inroot, "Output_Cloud_{0}\{1}\SitePixel\Polygon\{2}.csv".format(year, region, site_name))
            cloud_arr = np.loadtxt(cloud_file, delimiter = ",")             
            for j in range(len(_arr)):
                _arr[j,0] = arr[j,0]            
                if (int(arr[j,1]) != -9999) and (cloud_arr[j,1] < cloud_threshold):
                    print(arr[j,1], cloud_arr[j,1])
                    if _arr[j,1] == -9999:
                        _arr[j,1] = arr[j,1]
                    else:
                        _arr[j,1] = (_arr[j,1] + arr[j,1])/2                
        dist_arr = Merge_Ori(_arr, day_length)
        outPath = os.path.join(outroot, "Output_{0}\{1}\SitePixel\Polygon".format(dataset,region)) 
        if not os.path.exists(outPath):
            os.makedirs(outPath)
        outfile = os.path.join(outPath, "{0}.csv".format(site_name))
        np.savetxt(outfile, dist_arr, delimiter = ',')
        
