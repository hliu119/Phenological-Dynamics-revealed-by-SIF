# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 04:59:19 2020

@author: Administrator
"""

import pandas as pd
import os
import numpy as np
from datetime import datetime


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

site_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\sites.csv")
inroot = r"C:\Users\Administrator\Desktop\EVI&NIRv"
outroot = r"C:\Users\Administrator\Desktop"
day_length = 8

for i in range(len(site_pd)):
    site_name = site_pd.iloc[i,0] #获取站点名称
    if site_name.split("-")[0] == "US":
        region = "NorthAmerica"
    else:
        region = "Europe"
  
    file = os.path.join(inroot, "{0}.csv".format(site_name))
    dir_pd = pd.read_csv(file)
    list18 = []
    list19 = []
    for j in range(len(dir_pd)):
        Year = dir_pd.iloc[j,0].split(" ")[2]
        Month = dir_pd.iloc[j,0].split(" ")[0]
        Day = (dir_pd.iloc[j,0].split(" ")[1]).split(",")[0]
        Date = "{0}{1}{2}".format(Year, Month, Day)
        dt = datetime.strptime(Date, '%Y%b%d').date() 
        tt = dt.timetuple()
        if tt.tm_year == 2018:
            list18.append([tt.tm_yday, dir_pd.iloc[j,1], dir_pd.iloc[j,2]])
        if tt.tm_year == 2019:
            list19.append([tt.tm_yday, dir_pd.iloc[j,1], dir_pd.iloc[j,2]])

    arr_EVI = np.full((len(list18), 2), -9999.0) 
    arr_NIRv = np.full((len(list18), 2), -9999.0)         
    for k in range(len(list18)):
        arr_EVI[k,0] = list18[k][0]
        arr_NIRv[k,0] = list18[k][0]   
        if np.isnan(list18[k][1]):
            if np.isnan(list19[k][1]):
                continue
            else:
                arr_EVI[k,1] = list19[k][1]
        else:
            if np.isnan(list19[k][1]):
                arr_EVI[k,1] = list18[k][1]
            else:
                arr_EVI[k,1] = (list18[k][1] + list19[k][1])/2
                
        if np.isnan(list18[k][2]):
            if np.isnan(list19[k][2]):
                continue
            else:
                arr_NIRv[k,1] = list19[k][2]
        else:
            if np.isnan(list19[k][2]):
                arr_NIRv[k,1] = list18[k][2]
            else:
                arr_NIRv[k,1] = (list18[k][2] + list19[k][2])/2 
                
    arr_EVI_Merge = Merge_Ori(arr_EVI, day_length)
    arr_NIRv_Merge = Merge_Ori(arr_NIRv, day_length) 
    
    outfile_EVI = os.path.join(outroot, "Output_EVI\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name))
    outfile_NIRv= os.path.join(outroot, "Output_NIRv\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name))
    
    np.savetxt(outfile_EVI, arr_EVI_Merge, delimiter = ',')
    np.savetxt(outfile_NIRv, arr_NIRv_Merge, delimiter = ',')    