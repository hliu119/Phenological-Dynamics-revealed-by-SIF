# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 18:38:32 2019

@author: HaoranLiu
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')
    
def jd_to_time(fesc_18, fesc_19, fesc_month):
    fesc_18[:,0] = fesc_18[:,0] + 2018000
    fesc_19[:,0] = fesc_19[:,0] + 2019000   
    fesc = np.vstack((fesc_18, fesc_19))
    for m in range(1, 13):
        Count = 0 
        Sum = 0
        for i in range(len(fesc)):            
            dt = datetime.strptime(str(int(fesc[i,0])), '%Y%j').date()
            month = int(dt.strftime('%m'))           
            if month == m and not np.isnan(fesc[i,3]):
                Count = Count + 1
                Sum = Sum + fesc[i,3]
        dt = datetime.strptime("2018.{0}.15".format(m), '%Y.%m.%d')
        tt = dt.timetuple()
        fesc_month[m - 1, 0] = tt.tm_yday        
        if Count == 0:
            fesc_month[m - 1, 1] = np.nan             
        else:
            fesc_month[m - 1, 1] = Sum/Count 
    return fesc_month

def Merge_Ori(array, day_length):
    Year_Period = 365
    arr_merge = np.full((int(Year_Period/day_length), 2), np.nan)    
    for m in range(0, int(Year_Period/day_length)):
        count = 0
        value = 0
        for n in range(0, len(array[:,0])):       
            if array[n,0] >= m * day_length and array[n,0] < (m + 1) * day_length:  
                if not np.isnan(array[n,1]):
                    if array[n,1] != -9999:
                       count += 1
                       value += array[n,1]
        if count != 0:
            arr_merge[m,0] = (m + 1) * day_length
            arr_merge[m,1] = value/count 
    return  arr_merge 

plot_path = r"C:\Users\Administrator\Desktop\Paper\Plot_Fesc"
day_length = 1
day_length_points = 1
nrows = 1
ncols = 1   
figure_scale_row = nrows * 5      
figure_scale_col = ncols * 7     

def smooth(y, box_pts):
    box = np.ones(box_pts)/box_pts
    y_smooth = np.convolve(y, box, mode='same')
    return y_smooth

site_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\Paper\statistic_all.csv")
for i in range(len(site_pd)):
    site_name = site_pd.iloc[i,0] #获取站点名称
    print(site_name)
    #landtype = site_pd.iloc[i,2] #获取站点名称

    if site_name.split("-")[0] == "US":
        region = "NorthAmerica"
    else:
        region = "Europe"        

    """
    数据的原始数据
    """
    
    tropomi_ori_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Obs_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)
    tropomi_sc_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SC_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)
    tropomi_sr_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SR_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)   

    tropomi_ori_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Obs_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)
    tropomi_sc_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SC_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)
    tropomi_sr_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SR_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,site_name)   

    """
    将原始数据读入数组
    """  
             
    tropomi_ori_18_arr = np.loadtxt(tropomi_ori_18_file, delimiter = ",")
    tropomi_sc_18_arr = np.loadtxt(tropomi_sc_18_file, delimiter = ",")       
    #tropomi_sr_18_arr = np.loadtxt(tropomi_sr_18_file, delimiter = ",")       

    tropomi_ori_19_arr = np.loadtxt(tropomi_ori_19_file, delimiter = ",")
    tropomi_sc_19_arr = np.loadtxt(tropomi_sc_19_file, delimiter = ",")       
    #tropomi_sr_19_arr = np.loadtxt(tropomi_sr_19_file, delimiter = ",")

    merge_18 = np.hstack((tropomi_ori_18_arr, tropomi_sc_18_arr[:,1].reshape(-1,1)))                                                                                              
    merge_19 = np.hstack((tropomi_ori_19_arr, tropomi_sc_19_arr[:,1].reshape(-1,1)))    
    
    merge_18 = np.delete(merge_18, np.where(merge_18[:,1] == -9999), axis = 0) 
    merge_18 = np.delete(merge_18, np.where(merge_18[:,2] == -9999), axis = 0) 
    
    merge_19 = np.delete(merge_19, np.where(merge_19[:,1] == -9999), axis = 0)
    merge_19 = np.delete(merge_19, np.where(merge_19[:,2] == -9999), axis = 0)  
    
    fesc_18 = np.hstack((merge_18, (merge_18[:,1]/merge_18[:,2]).reshape(-1,1)))
    fesc_19 = np.hstack((merge_19, (merge_19[:,1]/merge_19[:,2]).reshape(-1,1))) 

    fesc_18 = np.delete(fesc_18, np.where(abs(fesc_18[:,3]) > 1), axis = 0)    
    fesc_19 = np.delete(fesc_19, np.where(abs(fesc_19[:,3]) > 1), axis = 0)  

    fesc = np.vstack((fesc_18, fesc_19)) 
    fesc_merge = Merge_Ori(fesc[:,[0,3]], 8)
    fesc_month = np.full((12,2), -9999.0)
    fesc_month = jd_to_time(fesc_18, fesc_19, fesc_month)
    #fesc排列顺序是 date, Obs, SC, GPP, fesc
    
    fig = plt.figure()
    fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize=(figure_scale_col, figure_scale_row))#, sharey=True)
    
    axes.scatter(fesc[:,0], fesc[:,3],  color='', marker='o', edgecolors='black', s = 100) # 把 corlor 设置为空，通过edgecolors来控制颜色 
    axes.scatter(fesc_merge[:,0], fesc_merge[:,1],  color='red', marker='x', s = 15) # 把 corlor 设置为空，通过edgecolors来控制颜色 
    axes.plot(fesc_month[:,0], fesc_month[:,1],  color='red') # 把 corlor 设置为空，通过edgecolors来控制颜色 
    # rft = np.fft.rfft(y)  
    axes.spines['left'].set_linewidth(2)
    axes.spines['right'].set_linewidth(2)
    axes.spines['top'].set_linewidth(2)
    axes.spines['bottom'].set_linewidth(2)
    axes.tick_params(axis='both', length = 3, width = 2, labelsize = 12)     
    #axes[0].legend(framealpha  = 0)    
    axes.set_xticks([15, 46, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349])
    axes.set_xticklabels(('Jan', 'Feb', 'Mar','Apr','May','Jun','Jul','Aug','Sep', 'Oct','Nov', 'Dec'))    
    axes.set_ylabel('Normalized Values', fontsize = 16)

    axes.set_xlabel("{0}".format(site_name), fontsize = 20)
                                  
    Plot_path = os.path.join(plot_path, "{0}_fesc.jpg".format(site_name))
    plt.show()
    fig.savefig(Plot_path, dpi=600, quality=100)
    fesc_merge[:,1] = 1/fesc_merge[:,1]
    #np.savetxt(r"C:\Users\Administrator\Desktop\Test\{0}.csv".format(site_name), fesc_merge, delimiter=',')
