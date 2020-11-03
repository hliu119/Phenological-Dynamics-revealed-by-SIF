# -*- coding: utf-8 -*-
"""
Created on Mon Dec  2 18:38:32 2019

@author: HaoranLiu
"""
import scipy.stats as stats
import matplotlib
import matplotlib.pyplot as plt
matplotlib.rc('font',family='Times New Roman')
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
                Sum = Sum + fesc[i,4]
        dt = datetime.strptime("2018.{0}.15".format(m), '%Y.%m.%d')
        tt = dt.timetuple()
        fesc_month[m - 1, 0] = tt.tm_yday        
        if Count == 0:
            fesc_month[m - 1, 1] = np.nan             
        else:
            fesc_month[m - 1, 1] = Sum/Count 
    return fesc_month
    
def Normalization(array):
    x = array[:,0]
    y = (array[:,1] - np.nanmin(array[:,1]))/(np.nanmax(array[:,1]) - np.nanmin(array[:,1]))
    return x,y
    
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


def Merge_Fit(array, day_length):
    Year_Period = 365
    arr_merge = np.full((int(Year_Period/day_length), 2), np.nan)    
    for m in range(0, int(Year_Period/day_length)):
        count = 0
        value = 0
        for n in range(0, len(array)):       
            if n >= m * day_length and n < (m + 1) * day_length:  
                if not np.isnan(array[n]):
                    if array[n] != -9999:
                       count += 1
                       value += array[n]
        if count != 0:
            arr_merge[m,0] = (m + 1) * day_length
            arr_merge[m,1] = value/count 
    return  arr_merge 


plot_path = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step3.plot\Polygon\plot_paper"
day_length = 1
day_length_points = 1
nrows = 4
ncols = 4   
figure_scale_row = nrows * 2     
figure_scale_col = ncols * 2
vline_width = 0.8
ftsize = 10
labelsize = 6
annotate_size = 8
font = {'family' : 'Times New Roman'
}
linewidth = 1.5
pointsize = 2.5

GPP_Color = 'black'
SC_Color = 'red'
SR_Color = 'darkorange'
Obs_Color = 'blue'

fig = plt.figure()
fig, axes = plt.subplots(nrows = nrows, ncols = ncols, figsize=(figure_scale_col, figure_scale_row))#, sharey=True)
site_pd = pd.read_csv(r"C:\Users\Administrator\Desktop\Paper\statistic.csv")

for i in range(len(site_pd)):
    siteName = site_pd.iloc[i,0] #获取站点名称
    siteType = site_pd.iloc[i,1] #获取站点地物类型
    #landtype = site_pd.iloc[i,2] #获取站点名称
    """
    用遥感值减去GPP的值
    """
    TRO_SC_SOS = site_pd.at[i, 'SC-GPP_SOS']
    TRO_SC_EOS = site_pd.at[i, 'SC-GPP_EOS']
    TRO_SC_GSL = site_pd.at[i, 'SC-GPP_GSL']
    TRO_SR_SOS = site_pd.at[i, 'SR-GPP_SOS']
    TRO_SR_EOS = site_pd.at[i, 'SR-GPP_EOS']
    TRO_SR_GSL = site_pd.at[i, 'SR-GPP_GSL']
    TRO_Obs_SOS = site_pd.at[i, 'Obs-GPP_SOS']	
    TRO_Obs_EOS = site_pd.at[i, 'Obs-GPP_EOS']	
    TRO_Obs_GSL = site_pd.at[i, 'Obs-GPP_GSL']	
    EVI_SOS	 = site_pd.at[i, 'EVI-GPP_SOS']
    EVI_EOS	 = site_pd.at[i, 'EVI-GPP_EOS']
    EVI_GSL	 = site_pd.at[i, 'EVI-GPP_GSL']
    NIRv_SOS = site_pd.at[i, 'NIRv-GPP_SOS']	
    NIRv_EOS = site_pd.at[i, 'NIRv-GPP_EOS']	
    NIRv_GSL = site_pd.at[i, 'NIRv-GPP_GSL']
    
    if siteName.split("-")[0] == "US":
        region = "NorthAmerica"
    else:
        region = "Europe"       
        
    """
    五种数据的原始数据
    """
    ori_root = r"C:\Users\Administrator\Desktop\Paper"    
    gpp_ori_file = os.path.join(ori_root, "Output_GPP\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName))          
    obs_ori_file = os.path.join(ori_root, "Output_Obs\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)) 
    sc_ori_file = os.path.join(ori_root, "Output_SC\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName))   
    sr_ori_file = os.path.join(ori_root, "Output_SR\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)) 
    nir_ori_file = os.path.join(ori_root, "Output_NIRv\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)) 
    evi_ori_file = os.path.join(ori_root, "Output_EVI\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)) 

    """
    Beck拟合数据
    """ 
    fit_root = r"C:\Users\Administrator\Desktop\Paper" 
    gpp_fit_file = os.path.join(fit_root, "Fitting_GPP\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName))          
    obs_fit_file = os.path.join(fit_root, "Fitting_Obs\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName)) 
    sc_fit_file = os.path.join(fit_root, "Fitting_SC\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName))   
    sr_fit_file = os.path.join(fit_root, "Fitting_SR\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName)) 
    nir_fit_file = os.path.join(fit_root, "Fitting_NIRv\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName)) 
    evi_fit_file = os.path.join(fit_root, "Fitting_EVI\{0}\SitePixel\Polygon\DoubleLogBeck\{1}.csv".format(region,siteName)) 

    """
    Beck物候指标
    """    
    phe_root = r"C:\Users\Administrator\Desktop\Paper"
    gpp_phe_file = os.path.join(phe_root, "Pheno_GPP\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 
    obs_phe_file = os.path.join(phe_root, "Pheno_Obs\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 
    sc_phe_file = os.path.join(phe_root, "Pheno_SC\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 
    sr_phe_file = os.path.join(phe_root, "Pheno_SR\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 
    nir_phe_file = os.path.join(phe_root, "Pheno_NIRv\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 
    evi_phe_file = os.path.join(phe_root, "Pheno_EVI\{0}\SitePixel\Polygon\DoubleLogBeck_Trs\{1}.csv".format(region,siteName)) 

    """
    将原始数据读入数组
    """               
    gpp_ori_arr = np.loadtxt(gpp_ori_file, delimiter = ",")
    obs_ori_arr = np.loadtxt(obs_ori_file, delimiter = ",")
    sc_ori_arr = np.loadtxt(sc_ori_file, delimiter = ",")       
    sr_ori_arr = np.loadtxt(sr_ori_file, delimiter = ",")       
    nir_ori_arr = np.loadtxt(nir_ori_file, delimiter = ",")   
    evi_ori_arr = np.loadtxt(evi_ori_file, delimiter = ",")   

    """
    将Beck拟合数据读入数组
    """               
    gpp_fit_arr = np.loadtxt(gpp_fit_file, delimiter = ",")
    obs_fit_arr = np.loadtxt(obs_fit_file, delimiter = ",")
    sc_fit_arr = np.loadtxt(sc_fit_file, delimiter = ",")       
    sr_fit_arr = np.loadtxt(sr_fit_file, delimiter = ",")       
    nir_fit_arr = np.loadtxt(nir_fit_file, delimiter = ",")   
    evi_fit_arr = np.loadtxt(evi_fit_file, delimiter = ",")   
                                                                                                   
    """
    将物候指标读入数组
    """               
    gpp_phe_pd = pd.read_csv(gpp_phe_file, header = None)                                   
    obs_phe_pd = pd.read_csv(obs_phe_file, header = None) 
    sc_phe_pd = pd.read_csv(sc_phe_file, header = None)    
    sr_phe_pd = pd.read_csv(sr_phe_file, header = None) 
    nir_phe_pd = pd.read_csv(nir_phe_file, header = None)                                                                                                    
    evi_phe_pd = pd.read_csv(evi_phe_file, header = None)
    
    GPP_SOS, GPP_EOS = int(gpp_phe_pd.loc[0,1]), int(gpp_phe_pd.loc[1,1])
    Obs_SOS, Obs_EOS = int(obs_phe_pd.loc[0,1]), int(obs_phe_pd.loc[1,1])
    SC_SOS, SC_EOS = int(sc_phe_pd.loc[0,1]), int(sc_phe_pd.loc[1,1])
    SR_SOS, SR_EOS = int(sr_phe_pd.loc[0,1]), int(sr_phe_pd.loc[1,1]) 
    NIRv_SOS, NIRv_EOS = int(nir_phe_pd.loc[0,1]), int(nir_phe_pd.loc[1,1])
    EVI_SOS, EVI_EOS = int(evi_phe_pd.loc[0,1]), int(evi_phe_pd.loc[1,1])

    """
    计算逃逸概率fesc
    """  
    tropomi_ori_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Obs_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)
    tropomi_sc_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SC_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)
    tropomi_sr_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SR_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)   
    tropomi_cloud_18_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Cloud_2018\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)   

    tropomi_ori_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Obs_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)
    tropomi_sc_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SC_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)
    tropomi_sr_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_SR_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)   
    tropomi_cloud_19_file = r"F:\Chlorophyll_Fluorescence\Process\Europe\Step2.Merge\Version5\1Day\Output_Cloud_2019\{0}\SitePixel\Polygon\{1}.csv".format(region,siteName)   

    station_ori_arr = np.loadtxt(gpp_ori_file, delimiter = ",")
             
    tropomi_ori_18_arr = np.loadtxt(tropomi_ori_18_file, delimiter = ",")
    tropomi_sc_18_arr = np.loadtxt(tropomi_sc_18_file, delimiter = ",")       
    #tropomi_sr_18_arr = np.loadtxt(tropomi_sr_18_file, delimiter = ",")       
    tropomi_cloud_18_arr = np.loadtxt(tropomi_cloud_18_file, delimiter = ",")  

    tropomi_ori_19_arr = np.loadtxt(tropomi_ori_19_file, delimiter = ",")
    tropomi_sc_19_arr = np.loadtxt(tropomi_sc_19_file, delimiter = ",")       
    #tropomi_sr_19_arr = np.loadtxt(tropomi_sr_19_file, delimiter = ",")
    tropomi_cloud_19_arr = np.loadtxt(tropomi_cloud_19_file, delimiter = ",")  

    merge_18 = np.hstack((tropomi_ori_18_arr, tropomi_sc_18_arr[:,1].reshape(-1,1), tropomi_cloud_18_arr[:,1].reshape(-1,1)))                                                                                              
    merge_19 = np.hstack((tropomi_ori_19_arr, tropomi_sc_19_arr[:,1].reshape(-1,1), tropomi_cloud_19_arr[:,1].reshape(-1,1)))    
    
    merge_18 = np.delete(merge_18, np.where(merge_18[:,1] == -9999), axis = 0) 
    merge_18 = np.delete(merge_18, np.where(merge_18[:,2] == -9999), axis = 0) 
    merge_18 = np.delete(merge_18, np.where(merge_18[:,3] > 0.1), axis = 0) 
    
    merge_19 = np.delete(merge_19, np.where(merge_19[:,1] == -9999), axis = 0)
    merge_19 = np.delete(merge_19, np.where(merge_19[:,2] == -9999), axis = 0)  
    merge_19 = np.delete(merge_19, np.where(merge_19[:,3] > 0.1), axis = 0) 
   
    fesc_18 = np.hstack((merge_18, (merge_18[:,1]/merge_18[:,2]).reshape(-1,1)))
    fesc_19 = np.hstack((merge_19, (merge_19[:,1]/merge_19[:,2]).reshape(-1,1))) 

    fesc_18 = np.delete(fesc_18, np.where(abs(fesc_18[:,4]) > 1), axis = 0)    
    fesc_19 = np.delete(fesc_19, np.where(abs(fesc_19[:,4]) > 1), axis = 0)  

    fesc = np.vstack((fesc_18, fesc_19))
    fesc = np.vstack((fesc, np.array([[ 8,  1,  1, 1, fesc_18[0,4]],
                                      [360, 1,  1, 1, fesc_19[-1,4]]])))
    fesc_month = np.full((12,2), -9999.0)
    fesc_month = jd_to_time(fesc_18, fesc_19, fesc_month)

    fesc[:,4] = 1/fesc[:,4]
    fesc = np.delete(fesc, np.where(abs(fesc[:,4]) > 5), axis = 0)       
    fesc_month[:,1] = 1/fesc_month[:,1]
    fesc_month = np.delete(fesc_month, np.where(np.isnan(fesc_month[:,1])), axis = 0)   
    
    """
    将数据合并到固定的时间分辨率
    """                  
    gpp_ori_arr_merge = Merge_Ori(gpp_ori_arr, day_length_points)
    obs_ori_arr_merge = Merge_Ori(obs_ori_arr, day_length_points) 
    sc_ori_arr_merge = Merge_Ori(sc_ori_arr, day_length_points)     
    sr_ori_arr_merge = Merge_Ori(sr_ori_arr, day_length_points)     
    nir_ori_arr_merge = Merge_Ori(nir_ori_arr, day_length_points)
    evi_ori_arr_merge = Merge_Ori(evi_ori_arr, day_length_points)
     
    #Beck_     
    gpp_fit_arr_merge = Merge_Fit(gpp_fit_arr, day_length)
    obs_fit_arr_merge = Merge_Fit(obs_fit_arr, day_length)    
    sc_fit_arr_merge = Merge_Fit(sc_fit_arr, day_length)        
    sr_fit_arr_merge = Merge_Fit(sr_fit_arr, day_length)        
    nir_fit_arr_merge = Merge_Fit(nir_fit_arr, day_length)
    evi_fit_arr_merge = Merge_Fit(evi_fit_arr, day_length)
    
    """
    数据标准化
    """  
    #Original
    x_ori_gpp, y_ori_gpp = Normalization(gpp_ori_arr_merge)
    x_ori_obs, y_ori_obs = Normalization(obs_ori_arr_merge)
    x_ori_sc, y_ori_sc = Normalization(sc_ori_arr_merge)    
    x_ori_sr, y_ori_sr = Normalization(sr_ori_arr_merge)    
    x_ori_nir, y_ori_nir = Normalization(nir_ori_arr_merge)
    x_ori_evi, y_ori_evi = Normalization(evi_ori_arr_merge)

    #Beck_     
    x_fit_gpp, y_fit_gpp = Normalization(gpp_fit_arr_merge)
    x_fit_obs, y_fit_obs = Normalization(obs_fit_arr_merge)
    x_fit_sc, y_fit_sc = Normalization(sc_fit_arr_merge)
    x_fit_sr, y_fit_sr = Normalization(sr_fit_arr_merge)
    x_fit_nir, y_fit_nir = Normalization(nir_fit_arr_merge)
    x_fit_evi, y_fit_evi = Normalization(evi_fit_arr_merge)

    
    """
    将物候指标读入数组
    """
    """
    x_GPP = [GPP_SOS, GPP_EOS]
    y_GPP = [y_fit_gpp[GPP_SOS - 1], y_fit_gpp[GPP_EOS - 1]]
    x_Obs = [Obs_SOS, Obs_EOS]
    y_Obs = [y_fit_obs[Obs_SOS - 1], y_fit_obs[Obs_EOS - 1]]
    x_SC = [SC_SOS, SC_EOS]
    y_SC = [y_fit_sc[SC_SOS - 1], y_fit_sc[SC_EOS - 1]]
    x_SR = [SR_SOS, SR_EOS]
    y_SR = [y_fit_sr[SR_SOS - 1], y_fit_sr[SR_EOS - 1]]
    x_NIR = [NIRv_SOS, NIRv_EOS]
    y_NIR = [y_fit_nir[NIRv_SOS - 1], y_fit_nir[NIRv_EOS - 1]]
    x_EVI = [EVI_SOS, EVI_EOS]
    y_EVI = [y_fit_evi[EVI_SOS - 1], y_fit_evi[EVI_EOS - 1]]  
    """  
    """
    分别对原始值和拟合值进行求相关系数
    """    
    y_obs_sta_ori = np.hstack((y_ori_gpp.reshape(-1,1), y_ori_obs.reshape(-1,1)))
    y_sc_sta_ori = np.hstack((y_ori_gpp.reshape(-1,1), y_ori_sc.reshape(-1,1)))    
    y_sr_sta_ori = np.hstack((y_ori_gpp.reshape(-1,1), y_ori_sr.reshape(-1,1)))    
    y_nir_sta_ori = np.hstack((y_ori_gpp.reshape(-1,1), y_ori_nir.reshape(-1,1)))
    y_evi_sta_ori = np.hstack((y_ori_gpp.reshape(-1,1), y_ori_evi.reshape(-1,1)))

    #Beck_     
    y_obs_sta_fit = np.hstack((y_fit_gpp.reshape(-1,1), y_fit_obs.reshape(-1,1)))
    y_sc_sta_fit = np.hstack((y_fit_gpp.reshape(-1,1), y_fit_sc.reshape(-1,1)))   
    y_sr_sta_fit = np.hstack((y_fit_gpp.reshape(-1,1), y_fit_sr.reshape(-1,1)))   
    y_nir_sta_fit = np.hstack((y_fit_gpp.reshape(-1,1), y_fit_nir.reshape(-1,1)))
    y_evi_sta_fit = np.hstack((y_fit_gpp.reshape(-1,1), y_fit_evi.reshape(-1,1)))

    y_obs_sta_ori = np.delete(y_obs_sta_ori, np.where(np.isnan(y_obs_sta_ori)), axis = 0)   #TROPOMI 
    y_sc_sta_ori = np.delete(y_sc_sta_ori, np.where(np.isnan(y_sc_sta_ori)), axis = 0)   #TROPOMI     
    y_sr_sta_ori = np.delete(y_sr_sta_ori, np.where(np.isnan(y_sr_sta_ori)), axis = 0)   #TROPOMI     
    y_nir_sta_ori = np.delete(y_nir_sta_ori, np.where(np.isnan(y_nir_sta_ori)), axis = 0)   #TROPOMI    
    y_evi_sta_ori = np.delete(y_evi_sta_ori, np.where(np.isnan(y_evi_sta_ori)), axis = 0)   #TROPOMI    

    #Beck_
    y_obs_sta_fit = np.delete(y_obs_sta_fit, np.where(np.isnan(y_obs_sta_fit)), axis = 0)   #TROPOMI 
    y_sc_sta_fit = np.delete(y_sc_sta_fit, np.where(np.isnan(y_sc_sta_fit)), axis = 0)   #TROPOMI 
    y_sr_sta_fit = np.delete(y_sr_sta_fit, np.where(np.isnan(y_sr_sta_fit)), axis = 0)   #TROPOMI 
    y_nir_sta_fit = np.delete(y_nir_sta_fit, np.where(np.isnan(y_nir_sta_fit)), axis = 0)   #TROPOMI    
    y_evi_sta_fit = np.delete(y_evi_sta_fit, np.where(np.isnan(y_evi_sta_fit)), axis = 0)   #TROPOMI    
           
    R_obs_ori = stats.pearsonr(y_obs_sta_ori[:,0], y_obs_sta_ori[:,1])[0]
    R_sc_ori = stats.pearsonr(y_sc_sta_ori[:,0], y_sc_sta_ori[:,1])[0]   
    R_sr_ori = stats.pearsonr(y_sr_sta_ori[:,0], y_sr_sta_ori[:,1])[0]   
    R_nir_ori = stats.pearsonr(y_nir_sta_ori[:,0], y_nir_sta_ori[:,1])[0]
    R_evi_ori = stats.pearsonr(y_evi_sta_ori[:,0], y_evi_sta_ori[:,1])[0]

    #Beck_ 
    import math
    n = 2
    R_obs_fit = math.floor(stats.pearsonr(y_obs_sta_fit[:,0], y_obs_sta_fit[:,1])[0]*10**n)/(10**n)
    R_sc_fit = math.floor(stats.pearsonr(y_sc_sta_fit[:,0], y_sc_sta_fit[:,1])[0]*10**n)/(10**n)
    R_sr_fit = math.floor(stats.pearsonr(y_sr_sta_fit[:,0], y_sr_sta_fit[:,1])[0]*10**n)/(10**n)
    R_nir_fit = stats.pearsonr(y_nir_sta_fit[:,0], y_nir_sta_fit[:,1])[0]
    R_evi_fit = stats.pearsonr(y_evi_sta_fit[:,0], y_evi_sta_fit[:,1])[0]
    site = i
    row = int(site / ncols) * 2
    col = int(site % ncols)
    
    axes[row,col].set_title("({0}) {1} {2}".format(chr(97+site), siteName, siteType), fontsize = ftsize, family="Times New Roman")  
    """
    标记物候指标
    """
    """
    axes[0,site].vlines(x_SC, 0, 1, colors = TROPOMI_Total_Color,linestyles = 'dashed', linewidths = vline_width)
    axes[0,site].vlines(x_SR, 0, 1, colors = TROPOMI_Soil_Color,linestyles = 'dashed', linewidths = vline_width)
    axes[0,site].vlines(x_GPP, 0, 1, colors = GPP_Color,linestyles = 'solid', linewidths = vline_width)
    """
   
    GPP, = axes[row,col].plot(x_fit_gpp, y_fit_gpp, GPP_Color, label = 'GPP', lw = linewidth)  
    Obs, = axes[row,col].plot(x_fit_obs, y_fit_obs, Obs_Color, label = 'TROPOMI_Original', lw = linewidth)     
    SC, = axes[row,col].plot(x_fit_sc, y_fit_sc, SC_Color,label = 'TROPOMI_Total', lw = linewidth) 
    #SR, = axes[0,site].plot(x_fit_sr, y_fit_sr, SR_Color,label = 'TROPOMI_Soil', lw = linewidth)  
    """
    axes[0,site].scatter(x_ori_gpp, y_ori_gpp, c = GPP_Color, marker = 'o', s = pointsize)   
    axes[0,site].scatter(x_ori_obs, y_ori_obs, c = Obs_Color, marker = 'x', s = pointsize) 
    axes[0,site].scatter(x_ori_sc, y_ori_sc, c = SC_Color, marker = '+', s = pointsize)  
    axes[0,site].scatter(x_ori_sr, y_ori_sr, c = SR_Color, marker = '*', s = pointsize)
    """
    axes[row,col].spines['left'].set_linewidth(0.7)
    axes[row,col].spines['right'].set_linewidth(0.7)
    axes[row,col].spines['top'].set_linewidth(0.7)
    axes[row,col].spines['bottom'].set_linewidth(0.7)
    axes[row,col].tick_params(axis='both', length = 1.5, width = 0.8, labelsize = 4.5)     
    #axes[0].legend(framealpha  = 0)    
    axes[row,col].set_xticks(np.linspace(0,395,13)[:-1])
    axes[row,col].set_xticklabels(('Jan', 'Feb', 'Mar','Apr','May','Jun','Jul','Aug','Sep', 'Oct','Nov', 'Dec'), family="Times New Roman") 
    if site%ncols == 0:  
        axes[row,col].set_ylabel('Normalized Values (days)', fontsize = ftsize/1.2, family="Times New Roman")
    axes[row,col].set_yticks(np.linspace(0,1.2,7)[:-1])
    axes[row,col].set_yticklabels([0.0,0.2,0.4,0.6,0.8,1.0], fontsize = labelsize, family="Times New Roman")                

    #axes[0,site].annotate('SOS', (x_TRO_SR[0] - 50, 0.5), color='green',size=7, family="Times New Roman")
    #axes[0,site].annotate('EOS', (x_TRO_SR[1] + 20, 0.5), color='red',size=7, family="Times New Roman")     
     #R$\mathregular{^2}$
    axes[row,col].annotate('MAE:\nSIF$_\mathdefault{SC}$\n' + ' = {0},{1}'.format(int(abs(TRO_SC_SOS)), int(abs(TRO_SC_EOS))),\
                     (0.015, 0.75), textcoords='axes fraction', color = GPP_Color, size=annotate_size, family="Times New Roman")    
    axes[row,col].annotate('SIF$_\mathdefault{Obs}$\n' + ' = {0},{1}'.format(int(abs(TRO_Obs_SOS)), int(abs(TRO_Obs_EOS))),\
                     (0.015, 0.55), textcoords='axes fraction', color = GPP_Color, size=annotate_size, family="Times New Roman")
        
    #axes[row,col].annotate('\u0394MAE$_G$$_S$$_L$ = {0}'.format(int(abs(TRO_SC_GSL)-abs(TRO_Obs_GSL))),\
    #                 (0.015, 0.65), textcoords='axes fraction', color = GPP_Color, size=annotate_size, family="Times New Roman")            
      
    

    
    """
    标记物候指标
    """     
    axes[row + 1,col].scatter(fesc[:,0], fesc[:,4],  color='black', marker='o',  s = pointsize) # 把 corlor 设置为空，通过edgecolors来控制颜色 
    axes[row + 1,col].plot(fesc_month[:,0], fesc_month[:,1],  color='red', lw = linewidth) # 把 corlor 设置为空，通过edgecolors来控制颜色 

    # rft = np.fft.rfft(y)  
    axes[row + 1,col].spines['left'].set_linewidth(0.7)
    axes[row + 1,col].spines['right'].set_linewidth(0.7)
    axes[row + 1,col].spines['top'].set_linewidth(0.7)
    axes[row + 1,col].spines['bottom'].set_linewidth(0.7)
    axes[row + 1,col].tick_params(axis='both', length = 1.5, width = 0.8, labelsize = 4.5)    
    #axes[0].legend(framealpha  = 0)    
    axes[row + 1,col].set_xticks(np.linspace(0,395,13)[:-1])
    axes[row + 1,col].set_xticklabels(('Jan', 'Feb', 'Mar','Apr','May','Jun','Jul','Aug','Sep', 'Oct','Nov', 'Dec'), family="Times New Roman") 
    if site%ncols == 0:    
        axes[row + 1,col].set_ylabel('1/Fesc', fontsize = ftsize/1.2, family="Times New Roman")
    axes[row + 1,col].set_yticks(np.linspace(1.0,6.0,6)[:-1])
    axes[row + 1,col].set_yticklabels([1.0,2.0,3.0,4.0,5.0], fontsize = labelsize, family="Times New Roman")                    
         
handles = [GPP, SC, Obs]
labels = ['GPP$_\mathdefault{EC}$','SIF$_\mathdefault{total-SC}$',\
          #'SIF$_T$$_R$$_O$$_P$$_O$$_M$$_I$$_\_$$_t$$_o$$_t$$_a$$_l$$_-$$_S$$_R$',\
          'SIF$_\mathdefault{Obs}$']  

fig.tight_layout() 
fig.subplots_adjust(bottom=0.12)
fig.legend(handles, labels, loc ='lower center', bbox_to_anchor=(0.5, 0.04),
          fancybox = False, shadow = False,frameon = False, ncol = 4, 
          handletextpad = 0.2, columnspacing = 1.2, prop={'family':"Times New Roman", 'size':12})                                 
Plot_path = os.path.join(plot_path, "Manuscript5.jpg")
plt.show()
fig.savefig(Plot_path, dpi=600, quality=100)

