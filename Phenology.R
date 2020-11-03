#关闭警告
options(warn = -1)

#计算物候期的两种方法
phenology = function(values, file, outFile, fitting_type){
  #阈值法计算物候期
  Trs = PhenoTrs(values, approach = "White", calc.pheno = TRUE, plot = FALSE)
  write.table(Trs, file = outFile, sep =",", row.names =TRUE, col.names =FALSE, quote =FALSE)
  
}

fitting = function(file, outFile){
  VI_arr = read.table(file, sep=',', header = FALSE)
  arr_day = as.matrix(VI_arr[1])
  arr = as.matrix(VI_arr[2])
  #将数据合并成Monthly
  arr_merge = matrix(0,12)
  for(i in 1:12){
    Month = i
    if (Month != 12){
      Start_Day = yday(as.Date(paste('2018-',Month,'-01', sep =""), format="%Y-%m-%d"))
      End_Day = yday(as.Date(paste('2018-',(Month + 1),'-01', sep =""), format="%Y-%m-%d"))      
    }
    if (Month == 12){
      Start_Day = yday(as.Date(paste('2018-',Month,'-01', sep =""), format="%Y-%m-%d"))
      End_Day = 365     
    }
    Value = 0
    Count = 0
    for(j in 1: length(arr_day)){
      if((floor(arr_day[j]) >= Start_Day)&&(floor(arr_day[j]) < End_Day)){
        if(floor(arr[j]) != -9999){
          Value = Value + arr[j]
          Count = Count + 1
        }
      }
    }
    if (Count != 0){
      arr_merge[i] = Value/Count      
    }
    else{
      arr_merge[i] = NA      
    }
  }

  #gapfilling process  
  arr_mulitiple = rbind(arr_merge, arr_merge, arr_merge)
  gaps = ts(data = arr_mulitiple, freq = 12, start = c(2018, 1))
  fill = FillPermanentGaps(gaps)
  
  #extract yearly data
  x = as.vector(window(fill, start=c(2018,1), end=c(2018, 12)))    
  # time steps for output (daily)
  tout <- seq(1, 12, length=365)	# time steps for output (daily)
  
  #DoubleLogBeck滤波平滑时间序列
  fit_Beck = FitDoubleLogBeck(x, tout = tout)
  phenology(fit_Beck$predicted, file, outFile, "DoubleLogBeck")
  
}

root = 'C:/Users/Administrator/Desktop/'

Datasets = c('SC','SR','Obs','EVI','NIRv','GPP')
Regions = c('Europe','NorthAmerica')
for (Dataset in Datasets){
  for (Region in Regions){
    inPath = paste(root,"Output_", Dataset,"/",Region,"/SitePixel/Polygon/" ,sep ="") 
    outPath = paste(root,"Pheno_",Dataset,"/",Region,"/SitePixel/Polygon/" ,sep ="")       
    fileLists = Sys.glob(paste(inPath, "*.csv", sep = ""))
    for (file in fileLists){
      outPath_Beck = paste(outPath, "/","DoubleLogBeck_Trs", sep ="")
      if(dir.exists(outPath_Beck) == FALSE){
        dir.create(outPath_Beck)
      }
      outFile = paste(outPath_Beck, "/",strsplit(file,"/")[[1]][9], sep ="")
      if (file.exists(outFile)){
        file.remove(outFile)
      }
      fitting(file, outFile)
      print(file)
    }      
  }
}

