library(ggplot2)
library(grid)
library(gridExtra)
library(raster)
library(rasterVis)
library(sf)
library(RColorBrewer)

#Plot in Manuscript, Using TROPOMI SC and EVI data 
SCEVI_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-EVI_1.tif")
SCEVI_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-EVI_3.tif")
SCEVI_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-EVI_5.tif")

SCNIR_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_1.tif")
SCNIR_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_3.tif")
SCNIR_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_5.tif")

shp_file = "F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.0.BasicData/North.shp"
# read the shapefile
world_shp <- read_sf(shp_file)
world_outline <- as(st_geometry(world_shp), Class="Spatial")

# Create a time series raster stack
raster_stack = stack(SCEVI_Delta_SOS,SCNIR_Delta_SOS,
                     SCEVI_Delta_EOS,SCNIR_Delta_EOS,
                     SCEVI_Delta_LOS,SCNIR_Delta_LOS)


#"dark green","red2",
ColorKey_Height = 0.9
ColorKey_Width = 0.7
ColorKey_tck = 0.5
Aixs_tck = 0.3
Aixs_size = 0.7
Main_Unit_x = 0.55 #调整标题到图像左右距离
Main_Unit_y = 0.6 #调整标题到图像上下距离

#cols = colorRampPalette(c("blue4","dodgerblue4","dodgerblue3","tomato1","red2","red4"))
#my.at = seq(-100, 100, length = 13)
cols = colorRampPalette(c(brewer.pal(9, "Reds")[8:4],brewer.pal(9, "Blues")[5:9]))
my.at = seq(-60, 60, length = 500)

#子标题的参数
p.strip = list(cex=1.2, fontfamily='serif')

l1 = levelplot(raster_stack, margin = FALSE, layout = c(2, 3), 
                cex = Main_size, at=my.at,
                col.regions = cols, xlab = NULL, ylab = NULL,
                names.attr = c(expression((a) ~Delta~SOS (SIF[total-SC]-EVI)),
                               expression((b) ~Delta~SOS (SIF[total-SC]-NIR[V])),                               
                               expression((c) ~Delta~EOS (SIF[total-SC]-EVI)),
                               expression((d) ~Delta~EOS (SIF[total-SC]-NIR[V])),                               
                               expression((e) ~Delta~GSL (SIF[total-SC]-EVI)),
                               expression((f) ~Delta~GSL (SIF[total-SC]-NIR[V]))
                               ),
                par.strip.text = p.strip,
                scales=list(x = list(alternating=3), y = list(alternating=3), 
                            tck = Aixs_tck, cex = Aixs_size),
                colorkey = list(space="right",tck = ColorKey_tck,
                                height = ColorKey_Height, width = ColorKey_Width,
                                at = my.at))
l1 = l1 + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.6))
l1$par.settings$axis.text$fontfamily = "serif"
plot(l1)
grid_plot = grid.arrange(l1, nrow=1, ncol=1) #use package gridExtra
grid.rect(gp = gpar(fill=NA, fontfamily="serif"))
ggsave("F:/Chlorophyll_Fluorescence/Process/Europe/Step13.Phenology_Raster/13.6.Plot_Delta/Rplot-EVI2.jpg",
        units = "in", dpi = 600, plot = grid_plot)

dev.off()