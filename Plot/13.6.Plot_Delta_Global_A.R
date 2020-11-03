SCNIR_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_1.tif")
SCNIR_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_3.tif")
SCNIR_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SC-NIRv_5.tif")

SREVI_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-EVI_1.tif")
SREVI_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-EVI_3.tif")
SREVI_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-EVI_5.tif")

SRNIR_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-NIRv_1.tif")
SRNIR_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-NIRv_3.tif")
SRNIR_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/SR-NIRv_5.tif")

dcEVI_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-EVI_1.tif")
dcEVI_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-EVI_3.tif")
dcEVI_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-EVI_5.tif")

dcNIR_Delta_SOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-NIRv_1.tif")
dcNIR_Delta_EOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-NIRv_3.tif")
dcNIR_Delta_LOS = raster("F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.2.Delta_Plot/dcSIF-NIRv_5.tif")

shp_file = "F:/Chlorophyll_Fluorescence/Process/Europe/Step14.Boreal_Map/14.0.BasicData/North.shp"
# read the shapefile
world_shp <- read_sf(shp_file)
world_outline <- as(st_geometry(world_shp), Class="Spatial")

# Create a time series raster stack
raster_stack = stack(SCNIR_Delta_SOS, SCNIR_Delta_EOS, SCNIR_Delta_LOS,
                     SREVI_Delta_SOS, SREVI_Delta_EOS, SREVI_Delta_LOS,
                     SRNIR_Delta_SOS, SRNIR_Delta_EOS, SRNIR_Delta_LOS,
                     dcEVI_Delta_SOS, dcEVI_Delta_EOS, dcEVI_Delta_LOS,
                     dcNIR_Delta_SOS, dcNIR_Delta_EOS, dcNIR_Delta_LOS)

cols = colorRampPalette(c("midnightblue","dodgerblue4","skyblue", "white","red"))
ColorKey_Height = 0.8
ColorKey_Width = 0.7
ColorKey_tck = 0
Aixs_tck = 0.3
Aixs_size = 0.7
Main_size = 0 #调整字体大小
Main_Unit_x = 0.55 #调整标题到图像左右距离
Main_Unit_y = 0.5 #调整标题到图像上下距离
my.at = seq(-150, 50, length = 500)
#子标题的参数
p.strip = list(cex=0.8, fontfamily='serif')

l1 = levelplot(raster_stack, margin = FALSE, layout = c(3, 5), 
                cex = Main_size, at=my.at,
                col.regions = cols, xlab = NULL, ylab = NULL,
                names.attr = c(expression((a) ~Delta~SOS (SIF[total-SC] - NIR[V])),
                               expression((b) ~Delta~EOS (SIF[total-SC] - NIR[V])),
                               expression((c) ~Delta~GSL (SIF[total-SC] - NIR[V])),
                               
                               expression((d) ~Delta~SOS (SIF[total-SR] - EVI)),
                               expression((e) ~Delta~EOS (SIF[total-SR] - EVI)),
                               expression((f) ~Delta~GSL (SIF[total-SR] - EVI)), 
                               
                               expression((g) ~Delta~SOS (SIF[total-SR] - NIR[V])),
                               expression((h) ~Delta~EOS (SIF[total-SR] - NIR[V])),
                               expression((i) ~Delta~GSL (SIF[total-SR] - NIR[V])),

                               expression((j) ~Delta~SOS (SIF[Obs] - EVI)),
                               expression((k) ~Delta~EOS (SIF[Obs] - EVI)),
                               expression((l) ~Delta~GSL (SIF[Obs] - EVI)),
                               
                               expression((m) ~Delta~SOS (SIF[Obs] - NIR[V])),
                               expression((n) ~Delta~EOS (SIF[Obs] - NIR[V])),
                               expression((o) ~Delta~GSL (SIF[Obs] - NIR[V]))                              
                               ),
                par.strip.text = p.strip,
                scales=list(x = list(alternating=3), y = list(alternating=3), 
                            tck = Aixs_tck, cex = Aixs_size),
                colorkey = list(space="right",tck = ColorKey_tck,
                                height = ColorKey_Height, width = ColorKey_Width,
                                at = my.at))
l1 = l1 + latticeExtra::layer(sp.lines(world_outline, col="black", lwd=0.6))
l1$par.settings$axis.text$fontfamily = "serif"

grid_plot = grid.arrange(l1, nrow=1, ncol=1) #use package gridExtra
grid.rect(gp = gpar(fill=NA))
ggsave("F:/Chlorophyll_Fluorescence/Process/Europe/Step13.Phenology_Raster/13.6.Plot_Delta/Rplot_Merge.jpg",
        units = "in", dpi = 600, plot = grid_plot)

#dev.off()