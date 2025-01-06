install.packages("spgwr")
install.packages("ggspatial")
install.packages("GISTools")
library(GISTools)
library(rgdal)
library(rasterVis) 
library(stars) 
library(raster)
library(biodivMapR)
library(labdsv)
library(tools)
library(ggplot2)
library(gridExtra)
library(stars)
library(ggspatial)


data_directory1 <- "/user/ziqitang/Documents/EMIT_L2A_RFL_001_20230119T114235_2301907_004_reflectance"
NameRaster1 <- 'EMIT_L2A_RFL_001_20230119T114235_2301907_004_reflectance'
raster_file1 <- file.path(data_directory1, NameRaster1)
hyperspectral_image <- brick(raster_file1)
print(hyperspectral_image)

library(raster)

shannon_index <- function(pixel_values) {
  pixel_values <- pixel_values / sum(pixel_values, na.rm = TRUE) 
  pixel_values <- pixel_values[pixel_values > 0] 
  -sum(pixel_values * log(pixel_values))
}

shannon_diversity <- calc(hyperspectral_image, fun = shannon_index)

png("shannon_diversity.png", width = 800, height = 600)
plot(shannon_diversity, main = "Shannon Index for Spectral Diversity")
dev.off()









hdr_file1 <- get_HDR_name(raster_file1, showWarnings = FALSE)
read_ENVI_header(hdr_file1)
raster_data1 <- read_stars(raster_file1)
Input_Image_File1 <- raster_data1

Input_Mask_File <- FALSE
# Define path for master output directory where files produced during the process are saved
Output_Dir <- '/user/ziqitang/Documents/tzq/output'
dir.create(path = Output_Dir,recursive = T,showWarnings = F)
# Define levels for radiometric filtering
NDVI_Thresh <- 0.02
Blue_Thresh <- 500
NIR_Thresh <- 1400
Continuum_Removal <- TRUE
TypePCA <- 'SPCA'
FilterPCA <- FALSE
window_size <- 10
nbCPU <- 4
MaxRAM <- 0.2
nbclusters <- 50

print("PERFORM RADIOMETRIC FILTERING")
Input_Mask_File <- perform_radiometric_filtering(Image_Path = raster_file1,
                                                 Mask_Path = Input_Mask_File,
                                                 Output_Dir = Output_Dir,
                                                 TypePCA = TypePCA,
                                                 NDVI_Thresh = NDVI_Thresh,
                                                 Blue_Thresh = Blue_Thresh,
                                                 NIR_Thresh = NIR_Thresh)




print("PERFORM DIMENSIONALITY REDUCTION")
Excluded_WL <- c(0, 400)
Excluded_WL <- rbind(Excluded_WL, c(895, 1005))
Excluded_WL <- rbind(Excluded_WL, c(1180, 1480))
Excluded_WL <- rbind(Excluded_WL, c(1780, 2040))


PCA_Output <- perform_PCA(Input_Image_File = raster_file1,
                          Input_Mask_File = FALSE,
                          Output_Dir = Output_Dir,
                          TypePCA = TypePCA,
                          FilterPCA = FilterPCA,
                          nbCPU = nbCPU,
                          MaxRAM = MaxRAM,
                          Continuum_Removal = Continuum_Removal)
PCA_Files <- PCA_Output$PCA_Files
# path for the updated mask
Input_Mask_File <- PCA_Output$MaskPath

data_directory2 <- "/user/ziqitang/Documents/tzq/EMIT_L2A_RFL_001_20230119T114235_2301907_004_reflectance/SPCA/PCA"
NameRaster2 <- 'OutputPCA_30_PCs'
raster_file2 <- file.path(data_directory2, NameRaster2)
hyperspectral_image2 <- brick(raster_file2)
plot(hyperspectral_image2)

print("MAP SPECTRAL SPECIES")
Kmeans_info <- map_spectral_species(Input_Image_File = raster_file1, 
                                    Input_Mask_File = PCA_Output$MaskPath,
                                    Output_Dir = Output_Dir,
                                    SpectralSpace_Output = PCA_Output, 
                                    nbclusters = nbclusters, 
                                    nbCPU = nbCPU, MaxRAM = MaxRAM)
print("MAP ALPHA DIVERSITY")
# Index.Alpha   = c('Shannon','Simpson')
Index_Alpha <- c('Shannon')
map_alpha_div(Input_Image_File = raster_file1, 
              Output_Dir = Output_Dir, 
              TypePCA = TypePCA,
              window_size = window_size,
              nbCPU = nbCPU,
              MaxRAM = MaxRAM,
              Index_Alpha = Index_Alpha, 
              nbclusters = nbclusters)

print("MAP BETA DIVERSITY")
map_beta_div(Input_Image_File = raster_file1, 
             Output_Dir = Output_Dir, 
             TypePCA = TypePCA,
             window_size = window_size, 
             nbCPU = nbCPU, 
             MaxRAM = MaxRAM,
             nbclusters = nbclusters)



library(terra)
raster_image1 <- rast('/user/ziqitang/Documents/tzq/output/EMIT_L2A_RFL_001_20230119T114235_2301907_004_reflectance/SPCA/BETA/BetaDiversity_BCdiss_PCO_10')
plot(raster_image1,main="Distribution map of BetaDiversity of the area")

raster_SpectralSpecies <- rast('/user/ziqitang/Documents/tzq/output/EMIT_L2A_RFL_001_20230119T114235_2301907_004_reflectance/SPCA/SpectralSpecies/SpectralSpecies_Distribution')
plot(raster_SpectralSpecies,main="SpectralSpecies_Distribution")




landcover_data <- rast("/user/ziqitang/Documents/SA_NLC_2022_GEO.tif/SA_NLC_2022_GEO.tif")
landcover_data_projected <- terra::project(landcover_data, raster_image1)
lon_min <- 18.3
lon_max <- 18.6
lat_min <- -34.4
lat_max <- -33.8
bbox_extent <- ext(lon_min, lon_max, lat_min, lat_max)
cropped_raster1 <- crop(raster_image1, bbox_extent)
plot(cropped_raster1,main='Spatial Distribution of Shannon Index (Cropped Area)')

# x1 - land cover type
landcover_data_cropped <- crop(landcover_data_projected,bbox_extent)
plot(landcover_data_cropped,main='landcover_data_cropped')
filtered_layer <- subst(landcover_data_cropped, "natural ocean & coastal", NA)
plot(filtered_layer, main = " Distribution map of filtered: Natural Ocean & Coastal")

# x2 - fire age
image_fire_age <- rast("/user/ziqitang/Documents/tzq/years_since_fire.tif")
plot(image_fire_age, main="Distribution map of fire age")
image_fire_age_clean <- classify(image_fire_age, cbind(NA, 0))
plot(image_fire_age_clean)
image_fire_age_projected <- terra::project(image_fire_age_clean, filtered_layer)
plot(image_fire_age_projected)
fire_age_cropped <- crop(image_fire_age_projected, filtered_layer)
plot(fire_age_cropped,main='Spatial Distribution of Fire Age (Cropped Area)')
#fire_age_resampled <- resample(fire_age_cropped,cropped_raster2)
#fire_age_intersection <- mask(fire_age_resampled,cropped_raster2)
#plot(fire_age_intersection,main='Distribution of fire age in an area')
# cropped_fire_age_intersection <- crop(fire_age_intersection, bbox_extent)
# plot(cropped_fire_age_intersection,main='cropped raster of fire age')


library(spgwr)
library(raster)

fire_age_cropped <- raster(fire_age_cropped)
points_raster1 <- rasterToPoints(fire_age_cropped, spatial = TRUE)  # independent
points_raster2 <- rasterToPoints(raster(cropped_raster1), spatial = TRUE)  # dependent response y

points_raster2_df <- as.data.frame(points_raster2[["Shannon_10"]])
points_raster1_df <- as.data.frame(points_raster1)

data <- cbind(points_raster1_df, points_raster2_df)
head(data)
#data <- cbind(points_raster1, points_raster2[["lyr.1"]])  
#names(data) <- c("x", "y", "response", "predictor") 
coords <- coordinates(points_raster1)
response <- points_raster1@data[, 1]  
predictor <- points_raster2@data[, "Shannon_10"]
data <- data.frame(x = coords[, 1], y = coords[, 2], response = response, predictor = predictor)
head(data)


bandwidth <- gwr.sel(response ~ predictor, data = data, coords = cbind(data$x, data$y))
gwr_model <- gwr(response ~ predictor, data = data, bandwidth = bandwidth, coords = cbind(data$x, data$y))
summary(gwr_model)








