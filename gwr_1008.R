install.packages("spgwr")
install.packages("terra")
install.packages("car")
install.packages("GWmodel")
install.packages("spgwr") 
install.packages("raster") 
install.packages("car")   
install.packages("dplyr") 
install.packages("leaflet")
install.packages("magrittr") 
install.packages("leaflet")  
install.packages("mapview")
install.packages("webshot2")
install.packages("spdep")
install.packages("terra") 
install.packages("randomForest")
install.packages("ggsci")
install.packages("httr")
install.packages("ggspatial")
install.packages("ggmap")
install.packages("GGally")
library(GGally)
library(httr)
library(randomForest)
library(ggsci)
library(spdep)
library(terra)
library(webshot2)
library(magrittr)
library(leaflet)
library(dplyr)
library(car)
library(raster)
library(spgwr)
library(caret)
library(sp)
library(mapview)
library(RColorBrewer)
library(ggspatial)
library(sf)
library(ggplot2)
library(tmap)
library(viridis)
library(ggmap)
library(tmap)
library(gridExtra)

#data
dependent_raster <- rast("/user/ziqitang/data_ziqi/clipped_shannon.tif")
dependent_clean <- classify(dependent_raster, cbind(0, NA))
plot(dependent_clean, main = "Alpha Diversity")
total_cells <- ncell(dependent_raster)
cat("Total number of cells:", total_cells, "\n")

#
dependent_clean <- project(dependent_clean, "EPSG:4326")
dependent_df <- as.data.frame(dependent_clean, xy = TRUE, na.rm = TRUE)
colnames(dependent_df) <- c("x", "y", "alpha_diversity") 
res(dependent_clean)
# Reproject the raster to a UTM projection (adjust EPSG code as needed)
dependent_clean_proj <- project(dependent_clean, "EPSG:32633")
res(dependent_clean_proj)
cell_area <- prod(res(dependent_clean_proj)) 
non_na_cells <- sum(!is.na(values(dependent_clean_proj)))
total_area <- non_na_cells * cell_area  # Area in m²
total_area_km2 <- total_area / 1e6      # Area in km²
cat("Total Area:\n")
cat("In square meters:", total_area, "\n")
cat("In square kilometers:", total_area_km2, "\n")


# Load and clean the raster
dependent_raster <- rast("/user/ziqitang/data_ziqi/clipped_shannon.tif")
dependent_clean <- classify(dependent_raster, cbind(0, NA))
plot(dependent_clean,main='dependent_clean')
alpha_df <- as.data.frame(dependent_raster, xy = TRUE, na.rm = TRUE)
names(alpha_df)[3] <- "layer" 
alpha_df$layer <- as.numeric(alpha_df$layer)
coordinates(alpha_df) <- ~x + y  
alpha_nb <- knn2nb(knearneigh(coordinates(alpha_df), k = 4))
alpha_listw <- nb2listw(alpha_nb, style = "W")
alpha_moran_result <- moran.test(alpha_df$layer, alpha_listw)
print(alpha_moran_result)
crs_dependent <- crs(dependent_clean)

#ndvi independent variable1
independent_raster0 <- rast("/user/ziqitang/data_ziqi/true_ndvi.tif") 
independent_raster0 <- project(independent_raster0, crs(dependent_clean))
independent0_resampled <- resample(independent_raster0, dependent_clean, method = "bilinear")
combined_mask0 <- mask(independent0_resampled, dependent_clean)
combined_independent0 <- mask(independent0_resampled, combined_mask0)
plot(combined_independent0, main = 'NDVI')
ndvi_df <- as.data.frame(combined_independent0, xy = TRUE, na.rm = TRUE)
names(ndvi_df)[3] <- "layer" 
ndvi_df$layer <- as.numeric(ndvi_df$layer)
coordinates(ndvi_df) <- ~x + y  
ndvi_nb <- knn2nb(knearneigh(coordinates(ndvi_df), k = 4))
ndvi_listw <- nb2listw(ndvi_nb, style = "W")
ndvi_moran_result <- moran.test(ndvi_df$layer, ndvi_listw)
print(ndvi_moran_result)

# convert NDVI raster to a data frame
ndvi_df <- as.data.frame(combined_independent0, xy = TRUE, na.rm = TRUE)
colnames(ndvi_df) <- c("x", "y", "ndvi_value")
ggplot(ndvi_df, aes(x = x, y = y, fill = ndvi_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,           
    option = "C",          
    guide = guide_colorbar(barwidth = 1, barheight = 15)  
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "NDVI Map"
  ) +
  theme(
    panel.grid = element_blank(),       
    plot.title = element_text(hjust = 0.5, size = 16),  
    axis.text = element_blank(),        
    axis.ticks = element_blank(),      
    axis.line = element_blank(),       
    axis.title = element_blank()      
  )
ggsave("/user/ziqitang/data_ziqi/NDVI_1124.png", width = 10, height = 8, dpi = 300)

#Moran I test under randomisation
#data:  ndvi_df$layer  
#weights: ndvi_listw    

#Moran I statistic standard deviate = 42.675, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#  Moran I statistic       Expectation          Variance 
#0.8604746167     -0.0008503401      0.0004073687 


raster_shannon_clean <- classify(dependent_clean, cbind(0, NA))
plot(raster_shannon_clean, main = "Map of shannon index")
#TAS independent variable
independent_raster2 <- rast("/user/ziqitang/data_ziqi/true_tas.tif")
raster_true_tas <- classify(independent_raster2, cbind(0, NA))
plot(raster_true_tas, main = "Map of true_tas")
true_tas_resampled <- resample(raster_true_tas,raster_shannon_clean, method = "bilinear")
plot(true_tas_resampled,main='true_tas_resampled')
combined_mask_tas <- mask(true_tas_resampled,raster_shannon_clean)
plot(combined_mask_tas, main = "combined_mask_tas")
combined_independent2 <- mask(true_tas_resampled,combined_mask_tas)
plot(combined_independent2, main = "TAS")

tas_df <- as.data.frame(combined_independent2, xy = TRUE, na.rm = TRUE)
colnames(tas_df) <- c("x", "y", "tas_value")
ggplot(tas_df, aes(x = x, y = y, fill = tas_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,         
    option = "C",          
    guide = guide_colorbar(barwidth = 1, barheight = 15)  
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "Surface Air Temperature (degree Celsius)"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove gridlines
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
    axis.text = element_blank(),        # Remove axis numeric labels
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.line = element_blank(),        # Remove axis lines
    axis.title = element_blank()        # Remove axis titles
  )
ggsave("/user/ziqitang/data_ziqi/TAS_1124.png", width = 10, height = 8, dpi = 300)



#dem/elevation independent variable
independent_raster3 <- rast("/user/ziqitang/data_ziqi/true_dem.tif") 
independent3_resampled <- resample(independent_raster3, dependent_clean, method = "bilinear")
combined_mask3 <- mask(independent3_resampled,dependent_clean)
combined_independent3 <- mask(independent3_resampled,combined_mask3)
plot(combined_independent3,main='DEM')
dem_df <- as.data.frame(combined_independent3, xy = TRUE, na.rm = TRUE)
colnames(dem_df) <- c("x", "y", "dem_value")  # Rename columns for clarity
ggplot(dem_df, aes(x = x, y = y, fill = dem_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,   # Label for the color bar
    option = "C",            # Viridis palette option
    guide = guide_colorbar(barwidth = 1, barheight = 15)  # Adjust color bar size
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "DEM (meters)"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove gridlines
    plot.title = element_text(hjust = 0.5, size = 16),  # Center title
    axis.text = element_blank(),        # Remove axis labels
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.line = element_blank(),        # Remove axis lines
    axis.title = element_blank()        # Remove axis titles
  )
ggsave("/user/ziqitang/data_ziqi/DEM_1124.png", width = 10, height = 8, dpi = 300)

elevation_df <- as.data.frame(combined_independent3, xy = TRUE, na.rm = TRUE)
names(elevation_df)[3] <- "layer" 
elevation_df$layer <- as.numeric(elevation_df$layer)
coordinates(elevation_df) <- ~x + y  
elevation_nb <- knn2nb(knearneigh(coordinates(elevation_df), k = 4))
elevation_listw <- nb2listw(elevation_nb, style = "W")
elevation_moran_result <- moran.test(elevation_df$layer, elevation_listw)
print(elevation_moran_result)
#Moran I test under randomisation
#data:  elevation_df$layer  
#weights: elevation_listw    
#Moran I statistic standard deviate = 43.837, p-value < 2.2e-16
#alternative hypothesis: greater
#sample estimates:
#  Moran I statistic       Expectation          Variance 
#0.8814582353     -0.0008460237      0.0004050920 

# calculate the tpi of dem
independent3_resampled <- resample(independent_raster3, dependent_clean, method = "bilinear")
combined_mask3 <- mask(independent3_resampled, dependent_clean)
combined_independent1 <- mask(independent3_resampled, combined_mask3)
tpi_raster <- terrain(combined_independent1, v = "TPI", neighbors = 4)
plot(tpi_raster, main = 'TPI from DEM')
# Divide each pixel by 100
combined_independent1 <- tpi_raster / 100
plot(combined_independent1, main = 'TPI from DEM')
tpi_df <- as.data.frame(combined_independent1, xy = TRUE, na.rm = TRUE)
colnames(tpi_df) <- c("x", "y", "tpi_value")  # Rename for clarity
ggplot(tpi_df, aes(x = x, y = y, fill = tpi_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,           # No label for the color bar
    option = "C",          # Viridis palette option (blue-green-yellow gradient)
    guide = guide_colorbar(barwidth = 1, barheight = 15)  # Adjust color bar size
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "TPI"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove gridlines
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
    axis.text = element_blank(),        # Remove axis numeric labels
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.line = element_blank(),        # Remove axis lines
    axis.title = element_blank()        # Remove axis titles
  )
ggsave("/user/ziqitang/data_ziqi/TPI_1124.png", width = 10, height = 8, dpi = 300)



# aspect independent variable66
independent_raster66 <- rast("/user/ziqitang/data_ziqi/true_aspect_area.tif") 
independent66_resampled <- resample(independent_raster66, dependent_clean, method = "bilinear")
combined_mask66 <- mask(independent66_resampled,dependent_clean)
combined_independent66 <- mask(independent66_resampled,combined_mask66)
plot(combined_independent66,main='Aspect')
aspect_df <- as.data.frame(combined_independent66, xy = TRUE, na.rm = TRUE)
names(aspect_df)[3] <- "layer" 
aspect_df$layer <- as.numeric(aspect_df$layer)
coordinates(aspect_df) <- ~x + y  
aspect_nb <- knn2nb(knearneigh(coordinates(aspect_df), k = 4))
aspect_listw <- nb2listw(aspect_nb, style = "W")
aspect_moran_result <- moran.test(aspect_df$layer, aspect_listw)
print(aspect_moran_result)

aspect_df <- as.data.frame(combined_independent66, xy = TRUE, na.rm = TRUE)
colnames(aspect_df) <- c("x", "y", "aspect_value")  # Rename columns for clarity

ggplot(aspect_df, aes(x = x, y = y, fill = aspect_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,           # No label for the color bar
    option = "C",          # Viridis palette option (blue-green-yellow gradient)
    guide = guide_colorbar(barwidth = 1, barheight = 15)  # Adjust color bar size
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "Aspect (radians)"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove gridlines
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
    axis.text = element_blank(),        # Remove axis numeric labels
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.line = element_blank(),        # Remove axis lines
    axis.title = element_blank()        # Remove axis titles
  )
ggsave("/user/ziqitang/data_ziqi/ASPECT_1124.png", width = 10, height = 8, dpi = 300)


# fire age  independent variable5
true_temp_fireAge <- rast("/user/ziqitang/data_ziqi/years_since_fire.tif")
raster_true_fireAge <- classify(true_temp_fireAge, cbind(0, NA))
crs_fireAge <- crs(raster_true_fireAge)
if (crs_dependent != crs_fireAge) {
  raster_true_fireAge_projected <- project(raster_true_fireAge, crs(dependent_clean))
} else {
  raster_true_fireAge_projected <- raster_true_fireAge
}
plot(raster_true_fireAge_projected, main = "Fire Age Raster Projected")
true_fireAge_resampled <- resample(raster_true_fireAge_projected,dependent_clean, method = "bilinear")
combined_mask_fireAge <- mask(true_fireAge_resampled,dependent_clean)
combined_true_fireAge <- mask(true_fireAge_resampled,combined_mask_fireAge)
plot(combined_true_fireAge, main = "Fire Age")
output_tif <- "/user/ziqitang/data_ziqi/combined_true_fireAge.tif"
writeRaster(combined_true_fireAge, filename = output_tif,overwrite = TRUE)
combined_independent4 <- rast("/user/ziqitang/data_ziqi/combined_true_fireAge.tif") 
plot(combined_independent4,main='combined_independent4')

fireage_df <- as.data.frame(combined_independent4, xy = TRUE, na.rm = TRUE)
colnames(fireage_df) <- c("x", "y", "fireage_value")  # Rename columns for clarity
ggplot(fireage_df, aes(x = x, y = y, fill = fireage_value)) +
  geom_tile() +
  scale_fill_viridis_c(
    name = NULL,   # Label for the color bar
    option = "C",            # Viridis palette option
    guide = guide_colorbar(barwidth = 1, barheight = 15)  # Adjust color bar size
  ) +
  theme_minimal() +
  coord_fixed() +
  labs(
    title = "Fire Age (years)"
  ) +
  theme(
    panel.grid = element_blank(),       # Remove gridlines
    plot.title = element_text(hjust = 0.5, size = 16),  # Center title
    axis.text = element_blank(),        # Remove axis labels
    axis.ticks = element_blank(),       # Remove axis ticks
    axis.line = element_blank(),        # Remove axis lines
    axis.title = element_blank()        # Remove axis titles
  )
ggsave("/user/ziqitang/data_ziqi/fireage_1124.png", width = 10, height = 8, dpi = 300)



plot(dependent_raster)
dependent_df <- as.data.frame(dependent_raster, xy = TRUE, na.rm = TRUE)
dependent_cell_count <- ncell(dependent_raster)
print(dependent_cell_count)
#Ndvi
plot(combined_independent0)
independent_df0 <- as.data.frame(combined_independent0, xy = TRUE, na.rm = TRUE)
independent_df0_cell_count <- ncell(combined_independent0)
print(independent_df0_cell_count)
#tpi
independent_df1 <- as.data.frame(combined_independent1, xy = TRUE, na.rm = TRUE)
independent_df1_cell_count <- ncell(independent_df1)
print(independent_df1_cell_count)
#tas
independent_df2 <- as.data.frame(combined_independent2, xy = TRUE, na.rm = TRUE)
independent_df2_cell_count <- ncell(independent_df2)
print(independent_df2_cell_count)
#dem/elevation
independent_df3 <- as.data.frame(combined_independent3, xy = TRUE, na.rm = TRUE)
independent_df3_cell_count <- ncell(independent_df3)
print(independent_df3_cell_count)
#fire age
independent_df4 <- as.data.frame(combined_independent4, xy = TRUE, na.rm = TRUE)
independent_df4_cell_count <- ncell(independent_df4)
print(independent_df4_cell_count)
#aspect
#independent_df6 <- as.data.frame(combined_independent6, xy = TRUE, na.rm = TRUE)
independent_df66 <- as.data.frame(combined_independent66, xy = TRUE, na.rm = TRUE)
independent_df66_cell_count <- ncell(independent_df66)
print(independent_df66_cell_count)





# Define a function to summarize raster data
summarize_raster <- function(raster, name) {
  values <- values(raster, na.rm = TRUE)
  data.frame(
    Name = name,
    Min = min(values, na.rm = TRUE),
    Max = max(values, na.rm = TRUE),
    Mean = mean(values, na.rm = TRUE),
    Median = median(values, na.rm = TRUE),
    StdDev = sd(values, na.rm = TRUE),
    Count = length(values)
  )
}


library(GGally)

# 设置主题和样式
custom_theme <- function() {
  theme_minimal() +
    theme(
      strip.background = element_blank(),
      axis.text = element_text(size = 8),
      strip.text = element_text(size = 10)
    )
}

ggpairs(combined_data[,c("Alpha_Diversity", "NDVI", "TPI", "TAS", "Elevation", "Fire_Age", "Aspect")],
        lower = list(continuous = "points"),  # 散点图
        diag = list(continuous = "density"),  # 密度图
        upper = list(continuous = "cor"),     # 相关系数
        title = "Scatterplot Matrix of Environmental Variables and Shannon Index") +
  theme_minimal() +
  theme(
    axis.text = element_text(size = 8),
    strip.text = element_text(size = 10)
  )

ggsave("/user/ziqitang/data_ziqi/correlation_matrix.png", width = 12, height = 10, dpi = 300)




# Perform descriptive analysis for each raster
ndvi_summary <- summarize_raster(combined_independent0, "NDVI")
tpi_summary <- summarize_raster(combined_independent1, "TPI")
tas_summary <- summarize_raster(combined_independent2, "Temperature (TAS)")
dem_summary <- summarize_raster(combined_independent3, "DEM/Elevation")
fire_age_summary <- summarize_raster(combined_independent4, "Fire Age")
aspect_summary <- summarize_raster(combined_independent66, "Aspect")

# Combine all summaries into a single data frame for easy comparison
summary_df <- bind_rows(ndvi_summary, tpi_summary, tas_summary, dem_summary, fire_age_summary, aspect_summary)
# Display the summary
print(summary_df)



#multicollinary analysis
combined_data1 <- merge(dependent_df, independent_df0, by = c("x", "y"))
combined_data1 <- merge(combined_data1, independent_df1, by = c("x", "y"))
combined_data1 <- merge(combined_data1, independent_df2, by = c("x", "y"))
combined_data1 <- merge(combined_data1, independent_df3, by = c("x", "y"))
#combined_data1 <- merge(combined_data1, independent_df4, by = c("x", "y"))
combined_data1 <- merge(combined_data1, independent_df66, by = c("x", "y"))
#combined_data <- merge(combined_data, independent_df6, by = c("x", "y"))
colnames(combined_data1) <- c("x", "y", "shannon_index", "var1", "var2", "var3","var4","var5")
combined_data1 <- combined_data1[, -c(1, 2)]  # Exclude x, y columns if they are coordinates
str(combined_data1)
colnames(combined_data1) <- c("Alpha Diversity", "NDVI", "TPI", "TAS", "Elevation", "Fire Age", "Aspect")
combined_data1_clean <- na.omit(combined_data1)
ggpairs(combined_data1_clean, 
        title = "Scatterplot Matrix of Environmental Variables and Shannon Index",
        mapping = aes(alpha = 0.5)) +
  theme_minimal()






#0236
combined_data <- merge(dependent_df, independent_df0, by = c("x", "y"))
combined_data <- merge(combined_data, independent_df2, by = c("x", "y"))
combined_data <- merge(combined_data, independent_df3, by = c("x", "y"))
combined_data <- merge(combined_data, independent_df66, by = c("x", "y"))
colnames(combined_data) <- c("x", "y", "shannon_index", "var1", "var2", "var3","var4")
# 构建 GWR 模型
combined_data$x_coord <- combined_data$x
combined_data$y_coord <- combined_data$y
coordinates(combined_data) <- ~ x + y
gwr_bandwidth <- gwr.sel(shannon_index ~ var1 + var2 + var3 + var4 , data = combined_data, adapt = TRUE)
print(gwr_bandwidth)
gwr_model <- gwr(shannon_index ~ var1 + var2 + var3 + var4, data = combined_data, bandwidth = gwr_bandwidth)
summary(gwr_model)






# validation rf交叉验证
# 1. 数据准备和分割
set.seed(123)
data_split <- initial_split(combined_data, prop = 0.8)
train_data <- training(data_split)
test_data <- testing(data_split)
# 2. 创建10折交叉验证
cv_folds <- vfold_cv(train_data, v = 10)
# 3. RF模型设置
rf_spec <- rand_forest(trees = 500) %>% 
  set_engine("ranger") %>% 
  set_mode("regression")
rf_workflow <- workflow() %>%
  add_formula(shannon_index ~ var1 + var2 + var3 + var4) %>%
  add_model(rf_spec)
# 4. GRF模型设置
grf_spec <- rand_forest(trees = 500) %>% 
  set_engine("ranger") %>%
  set_mode("regression")
grf_workflow <- workflow() %>%
  add_formula(shannon_index ~ var1 + var2 + var3 + var4 + x + y) %>%
  add_model(grf_spec)
# 5. 定义评估指标
metrics <- metric_set(rmse, rsq, mae)
# 6. 进行交叉验证
rf_results <- fit_resamples(
  rf_workflow,
  resamples = cv_folds,
  metrics = metrics
)
grf_results <- fit_resamples(
  grf_workflow,
  resamples = cv_folds,
  metrics = metrics
)
# 7. 收集并比较结果
model_comparison <- bind_rows(
  collect_metrics(rf_results) %>% mutate(model = "RF"),
  collect_metrics(grf_results) %>% mutate(model = "GRF")
)
# 打印比较结果
print(model_comparison)



library(ggplot2)
predicted_values <- gwr_model$SDF$pred
combined_data_df <- as.data.frame(combined_data)
combined_data <- as.data.frame(combined_data)
combined_data$predicted_shannon <- gwr_model$SDF$pred

if (!("shannon_index" %in% colnames(combined_data))) {
  combined_data$shannon_index <- combined_data_df$shannon_index  # Replace with the correct source
}

if (!("predicted_shannon" %in% colnames(combined_data)) || !("shannon_index" %in% colnames(combined_data))) {
  stop("Error: combined_data does not contain the required columns 'predicted_shannon' and 'shannon_index'")
}


# predicted value & actual value
lm_fit <- lm(shannon_index ~ predicted_shannon, data = combined_data)
r2 <- summary(lm_fit)$r.squared
equation <- sprintf("y = %.4fx + %.3f\nR² = %.3f", 
                    coef(lm_fit)[2], 
                    coef(lm_fit)[1],
                    r2)
# 生成图形
my_plot <- ggplot(combined_data, aes(x = predicted_shannon, y = shannon_index)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # 添加回归线
  annotate("text", 
           x = min(combined_data$predicted_shannon), 
           y = max(combined_data$shannon_index),
           label = equation,
           hjust = 0,
           vjust = 1) +
  labs(x = "Predicted Values (GWR)",
       y = "Actual Values (Shannon Index)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
plot(my_plot)
# 保存图形为 PNG 文件
ggsave("/user/ziqitang/data_ziqi/my_plot_gwr.png", plot = my_plot, width = 8, height = 6, dpi = 300)


# 构建 OLS 模型
ols_model <- lm(shannon_index ~ var1 + var2 + var3 + var4, data = combined_data)
summary(ols_model)


combined_data$predicted_shannon_ols <- predict(ols_model, combined_data)

lm_fit_ols <- lm(shannon_index ~ predicted_shannon_ols, data = combined_data)
r2_ols <- summary(lm_fit_ols)$r.squared
equation_ols <- sprintf("y = %.4fx + %.3f\nR² = %.3f", 
                        coef(lm_fit_ols)[2], 
                        coef(lm_fit_ols)[1],
                        r2_ols)
# 生成 OLS 图形
my_ols_plot <- ggplot(combined_data, aes(x = predicted_shannon_ols, y = shannon_index)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # 添加回归线
  annotate("text", 
           x = min(combined_data$predicted_shannon_ols), 
           y = max(combined_data$shannon_index),
           label = equation_ols,
           hjust = 0,
           vjust = 1) +
  labs(x = "Predicted Values (OLS)",
       y = "Actual Values (Shannon Index)") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
plot(my_ols_plot)
# 保存 OLS 图形为 PNG 文件
ggsave("/user/ziqitang/data_ziqi/my_plot_ols.png", plot = my_ols_plot, width = 8, height = 6, dpi = 300)
# Determine the global x and y range for both plots
x_range <- range(c(combined_data$predicted_shannon, combined_data$predicted_shannon_ols), na.rm = TRUE)
y_range <- range(combined_data$shannon_index, na.rm = TRUE)

# Generate GWR plot
my_plot <- ggplot(combined_data, aes(x = predicted_shannon, y = shannon_index)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add regression line
  annotate("text", 
           x = x_range[1], 
           y = y_range[2],
           label = equation,
           hjust = 0,
           vjust = 1) +
  labs(x = "Predicted Values (GWR)",
       y = "Actual Values (Shannon Index)") +
  coord_cartesian(xlim = x_range, ylim = y_range) +  # Ensure consistent axis range
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
plot(my_plot)
# Save GWR plot
ggsave("/user/ziqitang/data_ziqi/my_plot_gwr.png", plot = my_plot, width = 8, height = 6, dpi = 300)

# Generate OLS plot
my_ols_plot <- ggplot(combined_data, aes(x = predicted_shannon_ols, y = shannon_index)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # Add regression line
  annotate("text", 
           x = x_range[1], 
           y = y_range[2],
           label = equation_ols,
           hjust = 0,
           vjust = 1) +
  labs(x = "Predicted Values (OLS)",
       y = "Actual Values (Shannon Index)") +
  coord_cartesian(xlim = x_range, ylim = y_range) +  # Ensure consistent axis range
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid = element_blank(),
        axis.line = element_line(color = "black"))
plot(my_ols_plot)
# Save OLS plot
ggsave("/user/ziqitang/data_ziqi/my_plot_ols.png", plot = my_ols_plot, width = 8, height = 6, dpi = 300)
# 加载必要的库
library(ggplot2)
library(randomForest)

# 1. 构建 GRF 模型
set.seed(123)  # 设置随机种子，保证结果可复现

# 数据准备：增加坐标
combined_data$x_coord <- combined_data$x
combined_data$y_coord <- combined_data$y

# 构建 GRF 模型
grf_model <- randomForest(shannon_index ~ var1 + var2 + var3 + var4 + x_coord + y_coord, 
                          data = combined_data, 
                          ntree = 500,    # 决策树数量
                          importance = TRUE)  # 追踪变量重要性

# 2. 生成预测值
combined_data$predicted_shannon_grf <- predict(grf_model, newdata = combined_data)

# 3. 计算回归线方程和R²值
lm_fit_grf <- lm(shannon_index ~ predicted_shannon_grf, data = combined_data)
r2_grf <- summary(lm_fit_grf)$r.squared
equation_grf <- sprintf("y = %.4fx + %.3f\nR² = %.3f", 
                        coef(lm_fit_grf)[2], 
                        coef(lm_fit_grf)[1],
                        r2_grf)

# 4. 生成 GRF 散点图
my_grf_plot <- ggplot(combined_data, aes(x = predicted_shannon_grf, y = shannon_index)) +
  geom_point(color = "black") +  # 绘制散点
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +  # 添加45度线
  geom_smooth(method = "lm", color = "black", se = FALSE) +  # 添加回归线
  annotate("text", 
           x = min(combined_data$predicted_shannon_grf), 
           y = max(combined_data$shannon_index),
           label = equation_grf,
           hjust = 0, 
           vjust = 1) +  # 添加方程和R²值
  labs(x = "Predicted Values (GRF)",
       y = "Actual Values (Shannon Index)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),       # 移除网格线
    axis.line = element_line(color = "black")  # 添加坐标轴线
  )

# 5. 显示图形
plot(my_grf_plot)
ggsave("/user/ziqitang/data_ziqi/my_plot_grf.png", plot = my_grf_plot, width = 8, height = 6, dpi = 300)
# 加载必要的库








# 获取GWR模型的预测值
predictions_gwr <- gwr_model$SDF$pred  # 假设gwr_model存储了预测值在SDF$pred中
actual_values <- combined_data$shannon_index
tss <- sum((actual_values - mean(actual_values))^2)
rss <- sum((actual_values - predictions_gwr)^2)
r_squared_gwr <- 1 - (rss / tss)
print(paste("GWR R-squared:", r_squared_gwr))


# Extract R-squared values from the GWR model
combined_data$r_squared <- gwr_model$SDF$localR2

# Ensure combined_data has coordinates in standard columns (x and y)
combined_data <- as.data.frame(combined_data)  # Convert back to a data frame if needed

# Plot using ggplot2
ggplot(data = combined_data, aes(x = x, y = y, color = r_squared)) +
  geom_point(size = 3) +  # Plot points with size adjustment
  scale_color_gradientn(colors = c("blue", "green", "yellow", "red"),
                        name = "Local R-squared") +  # Custom gradient
  labs(title = "GWR Local R-squared Map",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +  # Clean theme
  theme(plot.title = element_text(hjust = 0.5, size = 16))  # Center and style title

ggplot(data = combined_data, aes(x = x, y = y, fill = r_squared)) +
  geom_tile() +  # Use tiles to create a heatmap-like plot
  scale_fill_gradientn(colors = c("purple", "blue", "yellow", "red"),
                       name = "Local R²") +  # Custom gradient
  labs(title = "Local R² Values of the GWR Model",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +  # Clean theme
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
    axis.title = element_text(size = 12),              # Adjust axis label size
    legend.title = element_text(size = 12),            # Adjust legend title size
    legend.text = element_text(size = 10)              # Adjust legend text size
  )

ggplot(data = combined_data, aes(x = x, y = y, fill = r_squared)) +
  geom_tile() +  # Use tiles for a heatmap-like plot
  scale_fill_gradientn(
    colors = c("purple", "blue", "green", "yellow"),  # Adjust to match your example
    limits = c(0, 1),  # Ensures the scale is from 0 to 1
    name = "Local R²"
  ) +
  labs(title = "Local R² Values of the GWR Model",
       x = "Longitude",
       y = "Latitude") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
#修正的r squared value of the map
# Correct assignment of the ggplot object
R_squared_gwr <- ggplot(data = combined_data, aes(x = x, y = y, fill = r_squared)) +
  geom_tile() +  # Use tiles for a heatmap-like plot
  scale_fill_gradientn(
    colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"),  # Custom palette
    limits = c(0, 1),  # Scale from 0 to 1
    name = "Local R²"
  ) +
  labs(title = "Local R² Values of the GWR Model",
       x = "Longitude",
       y = "Latitude") +
  theme_void() +  # 移除背景和多余的轴
  guides(
    fill = guide_colorbar(
      barwidth = 1,       # 设置颜色条宽度
      barheight = 15      # 设置颜色条高度
    )
  ) +
  theme(
    legend.position = "right",      # 确保颜色条在右侧
    legend.title = element_blank(), # 删除标题
    legend.text = element_text(size = 10)  # 设置刻度字体大小
  ) 
plot(R_squared_gwr)
# Save the plot to file
ggsave("/user/ziqitang/data_ziqi/R_squared_value.png", plot = R_squared_gwr, width = 8, height = 15, dpi = 300)





# 标准化自变量和因变量
combined_data$shannon_index <- scale(combined_data$shannon_index)
combined_data$var1 <- scale(combined_data$var1)
combined_data$var2 <- scale(combined_data$var2)
combined_data$var3 <- scale(combined_data$var3)
combined_data$var4 <- scale(combined_data$var4)
library(spgwr)

# Define spatial coordinates and projection
epsg_code <- 4326
coordinates(combined_data) <- ~x + y  
proj4string(combined_data) <- CRS(SRS_string = paste0("EPSG:", epsg_code))  # Define CRS
head(coordinates(combined_data))
proj4string(combined_data)
variables <- c("var1", "var2", "var3", "var4")
variable_names <- c("var1" = "NDVI", "var2" = "TAS", "var3" = "Elevation", "var4" = "Aspect")
for (var in variables) {
  
  # Step 1: Generate GWR Model
  formula <- as.formula(paste("shannon_index ~", var))
  gwr_bandwidth <- gwr.sel(formula, data = combined_data, adapt = TRUE)  # Select adaptive bandwidth
  
  gwr_model <- gwr(
    formula = formula,
    data = combined_data,
    adapt = gwr_bandwidth  # Use adaptive bandwidth
  )
  combined_data[[paste0("coef_", var)]] <- gwr_model$SDF[[var]]
  coef_raster <- rasterFromXYZ(as.data.frame(combined_data)[, c("x", "y", paste0("coef_", var))])
  crs(coef_raster) <- proj4string(combined_data)  # Set CRS from spatial data
  coef_df <- as.data.frame(coef_raster, xy = TRUE, na.rm = TRUE)
  colnames(coef_df) <- c("x", "y", "coef_value")  # Rename columns
  coef_df$coef_category <- cut(
    coef_df$coef_value,
    breaks = c(-Inf, -1, -0.5, 0.5, 1, Inf),  # Define intervals
    labels = c("< -1", "-1 to -0.5", "-0.5 to 0.5", "0.5 to 1", "> 1"),
    include.lowest = TRUE
  )
  plot_title <- paste("GWR Coefficient for", variable_names[var])  # Dynamic title based on variable
  
  gg <- ggplot(coef_df, aes(x = x, y = y, fill = coef_category)) +
    geom_tile() +
    scale_fill_manual(
      name = "Coefficient",   # Legend title
      values = c(
        "< -1" = "#d73027",    # Red for strong negatives
        "-1 to -0.5" = "#fdae61",  # Orange for weak negatives
        "-0.5 to 0.5" = "#ffffbf", # White for near zero
        "0.5 to 1" = "#abd9e9",    # Light blue for weak positives
        "> 1" = "#4575b4"       # Blue for strong positives
      )
    ) +
    coord_fixed() +
    labs(title = plot_title) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),       # Remove gridlines
      plot.title = element_text(hjust = 0.5, size = 16),  # Center the title
      axis.text = element_blank(),        # Remove axis numeric labels
      axis.ticks = element_blank(),       # Remove axis ticks
      axis.line = element_blank(),        # Remove axis lines
      axis.title = element_blank()        # Remove axis titles
    )
  
  # Save the Plot
  output_file <- paste0("/user/ziqitang/data_ziqi/", var, "_GWR_Diverging_Map.png")
  ggsave(output_file, plot = gg, width = 10, height = 8, dpi = 300)
  
  # Print message for progress
  print(paste("Completed GWR and plot for", variable_names[var], "- saved to", output_file))
}


#feature importance
varImpPlot(rf_model, main = "Feature Importance for Shannon Index Prediction")
print(rf_model$ntree)
feature_importance_df <- as.data.frame(importance(rf_model))
colnames(feature_importance_df) <- "FeatureImportance"
feature_importance_df$Feature <- rownames(feature_importance_df)
feature_importance_df <- feature_importance_df[order(feature_importance_df$FeatureImportance, decreasing = TRUE),]
library(dplyr)
feature_importance_df$Feature <- dplyr::recode(feature_importance_df$Feature,
                                               "var1" = "NDVI",
                                               "var2" = "Mean Daily Air Temperature(TAS)",
                                               "var3" = "Elevation",
                                               "var4" = "Aspect")










# Fit the OLS model
ols_model <- lm(shannon_index ~ var1 + var2 + var3 + var4, data = combined_data)
r_squared_ols <- summary(ols_model)$r.squared
adjusted_r_squared_ols <- summary(ols_model)$adj.r.squared
rmse_ols <- sqrt(mean(residuals(ols_model)^2))
aic_ols <- AIC(ols_model)
print(paste("OLS R-squared:", r_squared_ols))
print(paste("OLS Adjusted R-squared:", adjusted_r_squared_ols))
print(paste("OLS RMSE:", rmse_ols))
print(paste("OLS AIC:", aic_ols))

#regression coefficients of the ols model
# Get summary statistics
ols_summary <- summary(ols_model)
r_squared_ols <- ols_summary$r.squared
adjusted_r_squared_ols <- ols_summary$adj.r.squared
rmse_ols <- sqrt(mean(residuals(ols_model)^2))
aic_ols <- AIC(ols_model)

print(paste("OLS R-squared:", round(r_squared_ols, 3)))
print(paste("OLS Adjusted R-squared:", round(adjusted_r_squared_ols, 3)))
print(paste("OLS RMSE:", round(rmse_ols, 3)))
print(paste("OLS AIC:", round(aic_ols, 3)))
combined_data$predicted_shannon_ols <- predict(ols_model, newdata = combined_data)
library(ggplot2)
#title = "Scatter Plot of Actual vs Predicted OLS Values",
ggplot(combined_data, aes(x = predicted_shannon_ols, y = shannon_index)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
  labs(  x = "Predicted Values (OLS)",
         y = "Actual Values (Shannon Index)") +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),  
    panel.grid = element_blank(),           
    axis.line = element_line(color = "black") 
  )








# Extract coefficients and p-values
coefficients <- ols_summary$coefficients
coef_table <- data.frame(
  Term = rownames(coefficients),
  Estimate = coefficients[, "Estimate"],
  Std.Error = coefficients[, "Std. Error"],
  t.value = coefficients[, "t value"],
  p.value = coefficients[, "Pr(>|t|)"]
)




summary(combined_data@data)
any(is.na(combined_data@data))  


combined_data <- combined_data[complete.cases(combined_data@data), ]



# Fit the GWR model (assuming you've calculated gwr_bandwidth as shown earlier)
library(spgwr)
gwr_model <- gwr(shannon_index ~ var1 + var2 + var3 + var4, data = combined_data, bandwidth = gwr_bandwidth)
r_squared_gwr <- gwr_model$SDF$localR2
mean_r_squared_gwr <- mean(r_squared_gwr)  # Use the mean R-squared across locations
predictions_gwr <- gwr_model$SDF$pred  # Assuming gwr_model stores predictions in SDF$pred
rmse_gwr <- sqrt(mean((combined_data$shannon_index - predictions_gwr)^2))
aic_gwr <- gwr_model$results$AICc
n <- nrow(combined_data)
p <- 4  
adjusted_r_squared_gwr <- 1 - ((1 - mean_r_squared_gwr) * (n - 1) / (n - p - 1))
print(paste("GWR Mean R-squared:", mean_r_squared_gwr))
print(paste("GWR Adjusted R-squared:", adjusted_r_squared_gwr))
print(paste("GWR RMSE:", rmse_gwr))
print(paste("GWR AIC:", aic_gwr))

# Calculate RSS
rss <- sum((combined_data$shannon_index - predictions_gwr)^2)

# Manually calculate AICc
aicc_gwr <- n * log(rss / n) + 2 * p + (2 * p * (p + 1)) / (n - p - 1)
print(paste("Manually calculated GWR AICc:", aicc_gwr))



#GRF
library(randomForest)
combined_data$x_coord <- combined_data$x
combined_data$y_coord <- combined_data$y
combined_data_model <- combined_data[, c("shannon_index", "var1", "var2", "var3", "var4", "x_coord", "y_coord")]
set.seed(123) 
grf_model <- randomForest(shannon_index ~ var1 + var2 + var3 + var4 + x_coord + y_coord, 
                          data = combined_data_model, 
                          ntree = 500,  # Number of trees
                          importance = TRUE)  # Track variable importance
predictions_grf <- predict(grf_model, newdata = combined_data)
actual_values <- combined_data$shannon_index
tss_grf <- sum((actual_values - mean(actual_values))^2)
rss_grf <- sum((actual_values - predictions_grf)^2)
r_squared_grf <- 1 - (rss_grf / tss_grf)
print(paste("GRF R-squared:", r_squared_grf))

library(ggplot2)


combined_data_df <- data.frame(
  actual_shannon = actual_values,         
  predicted_shannon_grf = predictions_grf
)


ggplot(combined_data_df, aes(x = predicted_shannon_grf, y = actual_shannon)) +
  geom_point(color = "black") + 
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
  labs(
    x = "Predicted Values (GRF)",          
    y = "Actual Values (Shannon Index)",   
    title = "Scatter Plot: GRF Predictions vs Actual Values"  
  ) +
  theme_minimal() +  
  theme(
    plot.title = element_text(hjust = 0.5, size = 14),  
    panel.grid = element_blank(),                    
    axis.line = element_line(color = "black"),       
    axis.text = element_text(size = 10),             
    axis.title = element_text(size = 12)              
  )

x_range <- range(predictions_grf, na.rm = TRUE)
y_range <- range(actual_values, na.rm = TRUE)
print(predictions_grf)
print(actual_values)

lm_fit_grf <- lm(actual_values ~ predictions_grf)
r2_grf <- summary(lm_fit_grf)$r.squared
equation_grf <- sprintf("y = %.4fx + %.3f\nR² = %.3f", 
                        coef(lm_fit_grf)[2], 
                        coef(lm_fit_grf)[1],
                        r2_grf)


my_grf_plot <- ggplot(combined_data, aes(x = predictions_grf, y = actual_values)) +
  geom_point(color = "black") +
  geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") + 
  geom_smooth(method = "lm", color = "black", se = FALSE) +  
  annotate("text", 
           x = x_range[1], 
           y = y_range[2],
           label = equation_grf,
           hjust = 0,
           vjust = 1) +  
  labs(x = "Predicted Values (GRF)",
       y = "Actual Values (Shannon Index)") +
  coord_cartesian(xlim = x_range, ylim = y_range) +  
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),       
    axis.line = element_line(color = "black") 
  )

set.seed(123)
folds <- createFolds(combined_data$shannon_index, k = 5)
fit_gwr_cv <- function(train_data, test_data) {
  bandwidth <- gwr.sel(shannon_index ~ var1 + var2 + var3 + var4, data = train_data, adapt = TRUE)
  gwr_model <- gwr(shannon_index ~ var1 + var2 + var3 + var4, data = train_data, bandwidth = bandwidth)
  local_coefficients <- gwr_model$SDF
  predicted_test <- local_coefficients$`(Intercept)` + 
    local_coefficients$var1 * test_data$var1 + 
    local_coefficients$var2 * test_data$var2 +
    local_coefficients$var3 * test_data$var3 +
    local_coefficients$var4 * test_data$var4
  actual_test <- test_data$shannon_index
  rmse_test <- sqrt(mean((actual_test - predicted_test)^2))
  return(rmse_test)
}
rmse_values <- sapply(folds, function(train_indices) {
  train_data <- combined_data[train_indices, ]
  test_data <- combined_data[-train_indices, ]
  fit_gwr_cv(train_data, test_data)
})
print(rmse_values)
mean_rmse <- mean(rmse_values)
print(paste("Average RMSE across folds:", mean_rmse))



# calculate the rmse of ols model
ols_model <- lm(shannon_index ~var1 + var2 + var3 + var4, data = combined_data)
summary(ols_model)
predicted_ols <- predict(ols_model, newdata = combined_data)
residuals_squared_mean <- sqrt(mean((combined_data$shannon_index - predicted_ols)^2))
rmse_ols <- sqrt(residuals_squared_mean)
print(rmse_ols)
actual_values <- combined_data$shannon_index
predicted_ols <- predict(ols_model, newdata = combined_data)

rss <- sum((actual_values - predicted_ols)^2)  
tss <- sum((actual_values - mean(actual_values))^2)  

r_squared <- 1 - (rss / tss)
print(paste("Global R² for OLS Model:", r_squared))
plot_data_ols <- data.frame(
  Actual = combined_data$shannon_index,
  Predicted = predicted_ols
)
# Scatter Plot of Actual vs Predicted for OLS
ggplot(plot_data_ols, aes(x = Predicted, y = Actual)) +
  geom_point(alpha = 0.6) +
  geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "red") +  # y = x reference line
  labs(title = "Scatter Plot of Actual vs Predicted OLS Values",
       x = "Predicted Values (OLS)",
       y = "Actual Values") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))




coordinates(combined_data) <- ~x + y  # Specify spatial coordinates
proj4string(combined_data) <- CRS("+proj=longlat +datum=WGS84")  # Define the CRS (adjust as needed)

# Residual Plot: Residuals vs Predicted Values
ggplot(data = plot_data_ols, aes(x = Predicted, y = ols_residuals)) +
  geom_point(alpha = 0.6) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "blue") +  # Reference line at y = 0
  labs(title = "Residuals vs Predicted OLS Values",
       x = "Predicted Values (OLS)",
       y = "Residuals") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))

plot(combined_data$shannon_index, predicted_ols, xlab = "Actual Values", ylab = "Predicted Values", main = "Actual vs Predicted (OLS)")
abline(0, 1, col = "red")
print(paste("GWR RMSE:", rmse_gwr))
print(paste("OLS RMSE:", rmse_ols))
aic_ols <- AIC(ols_model)
print(aic_ols)
print(paste("OLS AIC:", aic_ols))
predicted_ols <- predict(ols_model, newdata = combined_data)
combined_data$predicted_ols <- predicted_ols
# Ensure both predicted values are available in combined_data
combined_data$predicted_gwr <- predicted_values
combined_data$predicted_ols <- predicted_ols

# Define color palette and range for consistency
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
min_val <- min(c(combined_data$predicted_ols, combined_data$predicted_gwr), na.rm = TRUE)
max_val <- max(c(combined_data$predicted_ols, combined_data$predicted_gwr), na.rm = TRUE)

# Plot OLS Predicted Values Map
spplot(combined_data, "predicted_ols", main = "Predicted OLS Values",
       col.regions = color_palette(100),  # Color palette for OLS
       at = seq(min_val, max_val, length.out = 100))  # Consistent range

# Plot GWR Predicted Values Map
spplot(combined_data, "predicted_gwr", main = "Predicted GWR Values",
       col.regions = color_palette(100),  # Color palette for GWR
       at = seq(min_val, max_val, length.out = 100))  # Consistent range

# Scatterplot Matrix for Independent Variables
pairs(combined_data[, c("var1", "var2", "var3", "var4")],
      main = "Scatterplot Matrix of Variables",
      pch = 19,       # Point shape
      col = "blue")   # Point color

# Create rasters for OLS and GWR predictions
ols_raster <- rasterFromXYZ(cbind(combined_data$x_coord, combined_data$y_coord, combined_data$predicted_ols))
gwr_raster <- rasterFromXYZ(cbind(combined_data$x_coord, combined_data$y_coord, combined_data$predicted_gwr))

# Define color palette and range
color_palette <- colorRampPalette(c("blue", "yellow", "red"))
min_val <- min(c(combined_data$predicted_ols, combined_data$predicted_gwr), na.rm = TRUE)
max_val <- max(c(combined_data$predicted_ols, combined_data$predicted_gwr), na.rm = TRUE)

# Set up side-by-side plotting
par(mfrow = c(1, 2))  # 1 row, 2 columns

# Plot OLS Predicted Values
plot(ols_raster, main = "Predicted OLS Values", 
     col = color_palette(100), zlim = c(min_val, max_val),
     legend.args = list(text = 'Predicted Values', side = 4, line = 2.5))

# Plot GWR Predicted Values
plot(gwr_raster, main = "Predicted GWR Values", 
     col = color_palette(100), zlim = c(min_val, max_val),
     legend.args = list(text = 'Predicted Values', side = 4, line = 2.5))

# Reset layout to default after plotting
par(mfrow = c(1, 1))

color_palette_leaflet <- colorNumeric(palette = c("blue", "yellow", "red"), 
                                      domain = c(min_val, max_val), na.color = "transparent")





#scatter plor diagram
# Step 1: Calculate VIF using a linear model
vif_model <- lm(shannon_index ~ var1 + var2 + var3 + var4, data = combined_data)
vif_values <- vif(vif_model)  # Calculate VIF for each predictor
print("VIF values for each variable:")
print(vif_values)

# 三个图片组合在一起
# OLS 地图
min_value <- min(c(
  min(combined_data_df$predicted_shannon_ols, na.rm = TRUE),
  min(combined_data_df$predicted_shannon_rf, na.rm = TRUE)
))
max_value <- max(c(
  max(combined_data_df$predicted_shannon_ols, na.rm = TRUE),
  max(combined_data_df$predicted_shannon_rf, na.rm = TRUE)
))

colnames(combined_data) <- make.names(colnames(combined_data), unique = TRUE)
combined_data_df <- as.data.frame(combined_data)
print(colnames(combined_data_df))
set.seed(123)  
rf_model <- randomForest(
  shannon_index ~ var1 + var2 + var3 + var4,
  data = combined_data_df,
  ntree = 500,
  importance = TRUE  
)
combined_data$predicted_shannon_rf <- predict(rf_model, newdata = combined_data_df)
plot_rf <- ggplot(combined_data_df, aes(x = x_coord, y = y_coord, fill = predicted_shannon_rf)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"), 
                       limits = c(min_value, max_value) ) +
  coord_fixed() +
  labs(title = "RF") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none" # 隐藏颜色条
  )
plot(plot_rf)

plot_ols <- ggplot(combined_data_df, aes(x = x_coord, y = y_coord, fill = predicted_shannon_ols)) +
  geom_tile() +
  scale_fill_gradientn( colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"),  # 自定义颜色
                        limits = c(min_value, max_value) )+
  coord_fixed() +
  labs(title = "OLS") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none" 
  )

combined_data_df <- as.data.frame(combined_data)
combined_data_df$x_coord <- combined_data$x
combined_data_df$y_coord <- combined_data$y
combined_data_df$predicted_shannon_gwr <- gwr_model$SDF$pred





min_value <- min(combined_data_df$predicted_shannon_gwr, na.rm = TRUE)
max_value <- max(combined_data_df$predicted_shannon_gwr, na.rm = TRUE)

plot_gwr <- ggplot(combined_data_df, aes(x = x_coord, y = y_coord, fill = predicted_shannon_gwr)) +
  geom_tile() +
  scale_fill_gradientn(
    colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"),
    limits = c(min_value, max_value)
  ) +
  coord_fixed() +
  labs(title = "GWR") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none"
  )


print(plot_gwr)


combined_data_df <- as.data.frame(combined_data)
combined_data_df$x_coord <- combined_data$x
combined_data_df$y_coord <- combined_data$y

combined_data_df$predicted_shannon_grf <- predictions_grf

print(head(combined_data_df))

plot_grf <- ggplot(combined_data_df, aes(x = x_coord, y = y_coord, fill = predicted_shannon_grf)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"),
                       limits = c(min(predictions_grf), max(predictions_grf))) +  
  coord_fixed() +
  labs(title = "GRF") +
  theme_minimal() +
  theme(
    panel.grid = element_blank(),
    plot.title = element_text(hjust = 0.5, size = 12),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.line = element_blank(),
    axis.title = element_blank(),
    legend.position = "none" 
  )
plot(plot_grf)
library(cowplot)


dummy_data <- data.frame(y = seq(1, 3, length.out = 100))

colorbar <- ggplot(dummy_data, aes(x = 1, y = y, fill = y)) +
  geom_tile(alpha = 0) + 
  scale_fill_gradientn(
    colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"),  
    limits = c(1, 3),  
    name = NULL 
  ) +
  theme_void() +  
  guides(
    fill = guide_colorbar(
      barwidth = 1,      
      barheight = 15     
    )
  ) +
  theme(
    legend.position = "right",      
    legend.title = element_blank(), 
    legend.text = element_text(size = 10)  
  ) 

plot(colorbar)



final_plot <- plot_grid(
  plot_ols,  
  plot_gwr,  
  plot_rf,  
  plot_grf,
  colorbar,  
  ncol = 4, 
  rel_widths = c(1, 1, 1.1, 0.3)
)

plot(final_plot)




colorbar <- ggplot() +
  scale_fill_gradientn(
    colors = c("#440154", "#3b528b", "#21908c", "#5dc863", "#fde725"), 
    limits = c(min_value, max_value), 
    guide = guide_colorbar(barwidth = 1, barheight = 15)
  ) +
  theme_void() +
  theme(legend.position = "right")

final_plot <- plot_grid(
  plot_grid(plot_ols, plot_gwr, ncol = 2, align = "v"),  
  plot_grid(plot_rf, plot_grf, ncol = 2, align = "v"), 
  ncol = 1,  
  rel_heights = c(1, 1) 
)

final_plot_with_colorbar <- plot_grid(
  final_plot, 
  colorbar, 
  ncol = 2, 
  rel_widths = c(0.9, 0.1)  
)

print(final_plot_with_colorbar)
ggsave("/user/ziqitang/data_ziqi/comparison_with_shared_colorbar.png", plot = final_plot_with_colorbar, width = 12, height = 8, dpi = 300)


