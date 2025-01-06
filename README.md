# Spatial Patterns of Species Diversity in the Cape Floristic Region

This repository contains the implementation and analysis of spatial patterns of species diversity (Shannon Index) within the Cape Floristic Region, South Africa, using spectral variation analysis. The study evaluates the relationships between environmental factors (NDVI, temperature, elevation, and aspect) and alpha diversity through four different models:

1. **Ordinary Least Squares (OLS)**
2. **Geographically Weighted Regression (GWR)**
3. **Random Forest (RF)**
4. **Geographically Random Forest (GRF)**

## Data Description  

### Satellite Imagery  
- **Source**: EMIT L2A product  
- **Acquisition Date**: January 19, 2023  
- **Spatial Resolution**: 60x60 meters  
- **Description**: Captured using the EMIT imaging spectrometer, installed on the International Space Station in 2022, to collect hyperspectral data across Earth's sunlit regions.  

### Digital Elevation Model (DEM)  
- **Source**: Copernicus satellite data  
- **Description**: Generated to represent terrain elevation.  

### Topographic Position Index (TPI)  
- **Calculation**: Derived from elevation data using QGIS software.  

### Mean Daily Air Temperature (TAS)  
- **Source**: CHELSA database (Climatologies at High Resolution for Earth's Land Surface Areas)  
- **Description**: Represents the mean daily air temperature across the study region.  

### Fire Age Dataset  
- **Description**: Records the time elapsed since the most recent fire event for each spatial pixel.  

### Normalized Difference Vegetation Index (NDVI)  
- **Source**: MODIS NDVI data  
- **Description**: Indicates vegetation health and coverage using remote sensing data.  

### Aspect Data  
- **Derived From**: Digital Elevation Model (DEM)  
- **Description**: Represents the slope direction (e.g., north, south) based on elevation.  

## Key Findings
- Coastal areas show higher alpha diversity due to favorable environmental conditions.
- Higher elevations correlate with reduced species diversity.
- NDVI exhibits a negative correlation with alpha diversity, with regional variations.
- GWR outperforms other models by capturing spatial heterogeneity, while GRF integrates spatial dependencies for robust predictive performance.

This work demonstrates the value of geographical statistical methods in biodiversity monitoring, highlighting the importance of spatial analysis in understanding ecological patterns.


## Maps Of Results

<div align="center">
  <img width="345" alt="f2d62c74bc173a589aba04019844eab" src="https://github.com/user-attachments/assets/fb7e60fe-2330-4a7f-a7e2-2018560baf90" /><br>
  <img width="338" alt="4767d9afbb230a64125dada361612c6" src="https://github.com/user-attachments/assets/a4745bd9-ea65-4ca3-80bf-0cab114b4a15" />
</div>

## References

1. Féret, J.-B., de Boissieu, F. (2019). [biodivMapR: an R package for α‐ and β‐diversity mapping using remotely‐sensed images](https://doi.org/10.1111/2041-210X.13310). *Methods in Ecology and Evolution*, 00:1-7.
   
2. Féret, J.-B., Asner, G.P. (2014). [Mapping tropical forest canopy diversity using high-fidelity imaging spectroscopy](https://doi.org/10.1890/13-1824.1). *Ecological Applications*, 24, 1289–1296.

