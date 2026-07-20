---
  title: "Solar PV Suitability Analysis: User Guide"
date: 2026-07-15
author: Kenzie Thomson

format:
  html:
  toc: true
toc-location: left
toc-depth: 4
theme: minty
title-block-banner: images/solarfarm.jpg
title-block-banner-color: white
---
  
  # About This Project
  
  As electricity demand grows across British Columbia, more renewable energy installations need to be developed to support a resilient, diversified, low-carbon future. In light of this challenge, a site suitability analysis was undergone to create a province-wide screening tool that determines the most optimal locations for solar PV development.

This analysis is broken down into 2 main sections:
  
  **Section 1** determines the geographical potential of BC - the total land area where solar farm development is possible.

**Section 2** determines the technical potential of BC - the total energy these sites could produce on a yearly basis.

By providing clear and accessible information on renewable energy potential, this project helps accelerate clean energy development, informs decision-making and supports BC's transformation to a more resilient and sustainable electricity system.

This document provides a detailed guide of the analysis, including all code used to derive results.

[Click for more info](https://www.uvic.ca/acet/projects/bc-renewable-energy-mapping/index.php)

[Click for first-draft publication](https://borealisdata.ca/dataset.xhtml?persistentId=doi:10.5683/SP3/MZYJFF)

[Click for raw data download](tbd)

**NOTE:** You will need to change each file path to a folder on your computer!

# Set-Up

First, each package must be installed and loaded to complete this analysis

```{r}
#| output: false
#| message: false

# Load packages
library(terra)
library(sf)
library(tidyverse)
library(exactextractr)
library(AHPtools) 
library(ggplot2)
library(scales) 
library(dplyr)
library(classInt)
library(bcmaps)
```

I recommend setting up an organized file system, similar to the one below, to save your intermediate outputs and final results

/Solar_Suitability
|- /raw_data 
|- /criteria
  |- /cleaned # for data layers that are reprojected, clipped, masked, resampled, and rasterized 
  |- /reclassified # for data layers that are reclassified to a 1-5 suitability scale
  |- /slope # to save the dem derived slope pieces 
  |- /aspect # to save the dem derived aspect pieces 
|- /constraints 
  |- /binary # for data layers brought to a binary classification 
  |- /mask # for the final constraints mask 
|- /suitability_surfaces
  |-/unmasked # suitability results without constraints mask 
  |- /masked # suitability results with constraints mask

# Step 1: Data Pre-Processing

## 1.1 Create Template

To first step of the analysis is the creation of a template raster. This layer will ensure that every input layer is the exact same extent and resolution. We will do create this by rasterizing a polygon representing BCs administrative boundary. 

```{r}
#| message: false
#| output: false
#| eval: false

# Create template  

# Read in layer
bc <- sf::st_read("~/Library/CloudStorage/OneDrive-UBC/RA/Data/Admin Boundaries/BC_RegionalDistricts.shp") 

# Project to BC Enviro Albers (EPSG:3005)
bc <- sf::st_transform(bc, crs = 3005)

# Create template spanning vector extent
# Set 30 m spatial resolution (each cell = 30m x 30m) 
template <- rast(bc, res = 30)

# Convert to raster and burning in the Regional District Name attribute column 
template <- rasterize(bc, template, field = "NAME_2")

# Save the output
writeRaster(template, "~/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/template/template.tif", overwrite = TRUE)

# Quick read
template <- rast("~/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/template/template.tif")
```

The next step of the analysis is data pre-processing. This section is divided in 2 parts - suitability criteria and land constraints.

## 1.2 Suitability Criteria Processing

There are 6 suitability criteria layers used in this analysis:

1.  Global Tilted Irradiance (GTI)

2.  Ground Slope

3.  Ground Aspect

4.  Proximity to Substations

5.  Proximity to Demand Centers

6.  Proximity to Major Roads

Each layer needs to be clipped to BC's bounding box, masked to a template raster, projected to a common CRS (EPSG:3005 - BC Environmental Albers), sampled to a consistent spatial resolution, and reclassified to a 1-5 suitability scale. We will process each layer individually to ease computational power.

### 1.2.1 Global Tilted Irradiance

Data source: [Global Solar Atlas GIS Data](https://globalsolaratlas.info/download/canada) 

Clean data:
  
  ```{r}
#| message: false 
#| eval: false

# Read in GTI raster layer 
GTI <- rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Existing Renewables Mapping Project/Raw solar and wind/Canada_GISdata_LTAy_YearlyMonthlyTotals_GlobalSolarAtlas-v2_GEOTIFF/GTI.tif")

# Project to template raster projection (EPSG:3005)
# Resampling method = bilinear interpolation (used for continuous data)
gti_proj <- project(GTI, template, method = "bilinear")

# Clip to BC bounding box
gti_crop <- crop(gti_proj, template)

# Mask to BC administrative boundary
gti_mask <- mask(gti_crop, template) 

# Save in cleaned data folder
writeRaster(gti_mask, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/gti_cleaned.tif", overwrite = TRUE)
```

Reclassify to 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using natural breaks (jenks) 

# Get data bounds
true_min <- global(gti_mask, "min", na.rm = TRUE)[1, 1]
true_max <- global(gti_mask, "max", na.rm = TRUE)[1, 1]

set.seed(123)
vals <- terra::spatSample(gti_mask, 500000, method = "regular", na.rm = TRUE)
jenks_obj <- classIntervals(as.numeric(vals[[1]]), n = 5, style = "fisher") 
breaks <- jenks_obj$brks

# Set data bounds 
breaks[1] <- true_min
breaks[length(breaks)] <- true_max 

# Create classification matrix
m_jenks <- matrix(c(
  head(breaks, -1),
  tail(breaks, -1),
  1:5
), ncol = 3)

# Execute classification 
gti_reclass <- terra::classify(gti_mask, m_jenks, include.lowest = TRUE)

# Save in reclassified data folder
writeRaster(gti_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/gti_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 1 - GTI 
#| 
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/gti_reclassified.tif"), col=suitability_cols)
```


### 1.2.2 Ground Slope

Data source: [bcmaps R Package](https://cran.r-project.org/web/packages/bcmaps/index.html) 

The ground slope raster is created by piecing together a DEM from the bcmaps package and then utilizing terrain functions to calculate slope. Initiate this for-loop to download the DEM, create the slope raster, and mask to our template raster. Expect this code chunk to take a few hours to finish. 

```{r}
#| message: false
#| eval: false

# Check available layers in bcmaps package
bc_maps_avail <- available_layers()

# Download DEM by ecoregion, check layer names
ecoprov <- ecoprovinces() 

# Set output directory
out_dir <- "Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/slope"

# Define ecoprovinces to process in a list
ecoprovs <- c(
  "SAL", "NBM", "TAP", "BOP", "SBI",
  "SIM", "SOI", "COM", "GED", "NEP", "CEI"
)

# Call template raster mask template
template <- rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/template/template.tif") 

# Initiate for-loop to create slope raster from DEM chunks 
for (prov in ecoprovs) {
  
  message("Processing ", prov)
  
  # Get ecoprov boundary
  bound_sf <- ecoprovinces() %>%
    filter(ECOPROVINCE_CODE == prov)
  
  # Convert to SpatVector and project to BC Albers (EPSG:3005)
  bound_vect <- vect(bound_sf)
  bound_vect <- project(bound_vect, "EPSG:3005")
  
  # BUFFER: Create a 500m buffer around the boundary to pull neighboring DEM data
  bound_buffered <- buffer(bound_vect, width = 500) 
  
  # Get DEM using the buffered boundary 
  dem <- cded_raster(sf::st_as_sf(bound_buffered))
  
  # Project DEM to template grid
  dem_3005 <- project(
    rast(dem),
    template,
    method = "bilinear")
  
  # Calculate slope on the buffered DEM first to ensure seamless edges
  slope <- terrain(
    dem_3005,
    v = "slope",
    neighbors = 8,
    unit = "degrees")
  
  # Crop and mask to the strict, original ecoprovince boundary
  slope_clipped <- crop(slope, bound_vect)
  slope_clipped <- mask(slope_clipped, bound_vect)
  
  # Write slope raster
  writeRaster(slope_clipped, file.path(out_dir, paste0(prov, "_slope.tif")), overwrite = TRUE)
  
  # Clean memory
  rm(dem, dem_3005, slope, slope_clipped, bound_vect, bound_buffered)
  gc()
}

# Now join all individual ecoprovince slope rasters together 

# List and load slope rasters
slope_files <- list.files(path = out_dir, pattern = "_slope\\.tif$", full.names = TRUE)
slope_list <- lapply(slope_files, rast)

# Merge all slopes into BC-wide raster
bc_slope <- do.call(merge, slope_list)

# Force exact match to template extent, then mask
bc_slope_aligned <- crop(bc_slope, template)
bc_slope_final <- mask(bc_slope_aligned, template)

# Save in cleaned data folder
writeRaster(bc_slope_final, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/bc_slope_cleaned.tif")
```

Reclassify to a 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using a defined matrix

# Create classes using format (from, to, new)
m <- c(0, 5, 5, #best
       5, 10, 4, 
       10, 15, 3, 
       15, 20, 2, 
       20, 100, 1) #worst 

# Create matrix
slope_matrix <- matrix(m, ncol = 3, byrow = TRUE)

# Apply reclassification 
slope_reclass <- terra::classify(bc_slope_final, slope_matrix, include.lowest = TRUE, right = FALSE)

# Save in reclassified data folder
writeRaster(slope_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/slope_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 2 - Ground Slope
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/slope_reclassified.tif"), col=suitability_cols)
```

### 1.2.3 Ground Aspect 

Data source: [bcmaps R Package](https://cran.r-project.org/web/packages/bcmaps/index.html) 

We will apply the same for-loop methodology to create a BC wide aspect layer from a DEM. Again, expect this code chunk to take a few hours to finish. 

```{r}
#| message: false
#| eval: false

# Check available layers in bcmaps package
bc_maps_avail <- available_layers()

# Download DEM by ecoregion, check layer names
ecoprov <- ecoprovinces() 

# Set output directory
out_dir <- "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/aspect"

# Define ecoprovinces to process in a list
ecoprovs <- c(
  "SAL", "NBM", "TAP", "BOP", "SBI",
  "SIM", "SOI", "COM", "GED", "NEP", "CEI"
)

# Call template raster mask template
template <- rast("~/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/template/template.tif") 

# Initiate loop
for (prov in ecoprovs) {
  
  message("Processing ", prov)
  
  # Get ecoprov boundary
  bound_sf <- ecoprovinces() %>%
    filter(ECOPROVINCE_CODE == prov)
  
  bound_vect <- vect(bound_sf)
  bound_vect <- project(bound_vect, "EPSG:3005")
  
  # BUFFER: Create a 500m buffer to pull overlapping DEM data
  bound_buffered <- buffer(bound_vect, width = 500)
  
  # Get DEM using the buffered boundary
  dem <- cded_raster(sf::st_as_sf(bound_buffered))
  
  # Project DEM to template
  dem_3005 <- project(
    rast(dem),
    template,
    method = "bilinear")
  
  # Calculate aspect on the buffered DEM
  aspect <- terrain(
    dem_3005,
    v = "aspect",
    neighbors = 8,
    unit = "degrees")
  
  # Crop and mask strictly to the original ecoprovince boundary
  aspect_clipped <- crop(aspect, bound_vect)
  aspect_clipped <- mask(aspect_clipped, bound_vect)
  
  # Write aspect raster
  writeRaster(
    aspect_clipped,
    file.path(out_dir, paste0(prov, "_aspect.tif")),
    overwrite = TRUE)
  
  # Clean memory
  rm(dem, dem_3005, aspect, aspect_clipped, bound_vect, bound_buffered)
  gc()
}

# Now join all individual ecoprovince aspect rasters together 

# List all aspect rasters
aspect_files <- list.files(
  path = out_dir,
  pattern = "_aspect\\.tif$",
  full.names = TRUE)

# Load rasters into a list
aspect_list <- lapply(aspect_files, rast)

# Merge into one surface
bc_aspect <- do.call(merge, aspect_list)

# Crop and mask to template 
bc_aspect_aligned <- crop(bc_aspect, template)
bc_aspect_final <- mask(bc_aspect_aligned, template)

# Save in cleaned data folder
writeRaster(bc_aspect_final, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/bc_aspect_cleaned.tif", overwrite = TRUE)
```

Reclassify to a 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using a defined matrix

# Create classes using format (from, to, new)
m <- c(315, 360, 1, #N
       0, 45, 1, #N (wrap-around)
       45, 135, 3, #E
       135, 225, 5, #S
       225, 315, 3) #W

# Create matrix
aspect_matrix <- matrix(m, ncol = 3, byrow = TRUE)

# Apply classification
aspect_reclass <- terra::classify(bc_aspect_final, aspect_matrix, include.lowest = TRUE, right = TRUE)

# Save in reclassified data folder
writeRaster(aspect_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/aspect_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 3 - Ground Aspect
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/aspect_reclassified.tif"), col=suitability_cols)
```

### 1.2.4 Proximity to Substations

Data source: Private 

Substation locations are point data, therefor, we must first create a euclidean distance raster where each cell represents a distance in meters away from the closest substation. 

```{r}
#| message: false
#| eval: false

# Load in vector layer
subs <- sf::st_read("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Criteria/Raw/Substations/substations.shp")

# Transform to spatvect for processing
subs_v <- vect(subs)

# Create Euclidean distance raster showing proximity to substations that spatially conforms to our template
subs_dist <- terra::distance(template, subs_v, rasterize = TRUE)

# Mask to template boundary
subs_mask <- mask(subs_dist, template)

# Save in cleaned data folder
writeRaster(subs_mask, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/subs_cleaned.tif", overwrite = TRUE)
```

Reclassify to a 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using a defined matrix

# Units = meters because our crs = projected to BC Albers
m <- c(0, 1000, 5, # 1km surrounding substations is best
       1000, 2000, 4,
       2000, 5000, 3,
       5000, 7500, 2,
       7500, Inf, 1) # Anything beyond 7.5km is worst 

# Create matrix
subs_m <- matrix(m, ncol = 3, byrow = TRUE)

# Apply classification
subs_reclass <- terra::classify(subs_mask, subs_m, include.lowest = TRUE, right = FALSE)

# Save in reclassified data folder
writeRaster(subs_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/subs_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 4 - Proximity to Substations
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/subs_reclassified.tif"), col=suitability_cols)
```

### 1.2.5 Proximity to Demand Centers

Data source: [BC Data Catalogue](https://catalogue.data.gov.bc.ca/dataset/bc-major-cities-points-1-2-000-000-digital-baseline-mapping) 

Repeat process using point locations of major cities.

```{r}
#| message: false
#| eval: false

# Load in vector layer
cities <- sf::st_read("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Criteria/Raw/DBM_BC_7H_MIL_POPULATION_POINT/DBMBC7HML4_point.shp")

# Transform to spatvect for processing
cities_v <- vect(cities)

# Create Euclidean distance raster showing proximity to substations that spatially conforms to our template
cities_dist <- terra::distance(template, cities_v, rasterize = TRUE)

# Mask to template boundary
cities_mask <- mask(cities_dist, template)

# Save in cleaned data folder
writeRaster(cities_mask, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/cities_cleaned.tif", overwrite = TRUE)
```

Reclassify to a 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using a defined matrix

# Units = meters because our crs = projected to BC Albers
m <- c(0, 5000, 5, 
       5000, 15000, 4,
       15000, 30000, 3,
       30000, 45000, 2,
       45000, Inf, 1)  

# Create matrix
cities_m <- matrix(m, ncol = 3, byrow = TRUE)

# Apply classification
cities_reclass <- terra::classify(cities_mask, cities_m, include.lowest = TRUE, right = FALSE)

# Save in reclassified data folder
writeRaster(cities_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/cities_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 5 - Proximity to Demand Centers
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/cities_reclassified.tif"), col=suitability_cols)
```
### 1.2.6 Proximity to Major Roads

Data source: [BC Data Catalogue](https://catalogue.data.gov.bc.ca/dataset/7-5m-major-roads-the-atlas-of-canada-base-maps)

Repeat process using locations of major roads.

```{r}
#| message: false
#| eval: false

# Load in vector layer 
roads <- sf::st_read("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Criteria/Raw/BCGW_02001F02_1784319511692_11424/DBM_BC_7H_MIL_ROADS_LINE.gdb")

# Transform to spatvect for processing
roads_v <- vect(roads)

# Create Euclidean distance raster showing proximity to substations that spatially conforms to our template
roads_dist <- terra::distance(template, roads_v, rasterize = TRUE)

# Mask to template boundary
roads_mask <- mask(roads_dist, template)

# Save in cleaned data folder
writeRaster(roads_mask, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/cleaned/roads_cleaned.tif", overwrite = TRUE)
```

Reclassify to a 1-5 suitability scale:
  
  ```{r}
#| message: false
#| eval: false

# Reclassify using a defined matrix

# Units = meters because our crs = projected to BC Albers
m <- c(0, 1000, 5, 
       1000, 5000, 4,
       5000, 15000, 3,
       15000, 30000, 2,
       30000, Inf, 1)  

# Create matrix
roads_m <- matrix(m, ncol = 3, byrow = TRUE)

# Apply classification
roads_reclass <- terra::classify(roads_mask, roads_m, include.lowest = TRUE, right = FALSE)

# Save in reclassified data folder
writeRaster(roads_reclass, "/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/roads_reclassified.tif", overwrite = TRUE)
```

```{r}
#| echo: false
#| title: Visualize Criteria 6 - Proximity to Major Roads
suitability_cols <- c("#FFFFCC", "#FECC5C","#FD8D3C","#F03B20","#BD0026")
plot(rast("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/SolarModel_Guide/data/criteria/reclassified/roads_reclassified.tif"), col=suitability_cols)
```

## 1.3 Land Constraints Processing

This portion of the analysis creates the land constraints spatial mask used to eliminate land restricted from solar farm development due to technical, environmental, or regulatory factors. 

Land constraints include:
  
  1. Land covers - water, permanent snow and ice, wetlands, artificial and built environments 

2. All roads - including forest service roads and names roads

3. Densely forested areas (canopy density > 40%)

4. ALR land with a soil capability class of 1 and 2 (most arable land)

5. Protected and conserved areas

6. Archaeologically and culturally sensitive areas

Each layer will be transformed to a binary with 0 representing land to be included in the constraints mask and 1 representing land that is not included in the constraints mask and is therefor available for inclusion in the suitability model. Each layers '0' cells will then be merged to create a master restriction zone layer. 

### 1.3.1 Land Cover



```{r}
# View available layers in BC Digital Road Atlas geodatabase
st_layers("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Criteria/Raw/dgtl_road_atlas.gdb")

# Load in vector layer
roads <- sf::st_read("/Users/kenziethomson/Library/CloudStorage/OneDrive-UBC/RA/Data/Criteria/Raw/dgtl_road_atlas.gdb")
```


