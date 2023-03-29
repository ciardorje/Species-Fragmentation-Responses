rm(list=ls()); gc()

library(pacman) 
p_load(raster, rgdal, rgeos, gdalUtils, sp, tidyverse, EBImage, parallel, 
       geosphere, sf, s2, leaflet, mapview, rgee, reticulate, landscapemetrics,
       fasterRaster, smoothr, Rfast, googledrive, stars, rmapshaper,
       leafem, doSNOW, lwgeom, RANN)

setwd()

rasterOptions(maxmemory = 1e+10, chunksize = 5e+08)
memory.limit(size=56000)

#-----Initialize RGEE-----
rgee_environment_dir = ".conda\\envs\\rgee_py"

reticulate::use_python(rgee_environment_dir, required=T)
rgee::ee_install_set_pyenv(
  py_path = rgee_environment_dir, 
  py_env = "rgee_py" 
)
Sys.setenv(RETICULATE_PYTHON = rgee_environment_dir)
Sys.setenv(EARTHENGINE_PYTHON = rgee_environment_dir)

# Initialize the Python Environment
rgee::ee_Initialize(drive = T)


#-----Obtain Sites and Locality-----

#Read sample sites 
load('AF_DB_Clean_Data.RData'); rm(dbs_long, dbs_wide, effort)

sites <- frags %>% st_as_sf(coords = c("Core_Lon", "Core_Lat"), crs = 4326)

#Obtain Alta Floresta boundaries from online
  #download.file('http://stacks.stanford.edu/file/druid:cs711fk6471/data.zip', 
  #destfile = './Mato_Grosso_Municipalities.zip')
  #unzip('Mato_Grosso_Municipalities.zip')
MG <- st_read('51MUE250GC_SIR.shp')
AF <- MG[MG$NM_MUNICIP == 'ALTA FLORESTA',]
AF_area <- st_bbox(AF) %>% 
  st_as_sfc() %>% 
  st_as_sf() %>%
  st_buffer(dist = 5000) %>%
  st_bbox() %>% as.numeric()

#Convert to GEE objects
sites <- sf_as_ee(sites)
area <- ee$Geometry$Polygon(
  list(
    c(AF_area[1], AF_area[2]),
    c(AF_area[1], AF_area[4]),
    c(AF_area[3], AF_area[4]),
    c(AF_area[3], AF_area[2])
  )
)


#-----Get MapBiomas LandCover maps-----

mapbiomas <- ee$Image("projects/mapbiomas-workspace/public/collection6/mapbiomas_collection60_integration_v1")
  
#Mapbiomas color palette
idx_mapbiomas = c(1, 3, 4, 5, 49, 10, 11, 12, 32, 29, 13, 14, 15, 
                  18, 19, 39, 20, 40, 41,  36,  46,  47,  48,  9, 
                  21,  22,  23,  24,  30,  25,  26,  33,  31,  27) 
hex_mapbiomas = c("#129912", "#006400", "#00ff00", "#687537", "#6b9932", 
                  "#BBFCAC", "#45C2A5", "#B8AF4F", "#968c46", "#665a3a", 
                  "#f1c232", "#FFFFB2", "#FFD966", "#E974ED", "#D5A6BD", 
                  "#e075ad", "#C27BA0", "#982c9e", "#e787f8", "#f3b4f1",
                  "#cca0d4", "#d082de", "#cd49e4", "#ad4413", "#fff3bf", 
                  "#EA9999", "#DD7E6B", "#aa0000", "#af2a2a", "#ff3d3d", 
                  "#0000FF", "#0000FF", "#02106f", "#D5D5E5")
mapbiomas_palette = rep("#FFFFFF", 49)
mapbiomas_palette[idx_mapbiomas] = hex_mapbiomas

vis_mapbiomas = list(bands = c('classification_2008'), 
                     min = 1, max = 49, palette = mapbiomas_palette)

#Map
Map$centerObject(eeObject = area, 8)
Map$addLayer(eeObject = mapbiomas, 
             visParams = vis_mapbiomas, 
             'MapBiomas') +
  Map$addLayer(eeObject = area, 
               visParams = list(strokeWidth = 3,
                                opacity = 0), 
               'Sample_Area') +
  Map$addLayer(eeObject = sites, 
               visParams = list(pointRadius = 5, 
                                color = "FF0000"), 
               'Sampled_Fragments')

#-----Crop Mapbiomas Data to Alta Floresta-----

#Convert AF bounds to feature 
AF_ft <- ee$Feature(area)
AF_ft$getInfo()

#Crop Mabiomas
af_lc <- mapbiomas$clip(AF_ft)
ee_image_info(af_lc)

#Check crop
Map$centerObject(eeObject = area, 8)
Map$addLayer(eeObject = af_lc, 
             visParams = vis_mapbiomas, 
             'AF_MapBiomas') +
  Map$addLayer(eeObject = area, 
               visParams = list(strokeWidth = 3,
                                opacity = 0), 
               'Sample_Area')


#-----Fork to Drive and Download to Local Disk----- 

upload <- ee_image_to_drive(image = af_lc,
                            maxPixels = 1e12,
                            folder = "GEE",
                            description = paste0('Mapbiomas_Alta_Floresta_LU'),
                            fileNamePrefix = paste0('Mapbiomas_Alta_Floresta_LU'),
                            skipEmptyTiles = T,
                            timePrefix = FALSE)
upload$start()
ee_monitoring()

fname = ee_drive_to_local(task = upload)
file.rename(fname,
            paste0('./MapBiomas_Data/', basename(fname)))


#----Load and Recategorise Rasters-----

af_lc <- stack('./MapBiomas_Data/Mapbiomas_Alta_Floresta_LU.tif')
names(af_lc) <- paste0('X', 1985:2020)
names(af_lc@layers) <- paste0('X', 1985:2020)

#codes for MAPBIOMAS forest categories
forest <- c(1,3,4,5,49) 

#Raster compression options
tifoptions <- c("COMPRESS=DEFLATE", "PREDICTOR=2", "ZLEVEL=6") 

#Binarise and save
for(i in 1:nlayers(af_lc)){
  
  lyr <- af_lc[[i]]
  lyr[lyr %in% forest] <- 1
  lyr[!(lyr %in% forest)] <- 0
  
  writeRaster(lyr, 
              filename = paste0('./MapBiomas_Data/Binary/AF_Binary_LC_', names(lyr), '.tif'),
              format = 'GTiff', options = tifoptions, progress = 'text', overwrite = T)
  
  print(names(lyr))
  if(i == nlayers(af_lc)){rm(af_lc)}
  
}

rm(list=ls()); gc()

#-----Delineate Patches-----

source('./Code/BIOFRAG_Delineation.R')

#Read binary LU rasters 
af_lc <- stack(list.files('./MapBiomas_Data/Binary/', full.names = T))
years <- 1985:2020

#Apply watershed delineation and save
for(i in 1:nlayers(af_lc)){
  
  x <- delineate(af_lc[[i]], res_m = 30, edge_depth = 90, 
                gap_area = NA, corridor_width = 180)
  
  save(x, file = paste0('./Delineated_Patches/AF_Patches_', years[i], '.RData'))
  print(paste('\n', years[i]))
  gc()

}

