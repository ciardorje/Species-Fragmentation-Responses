rm(list=ls()); gc()

library(pacman) 
p_load(raster, rgdal, rgeos, gdalUtils, sp, tidyverse, parallel, 
       geosphere, sf, s2, landscapemetrics,fasterRaster, smoothr, 
       Rfast, stars, rmapshaper, doSNOW, lwgeom, sfheaders, purrr,
       exactextractr)

setwd()

rasterOptions(maxmemory = 1e+10, chunksize = 5e+08)
memory.limit(size=56000)

#Load delineated LULC map for AF in 2008
load('./Delineated_Patches/AF_Patches_2008.RData')
patch_raster <- x[[1]]
patch_polys <- x[[2]]
patch_polys <- st_transform(patch_polys, 'EPSG:32721')

colnames(patch_polys)[1] <- 'Fragment'
patch_polys$Fragment <- 1:nrow(patch_polys)

#-----Calculate Patch-Level Metrics-----
source('./Code/BIOFRAG_Delineation.R') #Load functions to derive shape and buffer metrics

#Select just patches sampled 
#Optimum buffer size will later be determined and extracted for all patches
frags <- read.csv('./Raw_Data/FragData.csv')

our_patches <- frags %>% 
  st_as_sf(coords = c("Core_Lon", "Core_Lat"), crs = 4326, remove = F) 
our_patches$Polygon <- as.numeric(st_intersects(st_transform(our_patches, 'EPSG:32721'), 
                                                patch_polys))
our_polys <- patch_polys[patch_polys$Fragment %in% our_patches$Polygon,] 

our_patches <- our_patches[order(our_patches$Polygon),]
our_polys <- our_polys[order(our_polys$Fragment),]
our_polys$Fragment <- our_patches$Fragment

#Find area of fragments
patch_metrics <- data.frame(Fragment = our_polys$Fragment, 
                            CF = our_patches$Continuous_Forest)
patch_metrics$Area_m2 <- as.numeric(st_area(our_polys))
biggest_patch <- max(patch_metrics$Area_m2[patch_metrics$CF == 0]) #Find largest non-continuous forest patch +
CF_area <- biggest_patch * 10 #derive assigned CF area value (1 order magnitude larger than largest sampled patch)
patch_metrics$Area_m2[patch_metrics$Area_m2 > biggest_patch] <- CF_area
patch_metrics$Area_Ha <- patch_metrics$Area_m2/10000 #Convert m2 to Ha

#Find shape (compactness) of fragments
for(i in as.numeric(row.names(patch_metrics[patch_metrics$CF == 0,]))){
  patch_metrics$Shape[i] <- compactness(our_polys[i,])
}
patch_metrics$Shape[patch_metrics$CF == 1] <- 1

#Find proximity (metric) of fragments to other forest patches
buffers <- c(100, 250, 500, 750, 1000, 2500) #buffer sizes

#Proximity Index Functions
#Function to create buffer around each fragment and remove patch area (only retain buffer)
outer_buffer <- function(x, size){
  x <- st_make_valid(x)
  y <- st_buffer(x, size) %>% st_make_valid()
  z <- st_difference(x = y, y = x)
  return(z)
}

#FRAGSTATS proximity metric 
prox_index <- function(focal, buff_size, polygons, CF_area){
  
  #Find all neighbouring patches within X m
  focal_id <- focal$Fragment
  dists <- gDistance(polygons, focal, byid = T)
  polygons$Distance <- as.numeric(dists)
  ngbs <- polygons[polygons$Distance <= buff_size,]
  ngbs <- ngbs[ngbs$Fragment != focal_id,]
  
  if(!is_empty(ngbs)){
    
    area <- as.numeric(area(ngbs))
    ngbs$Distance[ngbs$Distance == 0] <- 1
    prox <- ngbs$Area/(ngbs$Distance^2) 
    prox <- sum(prox)
    
  } else { prox <- 0 }
  
  return(prox)
  
}


#Add prox columns to results df
buffer_metrics_df <- as.data.frame(matrix(NA, nrow = nrow(patch_metrics), ncol = length(buffers)))
colnames(buffer_metrics_df) <- c(paste0('Prox_', buffers, 'm'))
patch_metrics <- cbind(patch_metrics, buffer_metrics_df)

patch_polys <- st_transform(patch_polys, 'EPSG:32721')

for(i in as.numeric(row.names(patch_metrics[patch_metrics$CF == 0,]))){
  cat(paste0('\n', patch_metrics$Fragment[i], ':'))
  for(n in 1:length(buffers)){
    buffer <- st_buffer(our_polys[i,], buffers[n]) %>% st_transform('EPSG:32721')
    patch_metrics[i,(n+5)] <- prox_index(our_polys[i,], buffer, patch_polys, CF_area)
    cat(paste(' ', buffers[n]))
  }
}

#Set CF site values
CF_prox <- patch_metrics[patch_metrics$Area_m2 == biggest_patch, c(6:11)] * 10
patch_metrics[patch_metrics$CF == 1, 6:11] <- CF_prox

#Find minimum distance between each sampled patch for spatial autocorrelation analyses

dist_mat <- as.data.frame(matrix(nrow = nrow(our_polys), ncol = nrow(our_polys)))
colnames(dist_mat) <- our_polys$Fragment
rownames(dist_mat) <-  our_polys$Fragment

for(i in 1:nrow(our_polys)){
  for(j in 1:nrow(our_polys)){
    
    dist_mat[i,j] <- as.numeric(st_distance(our_polys[i,], our_polys[j,]))
    
  }
}

#-----Save-----
patch_metrics <- merge(frags, patch_metrics[,-2])

write.csv(patch_metrics, 'AF_DB_Fragment_Metrics.csv')
write.csv(dist_mat, 'AF_DB_Patch_Distance_Matrix.csv', row.names = T)
save(patch_metrics, file = 'AF_2008_FragMetrics.RData')

