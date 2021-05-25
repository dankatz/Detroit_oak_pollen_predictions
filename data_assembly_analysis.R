# Manuscript: Modeling airborne pollen at municipal scales
# This is the main script for the oak pollen manuscript, and includes: data assembly and analysis 

### data assembly: airborne pollen ######################################################################

# set up work environment
library(lubridate)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(sf)
library(raster)

## data assembly: spring 2017 data =====================================================================================
# Airborne pollen data was collected spring of 2017, and is described in Katz et al. 2019:
# Effect of intra-urban temperature variation on tree flowering phenology, airborne pollen, and 
# measurement error in epidemiological studies of allergenic pollen
p17 <- read_csv("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pollen_2017.csv") %>% 
    dplyr::select(date, jday, long, lat, oak_mean, oak_sd) %>% 
    mutate(dates = ymd(date),
           taxa = "p_Quercus",
           taxa2 = "Quercus",
           taxa3 = "italic(Quercus)")
p17_utm <- st_as_sf(p17, coords = c("long", "lat"), crs = 4326)
p17_utm <- st_transform(p17_utm, 32617)

##add in SCS data that are available for 2017
scs <- read.csv("C:/Users/dsk856/Box/MIpostdoc/LakeshoreENT_pollen_counts_2009_2017_compiled.csv") %>% 
  mutate(dates = mdy(Date),
         years = year(dates),
         jday = yday(dates)) #names(scs)
scs_c <- filter(scs, years == 2017) %>%
  dplyr::select(dates, Quercus, jday) %>%
  rename(scs_oak = Quercus)

#p2017_sf_utm17 <- left_join(p2017_sf_utm17, scs_c)


## data assembly: spring 2018 data =====================================================================================
# Airborne pollen data was collected spring of 2018, and is described in Katz and Batterman 2020:
# Urban-scale variation in pollen concentrations: 
# a single station is insufficient to characterize daily exposure

p18 <- read_csv("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/airborne_p_2018_200228.csv")  %>% #from pollen_heterogeneity_200106.R
       filter(taxa2 == "Quercus") %>% 
       mutate(dates = ymd(date))
p18_utm <- st_as_sf(p18, coords = c("long", "lat"), crs = 4326) 
p18_utm <- st_transform(p18_utm, 32617) #set as UTM 17



### data assembly: tree location and pollen production map ##############################################
# a map of oak trees in the central portion of Detroit is described in Katz et al. 2020:
# Improved Classification of Urban Trees Using a Widespread Multi-Temporal Aerial Image Dataset
# original shapefile: 

# tree_pred <- st_read("E:/tree_classification/predictions/pred190715.shp")
# p_Quru <- filter(tree_pred, prdctd_ == "Quercus") %>%  dplyr::select(area)

# Pollen production is described in Katz et al. 2020:
# Pollen production for 13 urban North American tree species: allometric equations for tree trunk diameter and crown area

# p_Quru$predpollen <- p_Quru$area * 0.97 + 17.02 #using equation for red oak
# write_sf(p_Quru, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pred_qusp_pol_prod200310.shp")
# 
# d_rast <- raster(ncol = 145, nrow = 212) #approximately 100 x 100 m pixel size
# extent(d_rast) <- extent(tree_pred)
# 
# trees_as_points <- p_sp %>% st_cast("POINT", do_split = FALSE)
# p_rast <- rasterize(trees_as_points, d_rast, field = "predpollen", fun = sum, na.rm = TRUE)
# pixel_area <- (res(p_rast)[1] * res(p_rast)[2]) #get pixel area in m2
# p_rast <- (p_rast/ pixel_area) * 1000000 #convert to pol/m2 (was originally in millions)
# p_rast[p_rast < 0] <- 0 #make sure that none of the trees have values below 0
# crs(p_rast) <- CRS('+init=EPSG:32617')
# plot(p_rast)
# projection(p_rast)

#writeRaster(p_rast, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_200302.tif", format="GTiff")
p_rast <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_200302.tif")
#plot(p_rast)

### data assembly: phenology model data #################################################################
# oak phenology data was collected spring of 2017, and is described in Katz et al. 2019:
# Effect of intra-urban temperature variation on tree flowering phenology, airborne pollen, and 
# measurement error in epidemiological studies of allergenic pollen
# see: PHENO_ANALYSES_FIGS_V4_181025.R

#create rasters of percent active flowers for oaks in 2017
#this raster was projected to UTM 17N and cell size was adjusted to match the pollen production surface in GIS
d_peak <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/peakflowering_d_UTM17c.tif") 
d_peak <- raster::resample(d_peak, p_rast)
d_peak <- crop(d_peak, p_rast)
#plot(d_peak)

# peakflow_oak_d <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/peakflowering_d_UTM17.tif") #plot(peakflow_oak_d)
# d <- extent(-83.3, -82.86, 42.30, 42.6)
# d_peak <- crop(peakflow_oak_d, d)
# plot(peakflow_oak_d)

#Add a lag to the flowering intensity model to account for differences between street trees and non-street trees
d_peak <- d_peak + 3 
d_peak[d_peak[] == -996] <- NA

##FOR ANIMATION SWAP THIS IN- IT GETS RID OF WHITE PIXELS: d_peak <- crop(peakflow_oak_d, d); plot(d_peak)
#d_peak_big <- aggregate(d_peak, 10, fun = mean, na.rm = TRUE) #plot(d_peak_big) #plot(peakflow_oak_d_water999)
#vmap$peak_oak <- raster::extract(x = peakflow_oak_d, y = cbind(vmap$long, vmap$lat), buffer = 500, fun = mean, na.rm = TRUE)

#I created an empirical look-up table of the average percent of mature flowers on days away from peak
#load in lookup table
pflow_lookup_summary <- read_csv("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/pflow_lookup_summary.csv") %>% 
  mutate(days_from_site_mean_low = days_from_site_mean - 0.5,
         days_from_site_mean_hi = pflow_lookup_summary$days_from_site_mean + 0.5)
#apply this function to create a raster for each day with the percent flowering
pflow_lookup_summary2 <- cbind(pflow_lookup_summary$days_from_site_mean_low, pflow_lookup_summary$days_from_site_mean_hi,
                               pflow_lookup_summary$mean_flow)

#creating a stack that has percent flowering on each day
flower_stack <- stack()
for(i in 110:145){
  focal_day <- i
  test2 <- focal_day - d_peak 
  test3 <- reclassify(test2, pflow_lookup_summary2)
  names(test3) <- paste("flow_", focal_day, sep = "")
  plot(test3, main = focal_day)
  flower_stack <- stack(flower_stack, test3)
}

### predictions of pollen release over space and time ###################################################
plot(p_rast)
flower_stack <- flower_stack * 0.01 #Changing percentage to proportion

# for(i in 1:20){plot(flower_stack[[i]])}
# plot(flower_stack[[20]])
p_per_day_stack <- p_rast * flower_stack



# plot(flower_stack[[1]])
# plot(p_per_day_stack[[2]])
# for(i in 1:36){plot(p_per_day_stack[[i]])}
#p_per_day_stack_test <- p_per_day_stack
#p_per_day_stack[p_per_day_stack[] == 0] <- NA #trying this so I can get white cells when no pollen prod

library(rgdal)
setwd("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/detroit cartography/detroit_water")
d_water <-readOGR(".", layer="detroit_water1")
d_water_utm <- st_as_sf(d_water)
d_water_utm <- st_transform(d_water_utm, crs(d_peak))
d_water_utm <-sf:::as_Spatial(d_water_utm)

setwd("C:/Users/dsk856/Box/MIpostdoc/Detroit spatial data/Detroit boundary/Detroit_boundary_shapefile")
d_bound <-readOGR(".", layer="POLYGON")
d_bound_utm <- st_as_sf(d_bound)
d_bound_utm <- st_transform(d_bound_utm, crs(d_peak))
d_bound_utm <-sf:::as_Spatial(d_bound_utm)




p_per_day_stack[is.na(p_per_day_stack)] <- 0 #set na values to 0 (they weren't really NA)
p_per_day_stack <- mask(p_per_day_stack,  d_bound_utm) #set all values outside of D to NA
# plot(p_per_day_stack[[10]])
# plot(mask(p_per_day_stack[[1]], d_bound_utm))

library('viridis') #throws an error, but still works fine
# colr <- colorRampPalette(brewer.pal(10, 'RdYlGn'))
# myColorkey <- list( at= c(0, 100, 1000, 10000, 100000), ## where the colors change #max(p_per_day_stack)
#                    labels=c("0", "100", "1,000", "10,000", "100,000"), ## where to print labels
#                    title = "pollen grains/day")
# 
# levelplot(x = p_per_day_stack, layers = 36, margin=FALSE, main = "Julian day 125", #cuts = 10, #at = seq(0, 100, by = 20),
#           colorkey = myColorkey, col.regions = viridis(25), background = "black") + #par.settings=list(panel.background=list(col="black")) +
#   layer(sp.polygons(d_bound_utm, lwd = 1)) + layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue")) 
#   


my.at=c(0,10, 50,100,500, 1000,5000, 10000,50000, 100000, 500000)
my.brks=seq(0, 100, by=10)

myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=c(0,10, 50,100,500, "1,000","5,000", "10,000","50,000", "100,000", "500,000")), space="bottom")
levelplot(p_per_day_stack, layers = 13, margin=FALSE, main = "Julian day 125", at= my.at, 
          colorkey=myColorkey, col.regions = viridis(25, option = "A", direction = -1), 
          background = "black", ylab = "", xlab = "",
          scales=list(x=list(draw=FALSE), y=list(draw=FALSE)))


#ANIMATE
#saving a png file for each graph
setwd("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/animation")
p_per_day_stack
#month and day as.Date((109 + 1), origin = "2016-12-31")
for(i in 1:36){ #1:36
  date_title <- paste("pollen release on:", as.Date((109 + i), origin = "2016-12-31")) #convert Julian to date
  png(filename=paste("D_qf_",109 + i,".png", sep = ""), width = 600, height = 800, units = "px")
  print(levelplot(x = p_per_day_stack, layers = i, margin=FALSE, main = date_title, cuts = 10, at = my.at,
                  colorkey = myColorkey, 
                  col.regions = viridis(25, option = "A", direction = -1),
                  ylab = "", xlab = "", scales=list(x=list(draw=FALSE), y=list(draw=FALSE))) + #
          layer(sp.polygons(d_bound_utm, lwd = 1))) +
    layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue"))
  
  dev.off()
}

#creating an animation

list.files(path = ".", pattern = "*.png", full.names = T) %>%
  purrr::map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=1) %>% # animates, can opt for number of loops
  image_write("oak_pollen_prod_in_Detroit_v1.gif") # write to current dir


## create a version where each pixel is the average of all pixels within 1km of it ------------------------------
fwModel <- focalWeight(p_per_day_stack, 1000, type='circle')
fwModel[fwModel>0] <- 1
# km_mean_fun <- function(x){focal(x, w=fwModel, fun=mean, na.rm = TRUE)}
# plot(km_mean_fun(p_per_day_stack[[2]]))

#I'm having a hard time applying this function to the whole raster brick, so I'm just doing it with a loop
pol_prod_day_1km <- stack() # create an empty stack
for (i in 1:nlayers(p_per_day_stack)){
  mean_pollen_prod1km <- focal(p_per_day_stack[[i]], w=fwModel, fun=mean, na.rm = TRUE)
  pol_prod_day_1km <- stack(pol_prod_day_1km , mean_pollen_prod1km )
}#plot(x[[10]])

pol_prod_day_1km_m <- pol_prod_day_1km
pol_prod_day_1km_m[is.na(pol_prod_day_1km_m)] <- 0 #set na values to 0 (they weren't really NA)
pol_prod_day_1km_m <- mask(pol_prod_day_1km_m,  d_bound_utm) #set all values outside of D to NA


#month and day as.Date((109 + 1), origin = "2016-12-31")
setwd("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/animation_prod__per_day_within1km")

my.at=c(0, 100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 30000,40000, 60000)
my.brks=seq(0, 100, by=10)
myColorkey <- list(at=my.brks, 
                   labels=list(at=my.brks,
                               labels= my.at, #c(0,10, 50,100,500, "1,000","5,000", "10,000","50,000", "100,000", "500,000")
                               #space="right", 
                               #title = "pollen release within 1 km (grains/m2)"
                   ))


for(i in 1:36){ #1:36
  date_title <- paste("pollen release on:", as.Date((109 + i), origin = "2016-12-31")) #convert Julian to date
  png(filename=paste("D_qf_",109 + i,".png", sep = ""), width = 600, height = 800, units = "px")
  print(levelplot(x = pol_prod_day_1km_m, layers = i, margin=FALSE, xlim=c(320420,333200), ylim=c(4680362,4700550),
                  main = date_title, cuts = 500, at = my.at, 
                  colorkey = myColorkey,  
                  col.regions = viridis(500, option = "A", direction = -1),
                  ylab = "", xlab = "", scales=list(x=list(draw=FALSE), y=list(draw=FALSE))) + #
          layer(sp.polygons(d_bound_utm, lwd = 1)) 
        #layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue"))
  )
  dev.off()
}


list.files(path = ".", pattern = "*.png", full.names = T) %>%
  purrr::map(image_read) %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=1) %>% # animates, can opt for number of loops
  image_write("oak_pollen_release_within_1km_v1.gif") # write to current dir


### Fig 2: comparing pollen release over space with airborne pollen concentrations ######################
#comparing 2017 data

#comparing 2018 data



### Fig 3: comparing pollen release over space and time with airborne pollen concentrations #############
p2017_sf_utm17$p_prod_day_1km <- NA
for(i in 1:nrow(p2017_sf_utm17)){
  date_obs <- p2017_sf_utm17$jday[i]
  p2017_sf_utm17_focal_obs <- p2017_sf_utm17[i,]
  
  #pol_prod_day_1km_m has 36 layers for sequential days. layer 1 == jday 110, layer 2 == jday 111, ...
  relevant_layer <- p2017_sf_utm17_focal_obs$jday - 109  #select the right layer
  p2017_sf_utm17$p_prod_day_1km[i] <- raster::extract(pol_prod_day_1km_m[[relevant_layer]], 
                                                      p2017_sf_utm17_focal_obs)  
}


#airborne pollen in 2017 vs pollen production on that day and location in 2017
ggplot(p2017_sf_utm17, aes(x = p_prod_day_1km +1, y = oak_mean)) + geom_jitter(alpha = 0.5, width = 0.1) + theme_few() + 
  geom_smooth(method = "lm", se = FALSE) + #ylab(bquote(oak~pollen~(grains per~m^3))) + 
  xlab("oak pollen released (grains/m2/day)") +
  scale_x_log10()

fit <- lm(p2017_sf_utm17$oak_mean ~ p2017_sf_utm17$p_prod_day_1km)
summary(fit)

fit <- lm(p2017_sf_utm17$oak_mean ~ log10(p2017_sf_utm17$p_prod_day_1km + 1))
summary(fit)
#+ facet_wrap(~name)

#airborne pollen in D 2017 vs airborne pollen at SCS in 2017
ggplot(p17_utm, aes(x = scs_oak, y = oak_mean)) + geom_point(alpha = 0.5) + theme_bw() + 
  #geom_smooth(method = "lm", se = FALSE) + 
  geom_abline(slope = 1, intercept = 0, lty = 2) + 
  coord_cartesian(ylim = c(0,800)) 
fit <- lm(p2017_sf_utm17$oak_mean ~ p2017_sf_utm17$scs_oak)
summary(fit)


## compare SCS estimate with my pollen prod estimate for each obs
#predict airborne pollen for each obs based on my pollen production regression
fit <- lm(p2017_sf_utm17$oak_mean ~ log10(p2017_sf_utm17$p_prod_day_1km + 1))
summary(fit)
fit$coefficients[1]
p2017_sf_utm17$p_prod_day_1km_predQ <- NA
p2017_sf_utm17$p_prod_day_1km_predQ <- fit$coefficients[1] + fit$coefficients[2] * log10(p2017_sf_utm17$p_prod_day_1km + 1)
p2017_sf_utm17 <- arrange(p2017_sf_utm17, oak_mean) %>% mutate(row_n = 1:nrow(p2017_sf_utm17))

ggplot(p2017_sf_utm17, aes(x = row_n, y = oak_mean)) + geom_point() + theme_few()+
  geom_point(aes(x= row_n, y = p_prod_day_1km_predQ), color = "blue") +
  geom_point(aes(x= row_n, y = scs_oak), color = "red")

filter(p2017_sf_utm17, !is.na(oak_mean) & !is.na(scs_oak) & !is.na(p_prod_day_1km_predQ)) %>%
  ggplot(aes(x = row_n, y = abs(oak_mean - scs_oak))) + geom_point(size = 2) +
  geom_point(aes(x= row_n, y = abs(oak_mean - p_prod_day_1km_predQ)), color = "blue") 
