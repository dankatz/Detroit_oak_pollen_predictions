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
library(scales)
library(boot)


## data assembly: spring 2017 data =====================================================================================
# Airborne pollen data was collected spring of 2017, and is described in Katz et al. 2019:
# Effect of intra-urban temperature variation on tree flowering phenology, airborne pollen, and 
# measurement error in epidemiological studies of allergenic pollen
p17 <- read_csv("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pollen_2017.csv") %>% 
    dplyr::select(name, long, lat, date, jday, oak_mean, oak_sd) %>% 
    mutate(dates = ymd(date),
           taxa = "p_Quercus",
           taxa2 = "Quercus",
           taxa3 = "italic(Quercus)")
p17_utm <- st_as_sf(p17, coords = c("long", "lat"), crs = 4326)
p17_utm <- st_transform(p17_utm, 32617)
p17_utm %>% dplyr::select(name) %>% plot(cex = 5, pch = 17)
#write_sf(p17_utm, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pollen_2017_sf.shp")

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
#write_sf(p18_utm, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pollen_2018_sf.shp")

#load in a shapefile that has Detroit's boundary clipped to the middle tile of WV2 imagery (i.e., the area of tree ID project)
d_wv2_boundary <- read_sf("C:/Users/dsk856/Box/MIpostdoc/Detroit spatial data/Detroit boundary/Detroit_wv2_middle_boundary_v2.shp")
d_wv2_boundary_utm <- st_transform(d_wv2_boundary, 32617)

#only include observations from sites within the AOI
p18_utm_inAOI <- st_join(p18_utm, d_wv2_boundary_utm) %>%  #filter out the observations that are from outside of the oak ID area
  filter(!is.na(is_D))

### using linear interpolation to fill in three missing observations
# visualize missing data during the peak of the season (5 sampling days with highest concentrations)
p18_utm_inAOI %>% 
  filter(date > ymd("2018-05-05") & date < ymd("2018-05-25")) %>% 
  #ggplot(aes(x = date, y = site, color = log10(pollen + 1), label = round(pollen))) + geom_text(size = 5) + theme_bw() + scale_color_viridis_c()
  ggplot(aes(x = date, y = pollen, color = site)) + geom_point() + geom_line() + theme_bw() + scale_color_viridis_d()  

# p18_inAOI <- p18_utm_inAOI 
# p18_inAOI$geometry <- NULL
# test <- 
# p18 %>% 
# group_by(site) %>% 
#   filter(date > ymd("2018-05-05") & date < ymd("2018-05-29")) %>% 
#   mutate(dates = paste0(month(date),"_",mday(date))) %>% 
#   dplyr::select(site, dates, pollen) %>% 
#   pivot_wider(names_from = dates, values_from = pollen, names_prefix = "pol_") %>% 
#   mutate(may_17_14_perc = pol_5_17/pol_5_14)


#interpolation 
p18_inAOI_NAs <- tibble(date = c(ymd("2018-05-17"), ymd("2018-05-21"), ymd("2018-05-21")),  #observations I'll be interpolating
                        site = c("I", "T", "D"), 
                        pollen = c(NA,NA,NA))
p18_geo <- dplyr::select(p18_utm_inAOI, site) %>% unique() #getting the geometry for each site #str(p18_geo) #str(p18_inAOI_NAs)
p18_inAOI_NAs <- left_join(p18_inAOI_NAs, p18_geo) #add the geometry to the NA tibble
p18_utm_inAOI_NAs <- bind_rows(p18_utm_inAOI, p18_inAOI_NAs)%>% #add the NA values (with geometry) to the larger dataframe
  arrange(site, date) %>% 
  mutate(pollen_interp = imputeTS::na_interpolation(pollen)) #linear interpolation of 3 missing pollen observations

p18_utm_inAOI_NAs %>% #visualize some data
  filter(date > ymd("2018-05-01") & date < ymd("2018-05-29")) %>% 
  ggplot(aes(x = date, y = site, color = log(pollen + 1), label = round(pollen_interp))) + geom_text(size = 5) + theme_bw() + scale_color_viridis_c()
  #ggplot(aes(x = date, y = pollen_interp, color = site)) + geom_point() + geom_line() + theme_bw() + scale_color_viridis_d()  

p18_peak_season <- 
  p18_utm_inAOI_NAs %>% 
  filter(date > ymd("2018-05-05") & date < ymd("2018-05-22")) %>% 
  group_by(site) %>% 
  summarize(oak_season_mean = mean(pollen_interp))   



### data assembly: tree location and pollen production map ##############################################
# a map of oak trees in the central portion of Detroit is described in Katz et al. 2020:
# Improved Classification of Urban Trees Using a Widespread Multi-Temporal Aerial Image Dataset
# original shapefile: 

#tree_pred <- st_read("C:/Users/dsk856/Box/MIpostdoc/trees/tree_identificaiton/predictions/pred190715.shp")
#p_Quru <- filter(tree_pred, prdctd_ == "Quercus") %>%  dplyr::select(area)

# Pollen production is described in Katz et al. 2020:
# Pollen production for 13 urban North American tree species: allometric equations for tree trunk diameter and crown area

#p_Quru$predpollen <- p_Quru$area * 0.97 + 17.02 #using equation for red oak
#write_sf(p_Quru, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pred_qusp_pol_prod200310.shp")
p_Quru <- read_sf("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/pred_qusp_pol_prod200310.shp")
sum(p_Quru$predpollen) * 1000000 #convert to actual number (was originally in millions)


# #create a blank raster with approximately 100 x 100 m pixel size
# d_rast <- raster(ncol = 145, nrow = 212) 
# extent(d_rast) <- extent(p_Quru)

#create a blank raster with approximately 10 x 10 m pixel size
d_rast_10m <- raster(ncol = 1452, nrow = 2030) #approximately 10 x 10 m pixel size: ncol = 1452, nrow = 2030
extent(d_rast_10m) <- extent(p_Quru)
d_rast <- d_rast_10m

trees_as_points <- p_Quru %>% st_cast("POINT", do_split = FALSE)
p_rast <- rasterize(trees_as_points, d_rast, field = "predpollen", fun = sum, na.rm = TRUE)
pixel_area <- (res(p_rast)[1] * res(p_rast)[2]) #get pixel area in m2
p_rast <- (p_rast/ pixel_area) * 1000000 #convert to pol/m2 (was originally in millions)
p_rast[p_rast < 0] <- 0 #make sure that none of the trees have values below 0
#p_rast[is.na(p_rast)] <- 0 #set cells with no oak tree (but within extent of tree map) to zero
crs(p_rast) <- CRS('+init=EPSG:32617')
# plot(p_rast)
# projection(p_rast)
p_rast_export <- p_rast
writeRaster(p_rast_export, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_10m_210611.tif", 
            format="GTiff", overwrite = TRUE)

#set NA values in that area to 0
# p_rast
# plot(st_geometry(p17_total_season), add = TRUE)
# plot(st_geometry(d_wv2_boundary_utm), add = TRUE) #this was loaded earlier in script
p_rast <- raster::mask(x = p_rast, mask = d_wv2_boundary_utm)

plot(p_rast)

#writeRaster(p_rast, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_210611.tif", format="GTiff", overwrite = TRUE)
#p_rast <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_200302.tif") #old version




### Fig 2: comparing pollen release over space with airborne pollen concentrations ######################
# 2017 data: 
# ultimately, I don't think this is a good comparison, because the samples were not collected simultaneously, 
# and each one is over a short duration
#
# p17_total_season <- p17_utm %>% 
#   group_by(name) %>% 
#   summarize(oak_season_mean = mean(oak_mean))
# 
# plot(p_rast)
# plot(st_geometry(p17_total_season), add = TRUE)
# 
# distance_loop <- c(seq(from = 100, to = 2000, by = 100)) #, seq(from = 1000, to = 10000, by = 1000) #distance_loop <- 400
# 
# #create a table to hold results
# results_df <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
#                          RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))
# results_df_log10 <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
#                                RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))
# 
# for(i in 1:length(distance_loop)){
#   pollen_within_distance_x <- raster::extract(p_rast, p17_total_season, fun = mean, buffer = distance_loop[i], na.rm = TRUE)
#   
#   p17_total_season2 <- p17_total_season %>% 
#     mutate(pollen_within_dist_x = pollen_within_distance_x) %>% 
#     filter(!is.na(pollen_within_distance_x)) #%>% filter(oak_season_mean < 300)
#   
#   panel_a <- p17_total_season2 %>% 
#     ggplot(aes(x = pollen_within_dist_x, y = oak_season_mean)) + geom_point() + theme_bw() + 
#     geom_smooth(method = "lm", se = FALSE, formula = y~log10(x)) + 
#     geom_smooth(method = "lm", se = FALSE, color = "red") + 
#     ggtitle(distance_loop[i]) #+ xlim(0, 20000)
#   
#   panel_b <- p17_total_season2 %>% 
#     ggplot(aes(x = pollen_within_dist_x + 1, y = oak_season_mean)) +  theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
#     geom_point(aes(color = name)) + scale_x_log10(labels = comma) + scale_y_log10()+ annotation_logticks()+ theme(panel.grid.minor = element_blank())+
#     ggtitle(distance_loop[i]) #+ xlim(0, 20000)
#   
#   both_panels <- cowplot::plot_grid(panel_a, panel_b)
#   print(both_panels)
#   
#   fit <- glm(p17_total_season2$oak_season_mean ~ p17_total_season2$pollen_within_dist_x )
#   fit_log10 <- glm(log10(p17_total_season2$oak_season_mean) ~ log10(p17_total_season2$pollen_within_dist_x + 1))
#   
#   results_df$distance[i] <- distance_loop[i]
#   results_df$R2[i] <- caret::R2(fit$fitted.values, p17_total_season2$oak_season_mean )
#   results_df$MAE[i] <- caret::MAE(fit$fitted.values, p17_total_season2$oak_season_mean )
#   results_df$RMSE[i] <- caret::RMSE(fit$fitted.values, p17_total_season2$oak_season_mean )
#   results_df$AIC[i] <- summary(fit)$aic
#   
#   results_df_log10$distance[i] <- distance_loop[i]
#   results_df_log10$R2[i] <- caret::R2(fit_log10$fitted.values, p17_total_season2$oak_season_mean )
#   results_df_log10$MAE[i] <- caret::MAE(fit_log10$fitted.values, p17_total_season2$oak_season_mean )
#   results_df_log10$RMSE[i] <- caret::RMSE(fit_log10$fitted.values, p17_total_season2$oak_season_mean )
#   results_df_log10$AIC[i] <- summary(fit_log10)$aic
# }
# 
# 
# results_df
# results_df_log10





#comparing 2018 data
# plot(p_rast)
# plot(st_geometry(p18_total_season), add = FALSE, pch = 4)
p18_total_season <- p18_peak_season

distance_loop <- c(seq(from = 100, to = 2000, by = 100)) #, seq(from = 1000, to = 10000, by = 1000) #distance_loop <- 400

#create a table to hold results
results_df <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
                    RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))
results_df_log10 <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
                    RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))

for(i in 1:length(distance_loop)){
  pollen_within_distance_x <- raster::extract(p_rast, p18_total_season, fun = mean, buffer = distance_loop[i], 
                                              na.rm = TRUE, na.pad = TRUE)
  
  p18_total_season2 <- p18_total_season %>% 
    mutate(pollen_within_dist_x = pollen_within_distance_x) %>% 
    filter(!is.na(pollen_within_distance_x)) #%>% filter(oak_season_mean < 300)
  
  panel_a <- p18_total_season2 %>% 
    ggplot(aes(x = pollen_within_dist_x, y = oak_season_mean)) + geom_point() + theme_bw() + 
    geom_smooth(method = "lm", se = FALSE, formula = y~log10(x)) + 
    geom_smooth(method = "lm", se = FALSE, color = "red") + 
    ggtitle(distance_loop[i]) #+ xlim(0, 20000)
  
  panel_b <- p18_total_season2 %>% 
    ggplot(aes(x = pollen_within_dist_x + 1, y = oak_season_mean)) +  theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
    geom_point(aes(color = site)) + scale_x_log10(labels = comma) + scale_y_log10()+ annotation_logticks()+ theme(panel.grid.minor = element_blank())+
    ggtitle(distance_loop[i]) #+ xlim(0, 20000)
  
  both_panels <- cowplot::plot_grid(panel_a, panel_b)
  print(both_panels)
  
  fit <- glm(p18_total_season2$oak_season_mean ~ p18_total_season2$pollen_within_dist_x )
  fit_log10 <- glm(log10(p18_total_season2$oak_season_mean) ~ log10(p18_total_season2$pollen_within_dist_x + 1))
  
  results_df$distance[i] <- distance_loop[i]
  results_df$R2[i] <- caret::R2(fit$fitted.values, p18_total_season2$oak_season_mean )
  results_df$MAE[i] <- caret::MAE(fit$fitted.values, p18_total_season2$oak_season_mean )
  results_df$RMSE[i] <- caret::RMSE(fit$fitted.values, p18_total_season2$oak_season_mean )
  results_df$AIC[i] <- summary(fit)$aic
  
  results_df_log10$distance[i] <- distance_loop[i]
  results_df_log10$R2[i] <- caret::R2(fit_log10$fitted.values, log10(p18_total_season2$oak_season_mean + 1) )
  results_df_log10$MAE[i] <- caret::MAE(fit_log10$fitted.values, log10(p18_total_season2$oak_season_mean + 1) )
  results_df_log10$RMSE[i] <- caret::RMSE(fit_log10$fitted.values, log10(p18_total_season2$oak_season_mean + 1) )
  results_df_log10$AIC[i] <- summary(fit_log10)$aic #summary(fit_log10)
}


results_df
results_df_log10 #writeClipboard(as.character(round(results_df_log10, 3)))
# autoplot(fit_log10, which = 1:6, ncol = 3, label.size = 3)
# autoplot(fit, which = 1:6, ncol = 2, label.size = 3, colour = "steelblue") + theme_bw()

## boot strapping to get empirical CI values for the best regression (log x log at x m buffer)
pollen_within_distance_xm <- raster::extract(p_rast, p18_total_season, fun = mean, buffer = 700, na.rm = TRUE)
bootxm <- data.frame(pxm = log10(pollen_within_distance_xm + 1), oak_season_mean = log10(p18_total_season2$oak_season_mean + 1),
                      site = p18_total_season2$site)

# Bootstrap 95% CI for regression coefficients
bs <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(coef(fit))
}
# bootstrapping with 1000 replications
results <- boot(data=bootxm, statistic= bs, R=1000, formula = oak_season_mean ~ pxm)

# # view results
# results
# plot(results, index=1) # intercept
# plot(results, index=2) # slope
# 
# # get 95% confidence intervals
# boot.ci(results, type="bca", index=1) # intercept
# boot.ci(results, type="bca", index=2) # slope

resultsxm_reg <- as.data.frame(results$t) %>% 
  rename(inter = `V1`, slope = `V2`) %>% 
  mutate(x = min(bootxm$pxm),
         y = inter + slope * x,
         xend = max(bootxm$pxm),
         yend = inter + slope * xend,
         id = 1:1000) %>% 
  sample_frac(1) #reduce when doing exploratory figure generation so it renders faster

## figure 2: plotting
bootxm %>% 
  ggplot(aes(x = 10^pxm, y = 10^oak_season_mean)) +  theme_bw() + 
  geom_segment(data = resultsxm_reg, aes(x = 10^x, y = 10^y, xend = 10^xend, yend = 10^yend, group = id), color = "gray50", alpha = 0.01) +
  geom_smooth(method = "lm", se = FALSE) + geom_point(size = 2) + coord_cartesian(ylim = c(100, 3000)) +
  scale_x_log10(labels = comma) + scale_y_log10()+ annotation_logticks()+ theme(panel.grid.minor = element_blank()) +
  xlab(pollen~production~(pollen~grains/m^2)) + ylab(average~airborne~pollen~(pollen~grains/m^3)) 
  
# # Bootstrap 95% CI for R-Squared
# library(boot)
# function to obtain R-Squared from the data
rsq <- function(formula, data, indices) {
  d <- data[indices,] # allows boot to select sample
  fit <- lm(formula, data=d)
  return(summary(fit)$r.square)
}
# bootstrapping with 1000 replications
results <- boot(data=bootxm, statistic=rsq,
                R=1000, formula = oak_season_mean~pxm)

# view results
results
plot(results)

# get 95% confidence interval
boot.ci(results, type="bca")



# saving the raster with the average at the selected distance (700 m) -------------------------------------
# the focal function  doesn't work with a circular mean, using Terra::focal instead
#https://gis.stackexchange.com/questions/358923/how-can-i-get-the-correct-mean-from-focal-with-a-circular-window?noredirect=1&lq=1
library(terra)
focal_distance_selected <- 700 #based on distance selection table (above)

p_rast_AOI <- p_rast
p_rast_AOI[is.na(p_rast_AOI)] <- 0 #turning false NA values to 0
p_rast_AOI <- raster::mask(p_rast_AOI, d_wv2_boundary_utm) #true NA values  #plot(p_rast_AOI)

p_spat_rast <- terra::rast(p_rast_AOI) #convert to terra format #plot(p_rast) #plot(p_spat_rast)
p_spat_rast <- terra::aggregate(p_spat_rast, fact = 2, fun = "mean", na.rm =TRUE) #just to speed up getting to the visualization

#focal weight matrix
fwModel <- focalMat(x = p_spat_rast, d = focal_distance_selected, type = "circle")
fwModel[fwModel > 0] <- 1
fwModel[fwModel == 0] <- NA

p_spat_rast_focal1 <- terra::focal(x= p_spat_rast, w = fwModel, fun = "mean", na.rm = TRUE) #takes 2 min at r = 70 cells
#plot(p_spat_rast_focal1); plot(d_wv2_boundary_utm, add = TRUE, color = NA)
p_spat_rast_focal2 <- raster::raster(p_spat_rast_focal1)
p_spat_rast_focal3 <- raster::mask(p_spat_rast_focal2, d_wv2_boundary_utm)
plot(p_spat_rast_focal3)
writeRaster(p_spat_rast_focal3,  overwrite = TRUE, 
            "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_season_210614d.tif", format="GTiff")

p_spat_rast_focal4 <- log10(p_spat_rast_focal3 + 1)
plot(p_spat_rast_focal4)
writeRaster(p_spat_rast_focal4,  overwrite = TRUE, 
            "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/p_prod_quru_log10_season_210614d.tif", format="GTiff")





# ## testing out the siland package ------------------------------------------------
# library(siland)
# #data(dataSiland)
# #data(landSiland)
# #landSiland$L3 <- 1
# 
# #resF1=Fsiland( obs ~ L3, land=landSiland, data=dataSiland)
# #plotFsiland.sif(resF1)
# 
# p_Quru2 <- p_Quru %>%  mutate(placeholder = 1) #predpollen
# #p_Quru2 <- st_set_crs(p_Quru2, crs(p18_total_season))
# #p_Quru2 <- sample_frac(p_Quru2, 0.1)
# 
# p18_total_season2 <- p18_total_season %>%  mutate(X = sf::st_coordinates(.)[,1],
#                              Y = sf::st_coordinates(.)[,2])
# 
# p18_total_season2$geometry <- NULL
# p18_total_season2 <- p18_total_season2  %>%  filter(oak_season_mean < 300)
# 
# 
# ### buffer approach
# resF1= Bsiland(oak_season_mean ~ placeholder, land= p_Quru2, data= p18_total_season2) #takes ~15 min
# resF1 #best buffer size estimated at: 587.704 m, but stays pretty good from ~ 500 m - 1500 m
# 1 - (resF1$result$deviance/resF1$result$null.deviance)
# summary(resF1)
# str(resF1)
# buffer_selected <- resF1$coefficients[3]
# 
# #a visual comparison of what buffer size has lowest negative log likelihood
# likresB1 = Bsiland.lik(resF1,land= p_Quru2, data=p18_total_season2, varnames=c("placeholder")) #takes a few minutes
# 
# 
# #map of landscape variable
# #plotBsiland.land(x=resF1,land=p_Quru2,data=as.data.frame(p18_total_season2))
# 
# pollen_within_distance_x <- raster::extract(p_rast, p18_total_season, fun = sum, buffer = buffer_selected, na.rm = TRUE)
# 
# p18_total_season <- p18_total_season %>% 
#   mutate(pollen_within_dist_x = pollen_within_distance_x) 
# p18_total_season %>% 
#         ggplot(aes(x = pollen_within_dist_x, y = oak_season_mean)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
#         ggtitle(distance_loop[i]) #+ xlim(0, 20000)
# summary( lm(p18_total_season$oak_season_mean ~ p18_total_season$pollen_within_dist_x ))$r.squared
# summary(glm(p18_total_season$oak_season_mean ~ p18_total_season$pollen_within_dist_x ))$aic
# 
# 
# 
# ### continuous approach
# resF2 = Fsiland(oak_season_mean ~ placeholder, land= p_Quru2, data= p18_total_season2, sif = "gaussian", wd = 10) #+I (placeholder^2) 
# #calculating the Nagelkerke GOF measurement which is about the same as R2 when the GLM is basically the same as the LM
# #https://stats.stackexchange.com/questions/46345/how-to-calculate-goodness-of-fit-in-glm-r/46358
# 1 - (resF2$result$deviance/resF2$result$null.deviance)
# 
# trace(Fsiland, edit=TRUE)
# trace(FsilandMinusLoglik, edit=TRUE)
# #Sys.time()
# resF2
# summary(resF2)
# #str(resF2)
# #resF2$result$fitted.values
# plotFsiland.sif(resF2)
# #?plotFsiland.sif
# Fsiland.quantile(resF2, p = c(0.5, 0.95, 0.99))  #95% of effect is from within 1307 m
# # plotFsiland.land(x=resF2,land=p_Quru2, data=p18_total_season2) #takes a minute, but it's a fairly pretty map of risk
# # test <- plotFsiland.land(x=resF2,land=p_Quru2, data=p18_total_season2, plot = FALSE) #takes a minute, but it's a fairly pretty map of risk
# # ggplot(data=test, aes_string(x="X", y="Y", fill = "V")) + geom_raster(interpolate = F)+
# #   scale_fill_gradient2(low="#0000FF",mid="white",high="#CC0000",midpoint=0)+
# #   coord_fixed()+
# #   theme_classic() + theme(axis.title=element_blank(),legend.title=element_blank(),legend.position="bottom") 
# # 
# # 
# # test2 <- st_as_sf(test, coords = c("X", "Y"), crs = 32617)
# # #test2 <- rasterFromXYZ(test)
# # test3 <- raster(test2) #rasterize(trees_as_points, d_rast, field = "predpollen", fun = sum, na.rm = TRUE)
# # test4 <- rasterize(test3, test2)
# # 
# # 
# # example_points <- as(test2, "Spatial")
# # empty_raster <- raster(test2, ncol = (max(test$X) - min(test$X))/200, nrow = (max(test$Y) - min(test$Y))/200)
# # # Generate empty raster layer and rasterize points
# # example_raster <- #raster(test2) %>%
# #   rasterize(example_points, empty_raster)
# # example_raster2 <- example_raster$V
# # example_raster2[100:1100]
# # plot(example_raster2)
# # example_raster3 <- raster::mask(x = example_raster2, mask = d_wv2_boundary_utm)
# # plot(example_raster3)
# # 
# # writeRaster(example_raster3, "C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/qusp_Fsiland_210528.tif")
# 
# #str(resF2)
# p18_total_season2 %>% 
#   mutate(fitted.valuess = resF2$result$fitted.values,
#          resid = resF2$result$residuals) %>% 
#   ggplot(aes(x = oak_season_mean, y = fitted.valuess, color = resid)) + geom_point() + geom_abline(slope = 1, intercept = 0, lty = 2) + theme_bw() 
#   #xlim(0,600) + ylim(0, 600)
# 
# p18_total_season2 %>% 
#   mutate(landcontri = resF2$landcontri,
#          fitted.valuess = resF2$result$fitted.values) %>% 
#   ggplot(aes(x = landcontri , y = oak_season_mean, color = fitted.valuess)) + geom_point() + geom_abline(slope = 1, intercept = 0, lty = 2) + 
#   theme_bw() + geom_smooth(method = "lm", se = FALSE)
#   
# 
# summary(lm(p18_total_season2$oak_season_mean ~ resF2$landcontri))
# 
# fitted.Fsiland(resF2)
# 
# resF2$landcontri
# resF2$result$model$placeholder
# 
# resF2$result$linear.predictors
# resF2$result$fitted.values
# 
# resF2$result$effects
# 
# plot(resF2$landcontri, resF2$result$fitted.values)




### Fig 3: create pollen release surfaces for each day  #######################################
# setwd("C:/Users/dsk856/Box/MIpostdoc/Detroit spatial data/Detroit boundary/Detroit_boundary_shapefile")
# d_bound <- rgdal::readOGR(".", layer="POLYGON")
# d_bound_utm <- st_as_sf(d_bound)
# d_bound_utm <- st_transform(d_bound_utm, crs(d_peak))
# d_bound_utm <- sf:::as_Spatial(d_bound_utm)

#total pollen production in AOI as a 10m raster
p_rast_AOI

#day of peak pollen release in D
#this data is described in the Sci Tot Env paper and comes from PHENO_ANALYSES_FIGS_V4_181025.R

#create rasters of percent active flowers for oaks on each day in 2017
#this raster was projected to UTM 17N and cell size was adjusted to match the pollen production surface in GIS
d_peak <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/peakflowering_d_UTM17c.tif") 
d_peak <- raster::resample(d_peak, p_rast) #plot(p_rast)
d_peak <- crop(d_peak, p_rast) #plot(d_peak)


#Add a lag to the flowering intensity model to account for differences between street trees and non-street trees
#note: another potential option would be to use the Landsat LST, see: https://www.mdpi.com/2072-4292/12/9/1471/htm
d_peak <- d_peak + 3 # 3 > 2 > 4> 5 > 0 
#d_peak[d_peak[] == -996] <- NA


#FOR ANIMATION SWAP THIS IN- IT GETS RID OF WHITE PIXELS: d_peak <- crop(peakflow_oak_d, d); plot(d_peak)
#d_peak_big <- aggregate(d_peak, 10, fun = mean, na.rm = TRUE) ##plot(peakflow_oak_d_water999)

#I created an empirical look-up table of the average percent of mature flowers on days away from peak
#load in lookup table
pflow_lookup_summary <- read.csv("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/pflow_lookup_summary.csv") 
#apply this function to create a raster for each day with the percent flowering
pflow_lookup_summary$days_from_site_mean_low <-  pflow_lookup_summary$days_from_site_mean - 0.5
pflow_lookup_summary$days_from_site_mean_hi <-  pflow_lookup_summary$days_from_site_mean + 0.5
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
flower_stack <- flower_stack * 0.01 #Changing percentage to proportion

#pollen release on each day
p_per_day_stack <- p_rast * flower_stack #plot(p_per_day_stack[[20]])



## Fig 3: assessing different distance buffers -------------------------------------------------------------
p_per_day_stack2 <- p_per_day_stack
p_per_day_stack2[is.na(p_per_day_stack2)] <- 0 #set na values to 0 (they weren't really NA) #takes ~3 min
p_per_day_stack2 <- mask(p_per_day_stack2,  d_wv2_boundary_utm) #set all values outside of AOI to NA


p17_utm2 <- st_join(p17_utm, d_wv2_boundary_utm) %>% 
  filter(is_D == "d") #%>% filter(oak_season_mean < 300)

distance_loop <- c(seq(from = 100, to = 2000, by = 100)) #, seq(from = 1000, to = 10000, by = 1000) #distance_loop <- c(400, 600)

#create a table to hold results
results_df <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
                         RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))
results_df_log10 <- data.frame(distance = rep(NA, length(distance_loop)), R2 = rep(NA, length(distance_loop)), MAE = rep(NA, length(distance_loop)),
                               RMSE = rep(NA, length(distance_loop)), AIC = rep(NA, length(distance_loop)))

for(i in 1:length(distance_loop)){
   focal_distance <- distance_loop[i] #focal_distance <- 100
  
  #extract for the right day (raster layer) for each sample in 2017
  for(j in 1:nrow(p17_utm2)){
    date_obs <- p17_utm2$jday[j]
    p17_utm_focal_obs <- p17_utm2[j,] #    p17_utm_focal_obs <- p17_utm2[1,]
    
    #pol_prod_day_1km_m has 36 layers for sequential days. layer 1 == jday 110, layer 2 == jday 111, ...
    relevant_layer <- p17_utm_focal_obs$jday - 109  #select the right layer #relevant_layer <- 115 - 109
    p17_utm2$p_prod_day_xm[j] <- raster::extract(p_per_day_stack2[[relevant_layer]], p17_utm_focal_obs,
                                                 fun = mean, buffer = focal_distance, na.rm = TRUE)  
  } 
 
  panel_a <- p17_utm2 %>% 
    ggplot(aes(x = p_prod_day_xm, y = oak_mean)) + geom_point() + theme_bw() + 
    #geom_smooth(method = "lm", se = FALSE, formula = y~log10(x)) + 
    geom_smooth(method = "lm", se = FALSE, color = "red") + 
    ggtitle(focal_distance) #+ xlim(0, 20000)
  
  panel_b <- p17_utm2 %>% 
    ggplot(aes(x = p_prod_day_xm + 1, y = oak_mean + 1)) +  theme_bw() + geom_smooth(method = "lm", se = FALSE) + 
    geom_point(aes(color = name)) + scale_x_log10(labels = comma) + scale_y_log10()+ annotation_logticks()+ theme(panel.grid.minor = element_blank())+
    ggtitle(focal_distance) #+ xlim(0, 20000)
  
  both_panels <- cowplot::plot_grid(panel_a, panel_b)
  print(both_panels)
  
  fit <- glm(p17_utm2$oak_mean ~ p17_utm2$p_prod_day_xm )
  fit_log10 <- glm(log10(p17_utm2$oak_mean + 1) ~ log10(p17_utm2$p_prod_day_xm + 1)) #summary(fit_log10)
  
  results_df$distance[i] <- focal_distance
  results_df$R2[i] <- caret::R2(fit$fitted.values, p17_utm2$oak_mean )
  results_df$MAE[i] <- caret::MAE(fit$fitted.values, p17_utm2$oak_mean )
  results_df$RMSE[i] <- caret::RMSE(fit$fitted.values, p17_utm2$oak_mean )
  results_df$AIC[i] <- summary(fit)$aic

  results_df_log10$distance[i] <- focal_distance
  results_df_log10$R2[i] <- caret::R2(fit_log10$fitted.values, log10(p17_utm2$oak_mean + 1) ) #summary(fit_log10)
  results_df_log10$MAE[i] <- caret::MAE(fit_log10$fitted.values, log10(p17_utm2$oak_mean + 1) )
  results_df_log10$RMSE[i] <- caret::RMSE(fit_log10$fitted.values, log10(p17_utm2$oak_mean + 1) )
  results_df_log10$AIC[i] <- summary(fit_log10)$aic
}

results_df
round(results_df_log10, 3)



## creating fig 3A using the distance selected above ---------------------------------------------------
focal_distance_selected <- 1000

p_per_day_stack2 <- p_per_day_stack
p_per_day_stack2[is.na(p_per_day_stack2)] <- 0 #set na values to 0 (they weren't really NA) #takes ~3 min
p_per_day_stack2 <- mask(p_per_day_stack2,  d_wv2_boundary_utm) #set all values outside of AOI to NA

p17_utm3 <- st_join(p17_utm, d_wv2_boundary_utm) %>% 
  filter(is_D == "d") %>% 
  mutate(p_prod_day_xm = NA)  #%>% filter(oak_season_mean < 300)

#extract for the right day (raster layer) for each sample in 2017
for(j in 1:nrow(p17_utm3)){
  date_obs <- p17_utm3$jday[j]
  p17_utm_focal_obs <- p17_utm3[j,] #    p17_utm_focal_obs <- p17_utm3[1,]
  
  #pol_prod_day_1km_m has 36 layers for sequential days. layer 1 == jday 110, layer 2 == jday 111, ...
  relevant_layer <- p17_utm_focal_obs$jday - 109  #select the right layer #relevant_layer <- 115 - 109
  p17_utm3$p_prod_day_xm[j] <- raster::extract(p_per_day_stack2[[relevant_layer]], p17_utm_focal_obs,
                                               fun = mean, buffer = focal_distance_selected, na.rm = TRUE)  
} 

#airborne pollen in 2017 vs pollen production on that day and location in 2017
ggplot(p17_utm3, aes(x = p_prod_day_xm +1, y = oak_mean + 1)) + geom_point(alpha = 0.5) + theme_bw() + 
  geom_smooth(method = "lm", se = FALSE) + ylab(airborne~oak~pollen~(grains~per~m^3)) + 
  xlab(oak~pollen~released~daily~(grains/m^2)) +
  scale_x_log10(labels = comma) + scale_y_log10() + annotation_logticks() +
  theme(panel.grid.minor = element_blank())

fit <- lm(p17_utm3$oak_mean  ~ p17_utm3$p_prod_day_xm)
summary(fit)

fit <- lm(log10(p17_utm3$oak_mean + 1) ~ log10(p17_utm3$p_prod_day_xm + 1))
summary(fit)


## compare SCS estimate with my pollen prod estimate for each obs -----------------------------
##add in SCS data that are available for 2017
scs <- read.csv("C:/Users/dsk856/Box/MIpostdoc/LakeshoreENT_pollen_counts_2009_2017_compiled.csv") %>% 
  mutate(date = mdy(Date),
       year = year(date),
       jday = yday(date)) %>% 
  filter(year == 2017) %>%
  dplyr::select(date, Quercus, jday) %>%
  rename(scs_oak = Quercus)

p17_utm4 <- left_join(p17_utm3, scs_c)

#predict airborne pollen for each obs based on my pollen production regression
fit <- lm(p17_utm4$oak_mean ~ p17_utm4$scs_oak)
summary(fit)
fit <- lm(log10(p17_utm4$oak_mean + 1) ~ log10(p17_utm4$scs_oak + 1))
summary(fit)

#direct comparison
ggplot(p17_utm4, aes(x = scs_oak + 1, y = oak_mean + 1)) + geom_point() + theme_bw() + geom_smooth(method = "lm", se = FALSE) +
  scale_x_log10() + scale_y_log10()

#time series
ggplot(p17_utm3, aes(x = date, y = oak_mean + 1, group = name)) + geom_point() + theme_bw() + geom_line() +
  geom_line(data = scs, aes(x= date, y = scs_oak + 1, group = "test"), color = "red") + scale_y_log10()




# ## compare Sylvania/Toledo estimate with my pollen prod estimate for each obs ?
# tol <- read.csv("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/NAB_sylvania_Toledo_data.csv") %>% 
#   mutate(date = mdy(Date),
#          year = year(date),
#          jday = yday(date),
#          name = "toledo") %>% 
#   #filter(year == 2017) %>%
#   dplyr::select(date, Quercus, jday) %>%
#   rename(tol_oak = Quercus)
# 
# #2017
# p17_utm3 <- left_join(p17_utm, tol)
# fit <- lm(p17_utm3$oak_mean ~ p17_utm3$tol_oak)
# summary(fit)
# ggplot(p17_utm3, aes(x = tol_oak, y = oak_mean)) + geom_point() + theme_bw() + geom_smooth(method = "lm")
# 
# ggplot(p17_utm3, aes(x = date, y = oak_mean + 1, group = name)) + geom_point() + theme_bw() + geom_line() +
#   geom_line(data = tol, aes(x= date, y = tol_oak + 1, group = "test"), color = "red") + scale_y_log10()
# 
# #2018
# p18_utm3 <- left_join(p18_utm, tol) %>% 
#   mutate(year = year(date)) %>% 
#   filter(year == 2018)
# 
# #predict airborne pollen for each obs based on my pollen production regression
# fit <- lm(p18_utm3$pollen ~ p18_utm3$tol_oak)
# summary(fit)
# ggplot(p18_utm3, aes(x = tol_oak, y = pollen)) + geom_point() + theme_bw() + geom_smooth(method = "lm") +
#   scale_x_log10() + scale_y_log10() + geom_abline(slope =1, intercept = 0, lty = 2)
# 
# ggplot(p18_utm3, aes(x = date, y = pollen + 1, group = site)) + geom_point() + theme_bw() + geom_line() +
#   geom_line(data = tol, aes(x= date, y = tol_oak + 1, group = "test"), color = "red", lwd = 2) + scale_y_log10() +
#   scale_x_date(lim = c(ymd("2018-04-05"), ymd("2018-06-01")))


### old stuff ##########################################################################

### data assembly: phenology model data #################################################################
# # oak phenology data was collected spring of 2017, and is described in Katz et al. 2019:
# # Effect of intra-urban temperature variation on tree flowering phenology, airborne pollen, and 
# # measurement error in epidemiological studies of allergenic pollen
# # see: PHENO_ANALYSES_FIGS_V4_181025.R
# 
# #create rasters of percent active flowers for oaks in 2017
# #this raster was projected to UTM 17N and cell size was adjusted to match the pollen production surface in GIS
# d_peak <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/peakflowering_d_UTM17c.tif") 
# d_peak <- raster::resample(d_peak, p_rast)
# d_peak <- crop(d_peak, p_rast)
# #plot(d_peak)
# 
# # peakflow_oak_d <- raster("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/peakflowering_d_UTM17.tif") #plot(peakflow_oak_d)
# # d <- extent(-83.3, -82.86, 42.30, 42.6)
# # d_peak <- crop(peakflow_oak_d, d)
# # plot(peakflow_oak_d)
# 
# #Add a lag to the flowering intensity model to account for differences between street trees and non-street trees
# d_peak <- d_peak + 3 
# d_peak[d_peak[] == -996] <- NA
# 
# ##FOR ANIMATION SWAP THIS IN- IT GETS RID OF WHITE PIXELS: d_peak <- crop(peakflow_oak_d, d); plot(d_peak)
# #d_peak_big <- aggregate(d_peak, 10, fun = mean, na.rm = TRUE) #plot(d_peak_big) #plot(peakflow_oak_d_water999)
# #vmap$peak_oak <- raster::extract(x = peakflow_oak_d, y = cbind(vmap$long, vmap$lat), buffer = 500, fun = mean, na.rm = TRUE)
# 
# #I created an empirical look-up table of the average percent of mature flowers on days away from peak
# #load in lookup table
# pflow_lookup_summary <- read_csv("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/pflow_lookup_summary.csv") %>% 
#   mutate(days_from_site_mean_low = days_from_site_mean - 0.5,
#          days_from_site_mean_hi = days_from_site_mean + 0.5)
# #apply this function to create a raster for each day with the percent flowering
# pflow_lookup_summary2 <- cbind(pflow_lookup_summary$days_from_site_mean_low, pflow_lookup_summary$days_from_site_mean_hi,
#                                pflow_lookup_summary$mean_flow)
# 
# #creating a stack that has percent flowering on each day
# flower_stack <- stack()
# for(i in 110:145){
#   focal_day <- i
#   test2 <- focal_day - d_peak 
#   test3 <- reclassify(test2, pflow_lookup_summary2)
#   names(test3) <- paste("flow_", focal_day, sep = "")
#   plot(test3, main = focal_day)
#   flower_stack <- stack(flower_stack, test3)
# }
#
# ### predictions of pollen release over space and time ###################################################
# plot(p_rast)
# flower_stack <- flower_stack * 0.01 #Changing percentage to proportion
# 
# # for(i in 1:20){plot(flower_stack[[i]])}
# # plot(flower_stack[[20]])
# p_per_day_stack <- p_rast * flower_stack
# 
# 
# 
# # plot(flower_stack[[1]])
# # plot(p_per_day_stack[[2]])
# # for(i in 1:36){plot(p_per_day_stack[[i]])}
# #p_per_day_stack_test <- p_per_day_stack
# #p_per_day_stack[p_per_day_stack[] == 0] <- NA #trying this so I can get white cells when no pollen prod
# 
# # library(rgdal)
# # setwd("C:/Users/dsk856/Box/MIpostdoc/trees/Phenology and daily variation in pollen release/detroit cartography/detroit_water")
# # d_water <-readOGR(".", layer="detroit_water1")
# # d_water_utm <- st_as_sf(d_water)
# # d_water_utm <- st_transform(d_water_utm, crs(d_peak))
# # d_water_utm <-sf:::as_Spatial(d_water_utm)
# # 
# # d_bound <- read_sf("C:/Users/dsk856/Box/MIpostdoc/Detroit spatial data/Detroit boundary/Detroit_wv2_middle_boundary_v2.shp")
# # d_bound_utm <- st_transform(d_bound, crs(d_peak))
# 
# 
# 
# p_per_day_stack[is.na(p_per_day_stack)] <- 0 #set na values to 0 (they weren't really NA)
# p_per_day_stack <- mask(p_per_day_stack,  d_wv2_boundary_utm) #set all values outside of D to NA
# # plot(p_per_day_stack[[10]])
# # plot(mask(p_per_day_stack[[1]], d_bound_utm))
# 
# #library('viridis') #throws an error, but still works fine
# # colr <- colorRampPalette(brewer.pal(10, 'RdYlGn'))
# # myColorkey <- list( at= c(0, 100, 1000, 10000, 100000), ## where the colors change #max(p_per_day_stack)
# #                    labels=c("0", "100", "1,000", "10,000", "100,000"), ## where to print labels
# #                    title = "pollen grains/day")
# # 
# # levelplot(x = p_per_day_stack, layers = 36, margin=FALSE, main = "Julian day 125", #cuts = 10, #at = seq(0, 100, by = 20),
# #           colorkey = myColorkey, col.regions = viridis(25), background = "black") + #par.settings=list(panel.background=list(col="black")) +
# #   layer(sp.polygons(d_bound_utm, lwd = 1)) + layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue")) 
# #   
# 
# # 
# # my.at=c(0,10, 50,100,500, 1000,5000, 10000,50000, 100000, 500000)
# # my.brks=seq(0, 100, by=10)
# # 
# # myColorkey <- list(at=my.brks, labels=list(at=my.brks, labels=c(0,10, 50,100,500, "1,000","5,000", "10,000","50,000", "100,000", "500,000")), space="bottom")
# # levelplot(p_per_day_stack, layers = 13, margin=FALSE, main = "Julian day 125", at= my.at, 
# #           colorkey=myColorkey, col.regions = viridis(25, option = "A", direction = -1), 
# #           background = "black", ylab = "", xlab = "",
# #           scales=list(x=list(draw=FALSE), y=list(draw=FALSE)))
# # 
# # 
# # #ANIMATE
# # #saving a png file for each graph
# # setwd("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/animation")
# # p_per_day_stack
# # #month and day as.Date((109 + 1), origin = "2016-12-31")
# # for(i in 1:36){ #1:36
# #   date_title <- paste("pollen release on:", as.Date((109 + i), origin = "2016-12-31")) #convert Julian to date
# #   png(filename=paste("D_qf_",109 + i,".png", sep = ""), width = 600, height = 800, units = "px")
# #   print(levelplot(x = p_per_day_stack, layers = i, margin=FALSE, main = date_title, cuts = 10, at = my.at,
# #                   colorkey = myColorkey, 
# #                   col.regions = viridis(25, option = "A", direction = -1),
# #                   ylab = "", xlab = "", scales=list(x=list(draw=FALSE), y=list(draw=FALSE))) + #
# #           layer(sp.polygons(d_bound_utm, lwd = 1))) +
# #     layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue"))
# #   
# #   dev.off()
# # }
# # 
# # #creating an animation
# # 
# # list.files(path = ".", pattern = "*.png", full.names = T) %>%
# #   purrr::map(image_read) %>% # reads each path file
# #   image_join() %>% # joins image
# #   image_animate(fps=1) %>% # animates, can opt for number of loops
# #   image_write("oak_pollen_prod_in_Detroit_v1.gif") # write to current dir
# 
# 
# # ## create a version where each pixel is the average of all pixels within 1km of it ------------------------------
# # fwModel <- focalWeight(p_per_day_stack, 1000, type='circle')
# # fwModel[fwModel>0] <- 1
# # # km_mean_fun <- function(x){focal(x, w=fwModel, fun=mean, na.rm = TRUE)}
# # # plot(km_mean_fun(p_per_day_stack[[2]]))
# # 
# # #I'm having a hard time applying this function to the whole raster brick, so I'm just doing it with a loop
# # pol_prod_day_1km <- stack() # create an empty stack
# # for (i in 1:nlayers(p_per_day_stack)){
# #   mean_pollen_prod1km <- focal(p_per_day_stack[[i]], w=fwModel, fun=mean, na.rm = TRUE)
# #   pol_prod_day_1km <- stack(pol_prod_day_1km , mean_pollen_prod1km )
# # }#plot(x[[10]])
# # 
# # pol_prod_day_1km_m <- pol_prod_day_1km
# # pol_prod_day_1km_m[is.na(pol_prod_day_1km_m)] <- 0 #set na values to 0 (they weren't really NA)
# # pol_prod_day_1km_m <- mask(pol_prod_day_1km_m,  d_bound_utm) #set all values outside of D to NA
# # 
# # 
# # #month and day as.Date((109 + 1), origin = "2016-12-31")
# # setwd("C:/Users/dsk856/Box/MIpostdoc/trees/airborne_pollen/animation_prod__per_day_within1km")
# # 
# # my.at=c(0, 100, 500, 1000, 2500, 5000, 7500, 10000, 15000, 20000, 30000,40000, 60000)
# # my.brks=seq(0, 100, by=10)
# # myColorkey <- list(at=my.brks, 
# #                    labels=list(at=my.brks,
# #                                labels= my.at, #c(0,10, 50,100,500, "1,000","5,000", "10,000","50,000", "100,000", "500,000")
# #                                #space="right", 
# #                                #title = "pollen release within 1 km (grains/m2)"
# #                    ))
# # 
# # 
# # for(i in 1:36){ #1:36
# #   date_title <- paste("pollen release on:", as.Date((109 + i), origin = "2016-12-31")) #convert Julian to date
# #   png(filename=paste("D_qf_",109 + i,".png", sep = ""), width = 600, height = 800, units = "px")
# #   print(levelplot(x = pol_prod_day_1km_m, layers = i, margin=FALSE, xlim=c(320420,333200), ylim=c(4680362,4700550),
# #                   main = date_title, cuts = 500, at = my.at, 
# #                   colorkey = myColorkey,  
# #                   col.regions = viridis(500, option = "A", direction = -1),
# #                   ylab = "", xlab = "", scales=list(x=list(draw=FALSE), y=list(draw=FALSE))) + #
# #           layer(sp.polygons(d_bound_utm, lwd = 1)) 
# #         #layer(sp.polygons(d_water_utm, lwd = 2, fill = "skyblue"))
# #   )
# #   dev.off()
# # }
# # 
# # 
# # list.files(path = ".", pattern = "*.png", full.names = T) %>%
# #   purrr::map(image_read) %>% # reads each path file
# #   image_join() %>% # joins image
# #   image_animate(fps=1) %>% # animates, can opt for number of loops
# #   image_write("oak_pollen_release_within_1km_v1.gif") # write to current dir
