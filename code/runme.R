## ---- include=FALSE-------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message = FALSE, warning = FALSE)


## ----setup----------------------------------------------------------
library(RandomFields)
library(INLA)
library(inlabru)
library(sp)
library(spatstat)
library(sf)
library(scico)
library(tidyverse)
library(patchwork)
library(maptools)

theme_map =  
  theme_light() + 
    theme(axis.ticks.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.y = element_blank(),
          axis.title.x=element_blank(),
          axis.title.y=element_blank()) 

# function that computes the value of the covariate in any point of the area of interest.
# this is required by inlabru

f.fill = function(x,y, covar) 
{
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string = CRS(my_crs))
  v <- over(spp, covar) 
  v[is.na(v)] <- 0 # NAs are a problem! Remove them
  return(v[,1])
}




## ----read_files-----------------------------------------------------

# read shapefiles -------------------------------------------

poly = read_sf("shape_interest/shape_interest.shp")
poly_sp = as(poly,"Spatial")

ferry_sp = read_sf("ferry_lines/ferry_lines.shp") %>%
  st_set_crs(proj4string(poly_sp))
ferry_sp = as(ferry_sp, "Spatial")

my_crs = proj4string(poly_sp)

depth =  read.table("depth.csv") 
depth_SPDF = SpatialPixelsDataFrame(points = cbind(depth$x,depth$y), 
                                    data = data.frame(depth = depth$z),
                                    proj4string = CRS(my_crs)) 


int_social = read.table("intensity.dat", header = T)
int_social_sp = SpatialPixelsDataFrame(cbind(x = int_social$x, y = int_social$y), 
                                  data = data.frame(cov = int_social$z),
                                  proj4string = CRS(my_crs))


win <- as.owin(poly)


## ----mesh, fig.width = 4, fig.height=  4----------------------------
# CREATE THE MESH
boundary = as(st_simplify(poly, dTolerance = 1.3),"Spatial") 
max.edge = 10
bound.outer = 50
mesh = inla.mesh.2d(boundary = boundary,
                          max.edge = c(1,4)*max.edge,
                          # - use 5 times max.edge in the outer extension/offset/boundary
                          crs = CRS(my_crs),
                          offset = c(max.edge, bound.outer))

A = inla.spde.make.A(mesh = mesh, loc = as.matrix(depth[,c(1,2)]))


## ----fig_mesh, fig.width = 4, fig.height=4, echo = FALSE------------

ggplot() + geom_sf(data = poly) +
  gg(mesh) + theme_map


## ----simulate-------------------------------------------------------
# simulate point process -------------------------------------------------
RFoptions(seed=67356453)


# create the intensity function
# linear part
beta0 = -4
beta_depth = -1
lambda = beta0 + beta_depth * depth$z
im = as.im(data.frame(x = depth$x, y = depth$y, z = lambda), W = win)


# parameters for the Gaussian field
sigma2x <- .5^2
range <- 100

# Simulate the sightings process
data <- rLGCP(model="matern", 
              mu=im, 
              var=sigma2x, 
              scale=range/sqrt(8), 
              nu=1, 
              win = win,
              saveLambda = TRUE)
data.sp <- SpatialPoints(cbind(data$x,data$y),
                         proj4string = CRS(my_crs)) 
# Save the intensity surface for later reference
Lambda = raster::raster(attributes(data)$Lambda)



## ----int:plot, fig.width = 7, fig.height= 7, echo = FALSE, fig.cap="Left: observed point process, in the background is the density of the LGCP. Right: density for the SM observation process", echo = FALSE----
r1 = raster::mask(Lambda, poly)
r2 = raster::crop(r1, poly)
r3 <- raster::rasterToPoints(r2, spatial = TRUE, proj4string = CRS(my_crs)) 


p1 = ggplot() + gg(r3, aes(color = layer)) + scale_color_scico() +
 theme_map + geom_sf(data = poly, alpha = 0) + theme(legend.position = "none") + 
  geom_point(data = data.frame(x = data$x, y = data$y), aes(x,y))

p2 = ggplot() + gg(int_social_sp) + theme_map + scale_fill_scico() + theme(legend.position = "none") + 
  geom_sf(data = poly, alpha = 0)

p1 + p2


## ----detection_param, echo = FALSE----------------------------------
scale_detect1 = 0.5
location_detect1 = 0
sig_detect2 = 2


## ----thinning-------------------------------------------------------
# function to compute the distance between a point in space and the nearest ferry track
surf2 = function(x,y)
{
  
  d = SpatialPoints(cbind(x,y),
                    proj4string = CRS(my_crs))
  d1 = rgeos::gDistance(d, ferry_sp, byid=TRUE)
  z = apply(d1,2,min)
  return(z)
}

# thinning the data points

# complete dataset
all_data = data.frame(x = data$x,
                      y = data$y)
  

# dataset from SM observation process
data1 = all_data %>% 
  mutate(surf = f.fill(x,y,int_social_sp)) %>%
  mutate(pdetect = pnorm(surf, mean = location_detect1, sd = scale_detect1)) %>%
  mutate(u = runif(data$n)) %>%
  mutate(detect = u<pdetect) %>%
  dplyr::filter(detect)

# dataset from ferry observation process
data2 = all_data %>% 
  mutate(surf = surf2(x,y)) %>%
  mutate(pdetect = exp(-0.5*(surf/ sig_detect2)^2),
         u = runif(data$n)) %>%
  mutate(detect = u<pdetect) %>%
  dplyr::filter(detect)


## ---- echo = FALSE, fig.width = 6, fig.height=  6-------------------
p1 = ggplot() + geom_point(data = all_data, aes(x,y), alpha = 0.2) +
  geom_point(data = data1, aes(x,y)) +
  geom_sf(data = poly, alpha = 0) +
  ggtitle("SM Data ") + theme_map
 
p2 = ggplot() + geom_point(data = all_data, aes(x,y), alpha = 0.2) +
  gg(ferry_sp) + 
  geom_point(data = data2, aes(x,y))+
    geom_sf(data = poly, alpha = 0) +
  ggtitle("Ferry Data ")+ theme_map
p1 + p2


## -------------------------------------------------------------------
# create the datasets as SP and all maps as SPDF --------------------------
data1_sp = SpatialPointsDataFrame(cbind(x = data1$x, y = data1$y), 
                                           data = data.frame(d = rep(1,length(data1$x))),
                                           proj4string = CRS(my_crs)) 

data2_sp = SpatialPointsDataFrame(cbind(x = data2$x, y = data2$y), 
                                  data = data.frame(d = rep(1,length(data2$x))),
                                  proj4string =CRS(my_crs))  

# we need to add the distance from the ferry track for the ferry data
d1 = rgeos::gDistance(data2_sp, ferry_sp, byid=TRUE)
d1 = apply(d1,2,min)
data2_sp$dist_obs = d1


## ----prior0, echo = F-----------------------------------------------
rate_detect1 = 1
rate_detect2 = 0.5


## ----priors---------------------------------------------------------
# Create the SPDE model 
matern <-  inla.spde2.pcmatern(mesh,
                                   prior.sigma = c(.7, 0.01),
                                   prior.range = c(100, 0.5))


# Transformation function the the \xi paramters
qexppnorm <- function(x, rate) {
  qexp(pnorm(x, lower.tail = FALSE, log.p = TRUE),
       rate = rate, 
       log.p = TRUE, 
       lower.tail = FALSE)
}

# Define the two detection functions --------------------------------------
log_detect_social = function(x,y, t1, t2)
{
  dens = f.fill(x,y,int_social_sp)
  pnorm(dens/qexppnorm(t1, rate = rate_detect1) +
          t2,
        log.p = TRUE)
}

log_detect_ferry = function(dist, sig)
{ 
  -0.5*(dist/ sig)^2
}


## ----run-model, message=FALSE---------------------------------------
# Define the model components


cmp = ~ 
  SPDE(main = coordinates, model = matern) +
  depth(depth_SPDF) +
  Intercept(1)+
  sig_detect2(main = 1, model = "linear",
             mean.linear = 0,
             prec.linear = 1) +
  scale_detect1(1, 
                mean.linear = 0,
                prec.linear = 1)+
  location_detect1(1,  
                   mean.linear = 0,
                   prec.linear = 10) 

# Define the formula for the two observed processes
formula_data1 = coordinates   ~ 
  Intercept +
  depth  + 
  SPDE +
  log_detect_social(x,y,location_detect1, scale_detect1)
  
formula_data2 = coordinates + dist_obs   ~ 
  Intercept +
  depth  + 
  SPDE +
  log_detect_ferry(dist_obs,qexppnorm(sig_detect2,rate =rate_detect2))

# Define the likelihood for the two observed processes
lik_data1 <- like("cp",
                        formula = formula_data1,
                        samplers = poly_sp, 
                        domain = list(coordinates = mesh),
                        data = data1_sp)

lik_data2 <- like("cp",
                  formula = formula_data2,
                  samplers = ferry_sp,
                  domain = list(coordinates = mesh,
                                dist_obs = inla.mesh.1d(seq(0, 7, length.out = 5),
                                                        degree = 1)),
                  data = data2_sp)

# Finally set things together and run the model
bru_options_set(bru_verbose = 1)



fit <- bru(components = cmp,  
         lik_data1,
         lik_data2,
         options = list(verbose = F,
                        num.threads = "1:1"))






## ----res-linear, fig.width = 4, fig.height= 4, echo = FALSE---------

data.frame(rbind(fit$marginals.fixed$Intercept,
      fit$marginals.fixed$depth ,
      fit$marginals.hyperpar$`Range for SPDE`,
      fit$marginals.hyperpar$`Stdev for SPDE`)) %>%
  mutate(param = rep(c("Intercept","beta","range for SPDE", "sd for SPDE"), each = 75)) %>%
  dplyr::filter(x<1000) %>%
  dplyr::filter((x<5 &  param=="sd for SPDE") | param!="sd for SPDE" ) %>%
  ggplot() + geom_line(aes(x,y,group = param)) + 
  facet_wrap(.~param, scales = "free") +
  geom_vline(data = data.frame(x = c(beta0, beta_depth, range, sqrt(sigma2x)), 
                                          param = c("Intercept","beta","range for SPDE", "sd for SPDE")),
                        aes(xintercept = x))



## ----res-social,fig.width = 4, fig.height= 4------------------------
# detection function for ferry data
scale_social = function(...)
{
  qexppnorm(sig_detect2, rate = 1)
}

# we sample from the model
samples = inla.posterior.sample(500, fit)

scale1 = inla.posterior.sample.eval(fun = scale_social,samples = samples)

df = data.frame(scale = as.numeric(scale1),
           loc = as.numeric(inla.posterior.sample.eval("location_detect1", samples)))

 mapply(function(x,y){pnorm(q = ((seq(-4,3,0.1)/qexppnorm(x, rate = 1) + y)))}, x = df$scale, y = df$loc) %>%
  as.data.frame() %>%
  mutate(x = seq(-4,3,0.1)) %>%
  pivot_longer(-x) %>%
  group_by(x) %>% summarise(m = median(value),
                            q1 = quantile(value, 0.025),
                            q2 = quantile(value, 0.975)) %>%
  ggplot() + geom_line(aes(x,m)) +
  geom_ribbon(aes(x,ymin =  q1, ymax = q2), alpha = 0.2) +
  geom_line(data = data.frame(x = seq(-4,3,0.1)) %>%
              mutate(y = pnorm(seq(-4,3,0.1)/qexppnorm(scale_detect1, rate = 1)) +
          location_detect1),
            aes(x,y),color = "red") +
  coord_cartesian(ylim=c(0,1))


## ----res-ferry,fig.width = 4, fig.height= 4-------------------------
# detection function for ferry data
detection_ferry = function(...)
{
  sig = qexppnorm(sig_detect2, rate = rate_detect2)
  return(sig)
}


scale2 = inla.posterior.sample.eval(fun = detection_ferry,samples = samples)
data.frame(apply(scale2, 2, function(x) exp(-0.5*(seq(0,7,0.1)/ x)^2))) %>%
  mutate(x = seq(0,7,0.1)) %>%
  pivot_longer(-x) %>%
  group_by(x) %>% summarise(m = median(value),
                            q1 = quantile(value, 0.025),
                            q2 = quantile(value, 0.975)) %>%
  ggplot() + geom_line(aes(x,m)) +
  geom_ribbon(aes(x,ymin =  q1, ymax = q2), alpha = 0.2) +
  geom_line(data = data.frame(x = seq(0,7,0.1)) %>%
              mutate(y = exp(-0.5*(x/sig_detect2 )^2)),
            aes(x,y),color = "red")



## ----res-field,fig.width = 7, fig.height= 7-------------------------
pxl = pixels(mesh ,nx = 200, ny = 200, mask =poly_sp)

samples_fit = generate(fit,  data = pxl, ~ exp(Intercept + depth +  SPDE))

## compute the posterior median of the intensity
p1 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = apply(samples_fit,1,median)) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0) +
  ggtitle("post median")


## compute the  RWPCI of the intensity
p2 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           m = apply(samples_fit,1,median),
           q1 = apply(samples_fit,1,quantile, 0.25),
           q2 = apply(samples_fit,1,quantile, 0.75)) %>%
  mutate(z = (q2-q1)/m) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0)+
  ggtitle("RWPCI")


p1+p2


## ----res-field2,fig.width = 7, fig.height= 7------------------------
ss = sample(1:100,4)
p1 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = samples_fit[,ss[1]]) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0) 
p2 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = samples_fit[,ss[2]]) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0) 

p3 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = samples_fit[,ss[3]]) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0)

p4 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = samples_fit[,ss[4]]) %>%
  ggplot() + geom_tile(aes(x,y,fill = z)) + scale_fill_scico() +
  theme_map + theme(legend.position = "none") +
  geom_sf(data = poly, alpha = 0) 


p5 = ggplot() + gg(r3, aes(color = layer)) + scale_color_scico() +
 theme_map + geom_sf(data = poly, alpha = 0) + theme(legend.position = "none") +
  ggtitle("True intensity")


p1+p2+p4+p4 + p5  + plot_layout(ncol = 2)

