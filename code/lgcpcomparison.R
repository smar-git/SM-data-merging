
# Load libraries ----------------------------------------------------------


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


my.pixels = function (shape, nx = 150, ny = 150, mask = TRUE) 
{
  
  if (length(nx) == 1) {
    xx = st_coordinates(shape)[,1]
    x <- seq(min(xx), max(xx), length = nx)
  }
  else {
    x <- nx
  }
  if (length(ny) == 1) {
    yy = st_coordinates(shape)[,2]
    y <- seq(min(yy), max(yy), length = ny)
  }
  else {
    y <- ny
  }
  lattice <- INLA::inla.mesh.lattice(x = x, y = y)
  pixels <- data.frame(x = lattice$loc[, 1], y = lattice$loc[, 
                                                             2])
  coordinates(pixels) <- c("x", "y")
  pixels <- SpatialPixels(pixels, 
                          proj4string = inla.sp_get_crs(as(shape,"Spatial")))
  
  if (mask) {
    oo = over(pixels, as(shape,"Spatial"))
    if(!is.null(dim(oo)))
      inside = oo[,1]
    else inside = oo
    pixels <- pixels[!is.na(inside)]
  }
  
  pixels
}


f.fill = function(x,y, covar) 
{
  # turn coordinates into SpatialPoints object:
  spp <- SpatialPoints(data.frame(x=x,y=y),
                       proj4string = inla.sp_get_crs(as(poly,"Spatial"))) 
  # Extract elevation values at spp coords, from our elev SpatialGridDataFrame
  v <- over(spp, covar) 
  v[is.na(v)] <- 0 # NAs are a problem! Remove them
  return(v[,1])
}



# read shapefile and create mesh -------------------------------------------

poly = read_sf("shapefiles/shape_interest/shape_interest.shp")
depth =  read.table("SIMULATIONS/depth.csv") 
win <- as.owin(poly)


## from the sp vignette:
l1 <- cbind(c(724.3269, 626.9514), c(4676.017, 4593.306))
l2 <- cbind(c(787.5578, 733.6008 ), c(4612.691, 4449.854 ))
l3 <- cbind(c(900.5304, 786.2932), c(4570.690, 4442.100))

Sl1 <- Line(l1)
Sl2 <- Line(l2)
Sl3 <- Line(l3)

S1 <- Lines(list(Sl1), ID = "a")
S2 <- Lines(list(Sl2), ID = "b")
S3 <- Lines(list(Sl3), ID = "c")

Sl <- SpatialLines(list(S1, S2, S3),
                   proj4string = inla.sp_get_crs(as(poly,"Spatial")))
df <- data.frame(len = sapply(1:length(Sl), function(i) rgeos::gLength(Sl[i, ])))
rownames(df) <- sapply(1:length(Sl), function(i) Sl@lines[[i]]@ID)


## SpatialLines to SpatialLinesDataFrame
ferry_sp <- SpatialLinesDataFrame(Sl, data = df)
ferry_sp$weights = 12

ggplot() + geom_sf(data = poly) + gg(ferry_sp)


## CREATE THE MESH
boundary = as(st_simplify(poly, dTolerance = 1.3),"Spatial") 
max.edge = 10#15
bound.outer = 50
mesh = inla.mesh.2d(boundary = boundary,
                          max.edge = c(1,4)*max.edge,
                          # - use 5 times max.edge in the outer extension/offset/boundary
                          crs = inla.sp_get_crs(boundary),
                          offset = c(max.edge, bound.outer))

A = inla.spde.make.A(mesh = mesh, loc = as.matrix(depth[,c(1,2)]))

ggplot() + geom_sf(data = poly) +
  gg(mesh)
# simulate point process -------------------------------------------------

# intensity zfunction
beta0 = -4
beta_depth = -1.5
lambda = beta0 + beta_depth * depth$z

linear_SPDF =  SpatialPixelsDataFrame(points = cbind(depth$x,depth$y), 
                                                  data = data.frame(surf = exp(lambda)),
                                                  proj4string = inla.sp_get_crs(as(poly,"Spatial")))

# parameters for the Gaussian field
sigma2x <- 1^2
range <- 100

# Simulate the sightings process
im = as.im(data.frame(x = depth$x, y = depth$y, z = lambda), W = win)

data <- rLGCP(model="matern", 
              mu=im, 
              var=sigma2x, 
              scale=range/sqrt(8), 
              nu=1, 
              win = win,
              saveLambda = TRUE)
data.sp <- SpatialPoints(cbind(data$x,data$y),
                         proj4string = inla.sp_get_crs(as(poly,"Spatial")))

Lambda = raster::raster(attributes(data)$Lambda)
r1 = raster::mask(Lambda, poly)
r2 = raster::crop(r1, poly)


# projection matrix for the datapoints
A_data = inla.spde.make.A(mesh = mesh, loc = cbind(data$x, data$y))
plot(data.sp)



plot(r2)
plot(data.sp, add = T)
plot(ferry_sp, add = T)


m <- rgeos::gDistance(data.sp, ferry_sp, byid=TRUE)
data.frame(x = coordinates(data.sp)[,1],
           y = coordinates(data.sp)[,2],
           z = apply(m,2,min)) %>%
  ggplot()  + geom_point(aes(x,y,color = z)) + scale_color_scico()

    
# simulate observation processes (thinning) -------------------------------

surf1 = function(x,y)
{
  z =  -.5 * (x-557)/(970-557) + cos((y-4379)/(4790-4379) * 8)
  return(z)
}

surf2 = function(x,y)
{
  
  d = SpatialPoints(cbind(x,y),
                    proj4string = inla.sp_get_crs(as(poly,"Spatial")))
  d1 = rgeos::gDistance(d, ferry_sp, byid=TRUE)
  z = apply(d1,2,min)
  #z = ((x-754)/(970-754))^2 + ((y-4619)/(4790-4619))^2
  return(z)
}

#probability of keeping a sightings

## pnorm
prob1 = data.frame(x = depth$x, y = depth$y) %>%
  mutate(surf = surf1(x,y)) %>%
  mutate(prob = pnorm(surf, mean = 0, sd = 1)) %>%
  ggplot() + geom_tile(aes(x,y,fill = prob)) +
  scale_fill_scico(lim = c(0,1)) +
  geom_sf(data = poly, alpha = 0) +
  ggtitle("detection prob \\
          for observator 1")

# 
prob2 = data.frame(x = depth$x, y = depth$y) %>%
  mutate(surf = surf2(x,y)) %>%
  mutate(prob = exp(-0.5*(surf/ 1.5)^2)) %>%
  ggplot() + geom_tile(aes(x,y,fill = prob)) +
  scale_fill_scico() +
  geom_sf(data = poly, alpha = 0)+
  ggtitle("detection prob \\
          for observator 2")
prob1 + prob2

## thinning the datapoits
all_data = data.frame(x = data$x,
                      y = data$y)
  
 

## parameters for the detection function 1
scale_detect1 = 1
location_detect1 = 0

data1 = all_data %>% 
  mutate(surf = surf1(x,y)) %>%
  mutate(pdetect = pnorm(surf, mean = 0, sd = 1),
         u = runif(data$n)) %>%
  mutate(detect = u<pdetect) %>%
  dplyr::filter(detect)

## parameters for the detection function 2
sig_detect2 = 2

data2 = all_data %>% 
  mutate(surf = surf2(x,y)) %>%
  mutate(pdetect = exp(-0.5*(surf/ sig_detect2)^2),
         u = runif(data$n)) %>%
  mutate(detect = u<pdetect) %>%
  dplyr::filter(detect)


ggplot() + geom_point(data = all_data, aes(x,y)) +
  geom_point(data = data1, aes(x,y), color = "red")

ggplot() + geom_point(data = all_data, aes(x,y)) +
  gg(ferry_sp) + 
  geom_point(data = data2, aes(x,y), color = "red")


# create the datasets as SP and all maps as SPDF --------------------------


data1_sp = SpatialPointsDataFrame(cbind(x = data1$x, y = data1$y), 
                                           data = data.frame(d = rep(1,length(data1$x))),
                                           proj4string = inla.sp_get_crs(as(poly,"Spatial")))

data2_sp = SpatialPointsDataFrame(cbind(x = data2$x, y = data2$y), 
                                  data = data.frame(d = rep(1,length(data2$x))),
                                  proj4string = inla.sp_get_crs(as(poly,"Spatial")))

d1 = rgeos::gDistance(data2_sp, ferry_sp, byid=TRUE)
d1 = apply(d1,2,min)
data2_sp$dist_obs = d1

df = data.frame(x = depth$x, y = depth$y) %>%
  mutate(surf = surf1(x,y))
surf1_SPDF = SpatialPixelsDataFrame(points = cbind(df$x,df$y), 
                                    data = data.frame(surf = df$surf),
                                    proj4string = inla.sp_get_crs(as(poly,"Spatial")))

df = data.frame(x = depth$x, y = depth$y) %>%
  mutate(surf = surf2(x,y))
surf2_SPDF = SpatialPixelsDataFrame(points = cbind(df$x,df$y), 
                                    data = data.frame(surf = df$surf),
                                    proj4string = inla.sp_get_crs(as(poly,"Spatial")))


depth_SPDF = SpatialPixelsDataFrame(points = cbind(depth$x,depth$y), 
                                    data = data.frame(depth = depth$z),
                                    proj4string = inla.sp_get_crs(as(poly,"Spatial")))


# Define the two detection functions --------------------------------------


qexppnorm <- function(x, rate) {
  qexp(pnorm(x, lower.tail = FALSE, log.p = TRUE),
       rate = rate, 
       log.p = TRUE, 
       lower.tail = FALSE)
}




log_detect1 = function(x,y, t1, t2)
{
  dens = f.fill(x,y,surf1_SPDF)
  pnorm(dens/qexppnorm(t1, rate = 1) +
          t2,
        log.p = TRUE)
}

log_detect2 = function(x,y, sig)
{ 
  dist1 = f.fill(x,y,surf2_SPDF)
  -0.5*(dist1/ sig)^2
}





# set up model ------------------------------------------------------------

matern <-  inla.spde2.pcmatern(mesh,
                                   prior.sigma = c(1, 0.01),
                                   prior.range = c(150, 0.5))

cmp = ~ 
  SPDE(main = coordinates, model = matern) +
  depth(depth_SPDF) +
  Intercept(1)+
  sig_detect2(main = 1, model = "linear",
             mean.linear = 0,
             prec.linear = 1) +
  scale_detect1(1, prec.linear = 1)+
  location_detect1(1, prec.linear = 10) 

formula_data1 = coordinates   ~ 
  Intercept +
  depth  + 
  SPDE +
  log_detect1(x,y,location_detect1, scale_detect1)
  
formula_data2 = coordinates + dist_obs   ~ 
  Intercept +
  depth  + 
  SPDE +
  log_detect2(x,y,qexppnorm(sig_detect2,
                                  rate = 1))


lik_data1 <- like("cp",
                        formula = formula_data1,
                        samplers = as(poly, "Spatial"), 
                        domain = list(coordinates = mesh),
                        data = data1_sp)


lik_data2 <- like("cp",
                  formula = formula_data2,
                  samplers = ferry_sp,
                  domain = list(coordinates = mesh,
                                dist_obs = inla.mesh.1d(seq(0, 6, length.out = 5),
                                                        degree = 1)),
                  data = data2_sp)

c.c <- list(cpo = TRUE,
            dic = TRUE,
            waic = TRUE,
            config = TRUE)

bru_options_set(bru_verbose = 1)
fit= bru(components = cmp,  
         lik_data1,
         lik_data2,
         options = list(verbose = F,
                        num.threads = "1:1",
                        control.compute = c.c))

# check results -----------------------------------------------------------

# Linear effects
par(mfrow = c(1,2))
plot(fit$marginals.fixed$Intercept , type = "l")
abline(v = beta0)
plot(fit$marginals.fixed$depth , type = "l")
abline(v = beta_depth)


# parameters of the SPDE model
par(mfrow =c(1,2))
plot(fit$marginals.hyperpar$`Range for SPDE` , type = "l")
abline(v = range)
plot(fit$marginals.hyperpar$`Stdev for SPDE` , type = "l")
abline(v = sqrt(sigma2x))


# parameters of the detection functions
detection_ferry = function(...)
{
  sig = qexppnorm(sig_detect2, rate = 1)
  return(sig)
}

scale_det1 = function(...)
{
  qexppnorm(scale_detect1, rate = 1) 
}

samples = inla.posterior.sample(1000, fit)
scale2 = inla.posterior.sample.eval(fun = sd_det2,samples = samples)
data.frame(apply(scale2, 2, function(x) exp(-0.5*(seq(0,6,0.1)/ x)^2))) %>% 
  mutate(x = seq(0,6,0.1)) %>%
  pivot_longer(-x) %>%
  group_by(x) %>% summarise(m = median(value), 
                            q1 = quantile(value, 0.025),
                            q2 = quantile(value, 0.975)) %>%
  ggplot() + geom_line(aes(x,m)) +
  geom_ribbon(aes(x,ymin =  q1, ymax = q2), alpha = 0.2) +
  geom_line(data = data.frame(x = seq(0,6,0.1), y = exp(-0.5*(seq(0,6,0.1)/ 1.5)^2)),
            aes(x,y),color = "red")


pxl = my.pixels(shape = poly ,nx = 200, ny = 200, mask =T)

samples_fit = predict(fit,  data = pxl, ~ exp(Intercept + depth +  SPDE))

p1 = data.frame(x = coordinates(pxl)[,1],
           y = coordinates(pxl)[,2],
           z = samples_fit$median) %>%
  ggplot()  + geom_tile(aes(x,y,fill  = z)) +
  geom_sf(data = poly, alpha =0) + 
  scale_fill_scico()

p2 = data.frame(x = coordinates(pxl)[,1],
                y = coordinates(pxl)[,2],
                z = samples_fit$cv) %>%
  ggplot()  + geom_tile(aes(x,y,fill  = z)) +
  geom_sf(data = poly, alpha =0) + 
  scale_fill_scico()

r3 <- rasterToPoints(r2, spatial = TRUE)
p3 = ggplot() + gg(r3, aes(color = layer)) + scale_color_scico() +
  geom_sf(data = poly, alpha =0) 
  
p1 + p2 + p3
