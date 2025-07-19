rm(list=ls()) #----Clear the workspace
packages <- c("raster", "haven", "sf","sp", "tmap","tmaptools","tidyverse","terra",
              "lattice", "gridExtra", "devtools", "rlang", "viridis", "spdep",
              "car", "MASS", "maps", "spData", "ceramic", "basemaps", "ggmap")
if(length(setdiff(packages, rownames(installed.packages()))) > 0) {
  install.packages(setdiff(packages, rownames(installed.packages()))) }


#Installing INLA!!
if(length(setdiff("INLA", rownames(installed.packages()))) > 0){
  install.packages("INLA", type="binary", repos=c(getOption("repos"), INLA="https://inla.r-inla-download.org/R/stable"), dep=TRUE)
}
library(INLA)
lapply(packages, library, character.only = TRUE) ##--access the libraries

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

# set directories
path <- "C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts"
data_path <- paste0(path, "/data")
out_path <- paste0(path, "/output")

#
# install.packages("devtools")
#devtools::install_github("wpgp/jollofR")
library(jollofR)
#---Select AdmUnit Size

##--Cluster
cl = 30#--900 clusters with 36 clusters per province
clust <- as(raster(nrow=cl, ncol=cl, xmn=0, xmx=1, ymn=0, ymx=1), "SpatialPolygons")
proj4string(clust) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
plot(clust, lwd=1, col="green")


#-- grid cells
pg = 120 #--14400 prediction grid cells with 16 grid cells per cluster
#-- 16 times 36 (576) grid cells in each of the 25 provinces
grid = raster(nrow=pg, ncol=pg, xmn=0, xmx=1, ymn=0, ymx=1) #--prediction grid
res(grid) = c(0.0083, 0.0083) #--approx 1km by 1km resolution
ncell(grid)

##--Vectorize the prediction raster
coords.grid = xyFromCell(grid, 1:ncell(grid), spatial=FALSE)


##---Visualise
plot(coords.grid, col="green",
     pch=16, cex.axis=1.5)
#plot(clust, lwd=1.5)
plot(clust, add=T, lwd=1.5)



##---Convert all polygons and grid cells to points
coords.clust <- coordinates(clust)
coords.grid <- coordinates(grid)




####----Generate id of grid cells within the Clusters
c.g=SpatialPoints(coords.grid)
proj4string(c.g) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
gc <- rep(NA, nrow(coords.grid))
for(i in 1:length(clust)){
  gc[as.vector(which(!is.na(over(c.g, clust[i]))))] <- i
}
grid4clust <- gc
table(grid4clust)
unique(grid4clust)



##--rename
grid_id <- 1:ncell(grid)
clust_id <- grid4clust


##----Check data for consistency
df1 <- data.frame(grd_id = grid_id,
                  clst_id = clust_id,
                  lon = coords.grid[,1],
                  lat = coords.grid[,2])
df1[order(df1$clst_id),] #--passed
plot(df1$lon, df1$lat)
par(mfrow=c(1,1))

dat_all <- df1
dat.grid <- df1



####---Model fit metrics function
model_metrics <- function(obs, pred)
{
  residual = pred - obs
  MAE = mean(abs(residual), na.rm=T)#MAE
  MAPE = (1/length(obs))*sum(abs((obs-pred)/obs))*100#MAPE
  MSE = mean(residual^2, na.rm=T)
  RMSE = sqrt(MSE)
  BIAS = mean(residual, na.rm=T)
  corr = cor(obs[!is.na(obs)],pred[!is.na(obs)])

  output <- list(MAE  = MAE,
                 RMSE = RMSE,
                 corr = corr)
  return(output)
}
#model_metrics(obs, pred)


datt <- dat.grid
###----Build non-hull mesh (for combined)
coords <- cbind(datt$lon, datt$lat)
plot(coords)
#library(jollofR)
non_convex_bdry1 <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
plot(non_convex_bdry1$loc)
points(datt$lon, datt$lat)

#-------Extract boundary for mesh
mesh <- inla.mesh.2d(boundary = non_convex_bdry1, max.edge=c(0.02,0.1),
                     offset = c(0.02, 0.04),
                     cutoff = 0.005)
par(mfrow=c(1,1))
plot(mesh)
plot(clust, add=T, lwd=1.5)
plot(mesh, add=T)

mesh$n # 7250 mesh nodes


##---Specify SPDE parameters
r0 <- 0.45 #--range
nu <- 1  #--smooth parameter
sigma0 <- 1 #--marginal variance  C(0.01, 0.1, 1)
kappa0 <- sqrt(8*nu)/r0 #--scale parameters
(tau0 <- 1/(sqrt(4*pi)*kappa0*sigma0)) #--precision

#---SPDE
spde <- inla.spde2.matern(mesh,
                          B.tau = matrix(c(log(tau0), -1, +1), nrow=1, ncol=3),
                          B.kappa = matrix(c(log(kappa0), 0, -1), nrow=1, ncol=3),
                          theta.prior.mean = c(0,0),
                          theta.prior.prec = c(0.1, 0.1))


#--Precision Matrix
Q <- inla.spde2.precision(spde=spde, theta=c(0,0))


#---Simulate the GMRF
sam <- as.vector(inla.qsample(
  n = 1, Q = Q, seed=100))


N = nrow(datt) # number of EAs
###---Build projector matrix A
A <- inla.spde.make.A(mesh=mesh, loc=coords);dim(A)
spat <- as.vector(A %*% sam)# Spatially correlated random effects


##----specify the observation indices for estimation
iset <- inla.spde.make.index(name = "s", spde$n.spde)


set.seed(189)
b0 = 3.5 #intercept
b1 = 1.781 # coefficient of the first covariate
b2 = 0.891 # coefficient of the second covariate
b3 = 1.158 # coefficient of the 3rd covariate


# For building count
b01 = 1.5 #intercept
b11 = 0.921 # coefficient of the first covariate
b21 = 0.85 # coefficient of the second covariate
b31 = 0.68# coefficient of the 3rd covariate

EA_ID= paste0("EA", 1:N)

# Simulate the standardized covariates
stdize <- function(x){return(stdize <- (x-min(x, na.rm=T))/(max(x, na.rm=T)-min(x, na.rm=T)))} # normalize

# Simulate covariates
(x1 <- stdize(rnorm(N, 20, 16)))
(x2 <- stdize(rnorm(N, 5, 4)))
(x3 <- stdize(runif(N, 0, 100)+0.35*x2))


# Initial data frame
(sim.dat <- data.frame(EA_ID, x1, x2, x3))
sigma_e <- 0.001 ##--
eps <- rnorm(nrow(sim.dat), 0, sigma_e) # IID
table(sim.dat$set_typ <- sample(c("urban", "rural"), nrow(sim.dat), rep=T,
                                prob=c(0.35, 0.65)))

ms <- length(unique(sim.dat$set_typ))
set_ranef <- rnorm(rep(0,ms), diag(c(0.00021, 0.00015)))
sim.dat$set_ranef <- sim.dat$set_typ
sim.dat$set_ranef[sim.dat$set_typ == "urban"] = set_ranef[1]
sim.dat$set_ranef[sim.dat$set_typ == "rural"] = set_ranef[2]
sim.dat$set_ranef <- round(as.numeric(sim.dat$set_ranef),4)

# simulate EA total pop as a Poisson rv
# Add settlement types



(mu = with(sim.dat, b0 + b1*x1 + b2*x2  + b3*x3 + set_ranef + spat + eps))
(lam <- exp(mu)) # average EA pop


# building count
mm <- stdize(lam)
(mu2 = with(sim.dat, b01 + b11*x1 + b21*x2  + b31*x3  + mm + eps))
(lam2 <- exp(mu2)) # average EA pop

boxplot(sim.dat$total <- rpois(N, lambda=lam)) # pop for each EA
sim.dat$total

###
boxplot(sim.dat$bld <-  rpois(N, lambda=lam2)) # building count
sim.dat$bld
dat <- cbind(datt, sim.dat)


plot(dat$total, dat$bld)

# add age groups

  group <- c(8, 16, 24, 32, 40)# Number of age groups

  #tryCatch(
  for(g in 1:length(group))
  {
     m <- group[g]
     
    result_path1 <- paste0(out_path,"/outputs_for_", group[g],"_age_groups")
    if (file.exists(result_path1)){
      setwd(file.path(result_path1))
    } else {
      dir.create(file.path(result_path1))
      setwd(file.path(result_path1))
    }

  # Proportion of groups ( age groups)
  (rng <- c(sort(sample(1:((m/2)-2), rep=T)),(m/3)-1,m/3,m/4, ((m/3)-2),rev(sort(sample(1:((m/2)-2), rep=T)/3))))
  sum(a_prop <- round(rng/sum(rng),4)) # calculate age group proportion

  grp_mat <- matrix(0, N, m)
  for(i in 1:N)
  {
    sum(grp_prop <- a_prop)
    grp_mat[i,] <- t(rmultinom(n=1, size = sim.dat$total[i], prob = as.vector(grp_prop)))
  }
  (grp_pop <- data.frame(grp_mat))
  grp_pop$total2 <- apply(grp_pop, 1, sum) # add row sums


  N <- nrow(dat)
  f_mat <- matrix(0, N, m) # for female population
  m_mat <- matrix(0, N, m) # for male population
  for(i in 1:m)
  {
    fprop <- round(runif(N, 0.4, 0.65),2) # set the proportion to vary between 0.4 and 0.65 ()
    f_mat[,i] <- round(grp_mat[,i]*fprop) # female population
    m_mat[,i] <- grp_mat[,i] - f_mat[,i] # male population
  }
  (f_pop <- data.frame(f_mat)) # female population
  (m_pop <- data.frame(m_mat)) # male population

  f_pop[8,4]+m_pop[8,4];grp_pop[8,4]# check!



  # rename variables
  colnames(f_pop)[1:m] <- paste0("fage_", 1:m) # females
  colnames(m_pop)[1:m] <- paste0("mage_", 1:m)  # males
  colnames(grp_pop)[1:m] <- paste0("age_", 1:m) # all ages



  # join the datasets
  dat2 <- cbind(dat, grp_pop, f_pop, m_pop)



  # add  educational level

  # Create educational level data
  #f <- 4 # 4 categories

  edu_prop <- c(0.15, 0.3, 0.35,0.2) # proportion of EDUCATIONAL LEVELS groups to add
  # no edu = 0.15
  # primary = 0.3
  #secondary = 0.35
  # higher = 0.2
   dat2$edu_no <- round(dat2$total*0.15)
   dat2$edu_prim <- round(dat2$total*0.3)
   dat2$edu_sec <- round(dat2$total*0.35)
   dat2$edu_high <- round(dat2$total*0.2)



  #save simulated grid data
  dat2$admin_id <- dat2$clst_id
  #saveRDS(dat, paste0(data_path, "/toy_grid_data_full.rds"))


  #names(dat)
  vars2include <- c("admin_id", "grd_id", "total", "bld", "set_typ", "lon", "lat") # select fewer variables
  dtt <- dat2[, vars2include]
  #saveRDS(dtt, paste0(data_path, "/toy_grid_data.rds"))

  # Visualise the pyramid
  #devtools::install_github("wpgp/jollofR")
  #library(dplyr)
  #library(jollofR)
  female_pop <- data.frame(dat2 %>% dplyr::select(starts_with("fage_"))) # extract females age data
  names(female_pop) <- paste0("pp_", names(female_pop)) # rename the variables by adding "pp_" as suffix to the existing names

  male_pop <- data.frame(dat2 %>% dplyr::select(starts_with("mage_")))# extract males age data
  names(male_pop) <- paste0("pp_", names(male_pop))# rename the variables by adding "pp_" as suffix to the existing names

 # pyramid(female_pop,male_pop) # make the observed pyramid plot
  f_mat <- female_pop
  age_classes <- names(f_mat)
  f_mat$id <- 1:nrow(f_mat)
  (female_pop <- reshape2::melt(f_mat, id = c("id"), value.name = "Population",
                                variable.name = "Age"))
  female_pop$Age <- factor(female_pop$Age)
  levels(female_pop$Age) <- gsub("pp_fage_", "", age_classes)
  female_pop$Gender <- rep("Female", nrow(female_pop))
  m_mat <- male_pop
  age_classes <- names(m_mat)
  m_mat$id <- 1:nrow(m_mat)
  (male_pop <- reshape2::melt(m_mat, id = c("id"), value.name = "Population",
                              variable.name = "Age"))
  male_pop$Age <- factor(male_pop$Age)
  levels(male_pop$Age) <- gsub("pp_mage_", "", age_classes)
  male_pop$Gender <- rep("Male", nrow(male_pop))
  dim(pop_pyramid <- rbind(female_pop, male_pop))
  population_pyramid1 <- ggplot(pop_pyramid, aes(x = Age, fill = Gender,
                                                 y = ifelse(test = Gender == "Male", yes = -Population,
                                                            no = Population))) + geom_bar(stat = "identity") +
    scale_y_continuous(labels = abs) + labs(x = "Age", y = "Population",
                                            fill = "Gender") + theme_minimal() + coord_flip()

  pyramid1 <- ggpubr::ggpar(population_pyramid1, ylab = "Population Count",
                            xlab = "Age (years)", legend = "right", legend.title = "Gender",
                            size = 22, font.legend = c(16), palette = "lancet", font.label = list(size = 15,
                            face = "bold", color = "red"), font.x = c(16), font.y = c(16),
                            font.main = c(14), font.xtickslab = c(14), font.ytickslab = c(16))
  print(pyramid1)

  ##
  # Aggregate to culster level

  names(dat2)

  ## aggeregate covariates, building counts and pop counts
  # install.packages("DescTools")
  dat_clst1 <- dat2 %>% group_by(clst_id) %>%
    reframe(x1 = mean(x1, na.rm=T),
            x2 = mean(x2, na.rm=T),
            x3 = mean(x3, na.rm=T),
            bld = sum(bld, na.rm=T),
            edu_no = sum(edu_no, na.rm=T),
            edu_prim = sum(edu_prim, na.rm=T),
            edu_sec = sum(edu_sec, na.rm=T),
            edu_high = sum(edu_high, na.rm=T),
            set_typ = DescTools::Mode(as.character(set_typ), na.rm=TRUE)[1], # settlement type
            total = sum(total, na.rm=T))

  dim(dat_clst1)
  table(dat_clst1$set_typ)

  dat_clst1$set_typ<- factor(dat_clst1$set_typ)


  # Aggregate age data
  clst_age <- dat2 %>% group_by(clst_id) %>%
    summarise_at(vars(age_1:paste0("mage_", group[g])), sum, na.rm=T)
  
  
  # combined all data
   clst_dat <- merge(dat_clst1, clst_age, by = "clst_id")


   # Add lon lat to the cluster data
   clst_dat$lon <- coords.clust[,1]
   clst_dat$lat <- coords.clust[,2]

   plot(clst_dat$lon, clst_dat$lat)

   # Build the mesh, Projection matix and spde at admin level
   coords2 <- cbind(clst_dat$lon, clst_dat$lat)
   plot(coords2)

   #library(jollofR)
   non_convex_bdry2 <- inla.nonconvex.hull(coords2, -0.03, -0.05, resolution = c(100, 100))
   plot(non_convex_bdry2$loc)


   #-------Extract boundary for mesh
   mesh2 <- inla.mesh.2d(boundary = non_convex_bdry2, max.edge=c(0.02,0.1),
                        offset = c(0.02, 0.04),
                        cutoff = 0.005)
   par(mfrow=c(1,1))
   plot(mesh2)
   plot(clust, add=T, lwd=1.5)
   plot(mesh2, add=T)

   mesh2$n # 6892 mesh nodes


   ##---Create the SPDE
   spde2 <- inla.spde2.matern(mesh2, alpha=2)

   ###---Build projector matrix A
   A2 <- inla.spde.make.A(mesh=mesh2, loc=coords2);dim(A2)

   ##----specify the observation indices for estimation
   iset2 <- inla.spde.make.index(name = "s", spde2$n.spde)



  # Make pyramid plots
   # check with the grid pyramid
   female_pop2 <- data.frame(clst_dat %>% dplyr::select(starts_with("fage_"))) # extract females age data
   names(female_pop2) <- paste0("pp_", names(female_pop2)) # rename the variables by adding "pp_" as suffix to the existing names

   male_pop2 <- data.frame(clst_dat %>% dplyr::select(starts_with("mage_")))# extract males age data
   names(male_pop2) <- paste0("pp_", names(male_pop2))# rename the variables by adding "pp_" as suffix to the existing names

   #pyramid(female_pop,male_pop) # make the observed pyramid plot
   # pyramid(female_pop,male_pop) # make the observed pyramid plot
   f_mat2 <- female_pop2
   age_classes <- names(f_mat2)
   f_mat2$id <- 1:nrow(f_mat2)
   (female_pop2 <- reshape2::melt(f_mat2, id = c("id"), value.name = "Population",
                                 variable.name = "Age"))
   female_pop2$Age <- factor(female_pop2$Age)
   levels(female_pop2$Age) <- gsub("pp_fage_", "", age_classes)
   female_pop2$Gender <- rep("Female", nrow(female_pop2))
   m_mat2 <- male_pop2
   age_classes <- names(m_mat2)
   m_mat2$id <- 1:nrow(m_mat2)
   (male_pop2 <- reshape2::melt(m_mat2, id = c("id"), value.name = "Population",
                               variable.name = "Age"))
   male_pop2$Age <- factor(male_pop2$Age)
   levels(male_pop2$Age) <- gsub("pp_mage_", "", age_classes)
   male_pop2$Gender <- rep("Male", nrow(male_pop2))
   dim(pop_pyramid2 <- rbind(female_pop2, male_pop2))

   population_pyramid2 <- ggplot(pop_pyramid2, aes(x = Age, fill = Gender,
                                                  y = ifelse(test = Gender == "Male", yes = -Population,
                                                 no = Population))) + geom_bar(stat = "identity") +
     scale_y_continuous(labels = abs) + labs(x = "Age", y = "Population",
                                             fill = "Gender") + theme_minimal() + coord_flip()

   pyramid2 <- ggpubr::ggpar(population_pyramid2, ylab = "Population Count",
                             xlab = "Age (years)", legend = "right", legend.title = "Gender",
                             size = 22, font.legend = c(16), palette = "lancet", font.label = list(size = 15,
                            face = "bold", color = "red"), font.x = c(16), font.y = c(16),
                             font.main = c(14), font.xtickslab = c(14), font.ytickslab = c(16))
   print(pyramid2)


   table(pop_pyramid2$Age)
   # Save the cluster data
   clst_dat$admin_id <- clst_dat$clst_id
   #write.csv(clst_dat, paste0(data_path, "/toy_admin_data_full.csv"), row.names=F)

   ## Apply missingness


   # create missingness in the data

  # missn = c(0.1, 0.3, 0.5, 0.7,0.9)

   sample_size <- c(900, 300, 100, 30, 10)
   metrics <- list()

   for(k in 1:length(missn))
   {
     #k=1
     #print(k)
    # result_path2 <- paste0(result_path1,"/for_", missn[k]*100, "%", "_missingness")
     result_path2 <- paste0(result_path1,"/for_", sample_size[k], "_sample_size")
     if (file.exists(result_path2)){
       setwd(file.path(result_path2))
     } else {
       dir.create(file.path(result_path2))
       setwd(file.path(result_path2))
     }
     #print(paste0(m," age groups", " at ", missn[k]*100, "%", " missingness"))#---
     
     print(paste0(m," age groups", " at ", sample_size[k], "_sample_size"))#---

  # create missingness in the data
     #miss = missn[k]
     ssize = sample_size[k]
     M <- nrow(clst_dat)
     #ind <- sample(1:M, M*miss)
     ind <- sample(1:M, ssize)

      grp_popm <-clst_age
      grp_popm[-ind,-1] = NA #allow some missing  proportion




      clst_datm <- merge(dat_clst1,  grp_popm, by = "clst_id")

      clst_datm[-ind, c("edu_no", "edu_prim", "edu_sec", "edu_high")] = NA


      # Add lon lat to the cluster data
      clst_datm$lon <- coords.clust[,1]
      clst_datm$lat <- coords.clust[,2]

      clst_datm$total <- round(clst_datm$total)
  

      # Make pyramid plots
      # check with the grid pyramid

      female_pop3 <- data.frame(clst_datm %>% dplyr::select(starts_with("fage_"))) # extract females age data
      names(female_pop3) <- paste0("pp_", names(female_pop3)) # rename the variables by adding "pp_" as suffix to the existing names

      male_pop3 <- data.frame(clst_datm %>% dplyr::select(starts_with("mage_")))# extract males age data
      names(male_pop3) <- paste0("pp_", names(male_pop3))# rename the variables by adding "pp_" as suffix to the existing names

      pyramid(female_pop3,male_pop3) # make the observed pyramid plot

      ### Apply to multiple groups
      # Build data stack

      clst_datm$admin_id <- clst_datm$clst_id

      library(tidyr)
      library(dplyr)

      # Create grid outputs matrices
      admin_data <- clst_datm
      grid_data <- dat
      grid_data$admin_id <- grid_data$clst_id
      # Train model at admin level and predict at grid cell level


      # extract class names for age
      age_dat <-  admin_data %>% dplyr::select(starts_with("age_"))
      pred_dt <- prop_dt <- age_dat
      pred_dtL <- prop_dtL <- age_dat
      pred_dtU <- prop_dtU <- age_dat
      age_classes <-  names(age_dat)
      age_dat$total <- apply(age_dat[,age_classes], 1, sum, na.rm=T)
      age_dat$total[age_dat$total==0] = NA
      covs <- admin_data %>% dplyr::select(starts_with("x"))
      age_dat$admin_id <- admin_data$admin_id
      age_dat$pop <- admin_data$total
      age_dat$set_typ <- factor(admin_data$set_typ)
      age_dat <- cbind(age_dat,covs)

      # female age data
      f_dat = admin_data %>% dplyr::select(starts_with("fage_"))
      f_age_classes <-  names(f_dat)
      f.pred_dt <- f.prop_dt <-m.pred_dt <- m.prop_dt <- f_dat
      f.pred_dtL <- f.prop_dtL <-m.pred_dtL <- m.prop_dtL <- f_dat
      f.pred_dtU <- f.prop_dtU <-m.pred_dtU <- m.prop_dtU <- f_dat

      # educational group data
      edu_dat <-  admin_data %>% dplyr::select(starts_with("edu_"))
      pred_edu_dt <- prop_edu_dt <- edu_dat
      pred_edu_dtL <- prop_edu_dtL <- edu_dat
      pred_edu_dtU <- prop_edu_dtU <- edu_dat
      edu_classes <-  names(edu_dat)
      edu_dat$total <- apply(edu_dat[,edu_classes], 1, sum, na.rm=T)
      edu_dat$total[edu_dat$total==0] = NA
      edu_dat$admin_id <- admin_data$admin_id
      edu_dat$pop <- admin_data$total
      edu_dat <- cbind(edu_dat,covs)


     # grid outputs matrices

      # For age
      male_grd <- fem_grd <- prop_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(age_classes))
      male_grdL <- fem_grdL <- prop_grdL <- male_grd
      male_grdU<- fem_grdU <- prop_grdU <- male_grd
      prop_grd2 <- prop_grd3 <- prop_gr4   <- prop_grd

      male_grdt <- fem_grdt <- pop_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(age_classes))
      male_grdtL <- fem_grdtL <- pop_grdL <- male_grd
      male_grdtU<- fem_grdtU <- pop_grdU <- male_grd
      pop_grd2 <- pop_grd3 <- pop_gr4   <- pop_grd

      # For education
      male_edu_grd <- fem_edu_grd <- prop_edu_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(edu_classes))
      male_edu_grdL <- fem_edu_grdL <- prop_edu_grdL <- male_edu_grd
      male_edu_grdU<- fem_edu_grdU <- prop_edu_grdU <- male_edu_grd
      prop_edu_grd2 <- prop__edu_grd3 <- prop_edu_gr4   <- prop_edu_grd


      male_edu_grdt <- fem_edu_grdt <- pop_edu_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(edu_classes))
      male_edu_grdtL <- fem_edu_grdtL <- pop_edu_grdL <- male_edu_grd
      male_edu_grdtU<- fem_edu_grdtU <- pop_edu_grdU <- male_edu_grd
      pop_edu_grd2 <- pop__edu_grd3 <- pop_edu_gr4   <- pop_edu_grd

      f_dat$admin_id <- admin_data$admin_id

      cov_names <- names(age_dat%>% dplyr::select(starts_with("x")))


      ## fit age models
      for(i in 1:length(age_classes))
      {

        #  i=5
        # Disaggregate Each age group count by sex
        prior.prec <- list(prec = list(prior = "pc.prec",
                                       param = c(1, 0.01))) # using PC prior
        print(paste(paste0("(",i,")"),paste0(age_classes[i], " model is running")))

        age_dat[,colnames(age_dat)[i]] <- round(age_dat[,i])

        covars <- age_dat[,c(cov_names,"set_typ", "admin_id")]; dim(covars) ##---Population density

        # stack for age
        stk <- inla.stack(data=list(y=age_dat[,i], n=age_dat$total), #the response

                          A=list(A2,1),  #the A matrix; the 1 is included to make the list(covariates)

                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset2),  #the spatial index
                                       #the covariates
                                       list(covars)
                          ),
                          #this is a quick name so you can call upon easily
                          tag='est')


        print(form_age <- as.formula(paste0("y", " ~ ",
                                      paste(c("-1","Intercept", cov_names), collapse = " + "),
                                      " +   f(admin_id, model = 'iid', hyper = prior.prec)+
                                      f(set_typ, model = 'iid', hyper = prior.prec) +  f(s, model = spde2)"))) # This is how you include  a typical random effect.
        mod_age <-inla(form_age, #the formula
                    data=inla.stack.data(stk,spde=spde2),  #the data stack
                    family= 'binomial', Ntrials = n,  #which family the data comes from
                    control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                    control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                    verbose = FALSE)

        summary(mod_age)

        indx <-inla.stack.index(stk, "est")$data #---extract the data location indices
        prop_dt[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"mean"]),6)
        prop_dtL[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.025quant"]),6)
        prop_dtU[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.975quant"]),6)

        pred_dt[,i] = round(prop_dt[,i]*age_dat$pop)
        pred_dtL[,i] = round(prop_dtL[,i]*age_dat$pop)
        pred_dtU[,i] = round(prop_dtU[,i]*age_dat$pop)

        summary(mod_age)

        
  ## Route 1: Using grid cells' geospatial covariates values
        # Extract spatial random effect
        Apred <- inla.spde.make.A(mesh=mesh2, loc=coords);dim(Apred)

        # mean
        sfield_nodes <- mod_age$summary.random$s['mean']
        field <- (Apred %*% as.data.frame(sfield_nodes)[, 1])
        # Extract IID
        IID <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$mean[1])
        strf <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$mean[2]) # settlement type random effect (posterior)
        prop_grd[,i] <- round(plogis(mod_age$summary.fixed$mean[1] +
                                       mod_age$summary.fixed$mean[2]*grid_data$x1 +
                                       mod_age$summary.fixed$mean[3]*grid_data$x2 +
                                       mod_age$summary.fixed$mean[4]*grid_data$x3 +
                                       field[,1] +
                                       strf +
                                       IID),4)

        # lower
        sfield_nodesL <- mod_age$summary.random$s['0.025quant']
        fieldL <- (Apred %*% as.data.frame(sfield_nodesL)[, 1])

        # Extract IID
        IIDL <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$`0.025quant`[1])
        strfL <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$`0.025quant`[2])
        prop_grdL[,i] <- round(plogis(mod_age$summary.fixed$`0.025quant`[1] +
                                       mod_age$summary.fixed$`0.025quant`[2]*grid_data$x1 +
                                       mod_age$summary.fixed$`0.025quant`[3]*grid_data$x2 +
                                       mod_age$summary.fixed$`0.025quant`[4]*grid_data$x3 +
                                       fieldL[,1] +
                                       strfL +
                                       IIDL),4)

        # UPPER
        sfield_nodesU <- mod_age$summary.random$s['0.975quant']
        fieldU <- (Apred %*% as.data.frame(sfield_nodesU)[, 1])

        # Extract IID
        IIDU <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$`0.975quant`[1])
        strfU <- rnorm(nrow(grid_data), 0, 1/mod_age$summary.hyperpar$`0.975quant`[2])
        prop_grdU[,i] <- round(plogis(mod_age$summary.fixed$`0.975quant`[1] +
                                       mod_age$summary.fixed$`0.975quant`[2]*grid_data$x1 +
                                       mod_age$summary.fixed$`0.975quant`[3]*grid_data$x2 +
                                       mod_age$summary.fixed$`0.975quant`[4]*grid_data$x3 +
                                       fieldU[,1] +
                                        strfU +
                                       IIDU),4)

        form_sex <- as.formula(paste0(colnames(f_dat)[i], " ~ ",
                                      "1 +   f(admin_id, model = 'iid', hyper = prior.prec)"))# Adding the IID here


        mod_sex  <- inla(form_sex,
                         data = f_dat,
                         family = "binomial", Ntrials = age_dat[,i],
                         control.predictor = list(compute = TRUE),
                         control.compute = list(dic = TRUE, cpo = TRUE)
        )
        f.prop_dt[,i] = round(plogis(mod_sex$summary.linear.predictor$mean),4) # female proportion - mean
        m.prop_dt[,i] = 1- f.prop_dt[,i] # male proportion

        fem_grd[,i] = round(plogis(mod_sex$summary.fixed$mean[1] +
                                     rnorm(nrow(fem_grd), 0, 1/mod_sex$summary.hyperpar$mean)),4) # grid estimates of female proportions

        male_grd[,i] = 1 - fem_grd[,i] # grid estimates of male proportions

        f.pred_dt[,i] = round(f.prop_dt[,i]*pred_dt[,i])# female counts
        f.pred_dtL[,i] = round(f.prop_dt[,i]*pred_dtL[,i])# female counts - lower
        f.pred_dtU[,i] = round(f.prop_dt[,i]*pred_dtU[,i])# female counts  upper



        m.pred_dt[,i] = pred_dt[,i] - f.pred_dt[,i] # male cunts
        m.pred_dtL[,i] = pred_dtL[,i] - f.pred_dtL[,i] # male cunts
        m.pred_dtU[,i] = pred_dtU[,i] - f.pred_dtU[,i] # male cunts

      }

   # predicted grid age-sex proportions
      prop_grd <- data.frame(prop_grd)
      prop_grdL <- data.frame(prop_grdL)
      prop_grdU <- data.frame(prop_grdU)

      fem_grd <- data.frame(fem_grd)
      fem_grdL <- data.frame(fem_grdL)
      fem_grdU <- data.frame(fem_grdU)

      male_grd <- data.frame(male_grd)
      male_grdL <- data.frame(male_grdL)
      male_grdU <- data.frame(male_grdU)

      
  ## Route 2
      # disaggregate grid pop counts based on the grid props
      for (i in admin_data$admin_id) {

        dim(grid_df <- grid_data[grid_data$admin_id == i, ])
        ids <- which(grid_data$admin_id == i)


        for (j in 1:length(age_classes)) {

          print(paste(paste0("age class ", j, " of admin ",
                             i, " is running")))

          # type 1 - sprinkle
           # grid prop predictions
          pop_grd[ids, j] <- round(prop_grd[ids, j] * grid_df$total, 5)

           # admin disaggregated
          prop_grd2[ids, j] <- prop_dt[i, j]
          pop_grd2[ids, j] <- round(prop_grd2[ids, j] * grid_df$total,5)

          fem_grdt[ids, j] <- round(fem_grd[ids, j] * pop_grd[ids, j], 5)
          male_grdt[ids, j] <- pop_grd[ids, j] - fem_grdt[ids, j]
        }
      }



      # fit model for education
      for(i in 1:length(edu_classes))
      {
        # Disaggregate Each age group count by sex
        prior.prec <- list(prec = list(prior = "pc.prec",
                                       param = c(1, 0.01))) # using PC prior
        print(paste(paste0("(",i,")"),paste0(edu_classes[i], " model is running")))

        edu_dat[,colnames(edu_dat)[i]] <- round(edu_dat[,i])

        covars <- edu_dat[,c(cov_names, "admin_id")]; dim(covars) ##---Population density

        # stack for age
        stk <- inla.stack(data=list(y=edu_dat[,i], n=edu_dat$total), #the response

                          A=list(A2,1),  #the A matrix; the 1 is included to make the list(covariates)

                          effects=list(c(list(Intercept=1), #the Intercept
                                         iset2),  #the spatial index
                                       #the covariates
                                       list(covars)
                          ),
                          #this is a quick name so you can call upon easily
                          tag='est')


        print(form_edu <- as.formula(paste0("y", " ~ ",
                                            paste(c("-1","Intercept", cov_names), collapse = " + "),
                                            " +  f(admin_id, model = 'iid', hyper = prior.prec)+
                                f(s, model = spde2)"))) # This is how you include  a typical random effect.

        mod_edu <-inla(form_edu, #the formula
                       data=inla.stack.data(stk,spde=spde2),  #the data stack
                       family= 'binomial', Ntrials = n,  #which family the data comes from
                       control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                       control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                       verbose = FALSE)

        summary(mod_edu)

        indx <-inla.stack.index(stk, "est")$data #---extract the data location indices
        prop_edu_dt[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"mean"]),6)
        prop_edu_dtL[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.025quant"]),6)
        prop_edu_dtU[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.975quant"]),6)

        pred_dt[,i] = round(prop_edu_dt[,i]*edu_dat$pop)
        pred_dtL[,i] = round(prop_edu_dtL[,i]*edu_dat$pop)
        pred_dtU[,i] = round(prop_edu_dtU[,i]*edu_dat$pop)

        summary(mod_edu)

        # Extract spatial random effect

        
    ### Route 1: Using grid cells' covariates values
        Apred <- inla.spde.make.A(mesh=mesh2, loc=coords);dim(Apred)

        # mean
        sfield_nodes <- mod_edu$summary.random$s['mean']
        field <- (Apred %*% as.data.frame(sfield_nodes)[, 1])
        # Extract IID
        IID <- rnorm(nrow(grid_data), 0, 1/mod_edu$summary.hyperpar$mean[1])
        prop_edu_grd[,i] <- round(plogis(mod_edu$summary.fixed$mean[1] +
                                           mod_edu$summary.fixed$mean[2]*grid_data$x1 +
                                           mod_edu$summary.fixed$mean[3]*grid_data$x2 +
                                           mod_edu$summary.fixed$mean[4]*grid_data$x3 +
                                           field[,1] +
                                           IID),4)

        # lower
        sfield_nodesL <- mod_edu$summary.random$s['0.025quant']
        fieldL <- (Apred %*% as.data.frame(sfield_nodesL)[, 1])

        # Extract IID
        IIDL <- rnorm(nrow(grid_data), 0, 1/mod_edu$summary.hyperpar$`0.025quant`[1])
        prop_edu_grdL[,i] <- round(plogis(mod_edu$summary.fixed$`0.025quant`[1] +
                                            mod_edu$summary.fixed$`0.025quant`[2]*grid_data$x1 +
                                            mod_edu$summary.fixed$`0.025quant`[3]*grid_data$x2 +
                                            mod_edu$summary.fixed$`0.025quant`[4]*grid_data$x3 +
                                            fieldL[,1] +
                                            IIDL),4)

        # UPPER
        sfield_nodesU <- mod_edu$summary.random$s['0.975quant']
        fieldU <- (Apred %*% as.data.frame(sfield_nodesU)[, 1])

        # Extract IID
        IIDU <- rnorm(nrow(grid_data), 0, 1/mod_edu$summary.hyperpar$`0.975quant`[1])

        prop_edu_grdU[,i] <- round(plogis(mod_edu$summary.fixed$`0.975quant`[1] +
                                            mod_edu$summary.fixed$`0.975quant`[2]*grid_data$x1 +
                                            mod_edu$summary.fixed$`0.975quant`[3]*grid_data$x2 +
                                            mod_edu$summary.fixed$`0.975quant`[4]*grid_data$x3 +
                                            fieldU[,1] +
                                            IIDU),4)

      }

      # predicted grid age-sex proportions
      prop_edu_grd <- data.frame(prop_edu_grd)
      prop_edu_grdL <- data.frame(prop_edu_grdL)
      prop_edu_grdU <- data.frame(prop_edu_grdU)

  ## Route 2
      # disaggregate grid pop counts based on the grid props
      for (i in admin_data$admin_id) {

        dim(grid_df <- grid_data[grid_data$admin_id == i, ])
        ids <- which(grid_data$admin_id == i)


        for (j in 1:length(edu_classes)) {

          print(paste(paste0("educational level class ", j, " of admin ",
                             i, " is running")))

          # type 1 - sprinkle
          # grid prop predictions
          pop_edu_grd[ids, j] <- round(prop_edu_grd[ids, j] * grid_df$total, 2)

          # admin disaggregated
          prop_edu_grd2[ids, j] <- prop_edu_dt[i, j]
          pop_edu_grd2[ids, j] <- round(prop_edu_grd2[ids, j] * grid_df$total,2)

        }
      }

      # for age
      grid_data$pred_age1 <- apply(pop_grd, 1, sum)
      grid_data$pred_age2 <- apply(pop_grd2, 1, sum)

      # for education
      grid_data$pred_edu1 <- apply(pop_edu_grd, 1, sum)
      grid_data$pred_edu2 <- apply(pop_edu_grd2, 1, sum)


      par(mfrow=c(2,2))
      plot(grid_data$total, grid_data$pred_age1)
      plot(grid_data$total, grid_data$pred_age2)
      plot(grid_data$total, grid_data$pred_edu1)
      plot(grid_data$total, grid_data$pred_edu2)


      # Save simulated counts by age groups
      grid_data$ssize <- rep(sample_size[k], nrow(grid_data))
      grid_data$age_grp <- rep(group[j], nrow(grid_data))
      write.csv(grid_data, file=paste0(result_path2,"/grid_data_with_totals.csv")) # save all predicted by age group

      # run model fit checks
      mets_age1  <- unlist(model_metrics(grid_data$total, grid_data$pred_age1)) # all ages - with covs
      mets_age2  <- unlist(model_metrics(grid_data$total, grid_data$pred_age2)) # all ages - no covs
      mets_edu1 <- unlist(model_metrics(grid_data$total, grid_data$pred_edu1)) # all educational levels combined - with covs
      mets_edu2 <- unlist(model_metrics(grid_data$total, grid_data$pred_edu2)) # all educational levels combined - no covs

      # save model fit metrics
      mets <- data.frame(rbind(mets_age1, mets_age2, mets_edu1,  mets_edu2))
      mets$ssize <- rep(sample_size[k], nrow(mets))
      mets$age_grp <- rep(group[g], nrow(mets))
      mets$data <- rownames(mets)
      print(mets)
      write.csv(mets, file=paste0(result_path2,"/fit_metrics.csv"), row.names=F) # save all predicted by age group
      #error = function(e){
      #print(e)
      #}
   }
  }

  
  
  
  
  #### plots
  # set directories
  path <- "C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts"
  data_path <- paste0(path, "/data")
  out_path <- paste0(path, "/output")
  
  
  group <- c(8, 16, 24, 32, 40)# Number of age groups
  sample_size <- c(900, 300, 100, 30, 10) # sample size
  metric <- c("MAPE", "MAE", "RMSE", "corr")
  
  
  n.age_grp <- length(group)
  n.size <- length(sample_size)
  
  
  
  dat_met1 <- list()
  dat_met2 <- list()
  for(j in 1:n.age_grp)
  {
    
    result_path1 <- paste0(out_path,"/outputs_for_", group[j],"_age_groups")
    for(k in 1:n.size)
    {
      
      result_path2 <- paste0(result_path1,"/for_", sample_size[k], "_sample_size")
      met0 <- read.csv(paste0(result_path2, "/fit_metrics.csv"))
      met0$data2 <- rep(c("age", "edu"), each=2)
      met0$age_grp <- rep(group[j], 4)
      dat_met1[[k]] = met0
    }
    dat_met2[[j]] = dat_met1
  }
  unnest_list <- unlist(dat_met2, recursive = FALSE)  #--unnest the list
  metrics <- do.call(rbind, unnest_list)
  
  library(dplyr)
  
  setwd(out_path)
  write.csv(metrics, "combined_fit_metrics.csv", row.names=FALSE)
  
  
  # Convert to long format for plotting
  # install.packages("reshape2")
  library(reshape2)
  dim(met_long <- melt(metrics, id.vars=c("data","data2","age_grp","ssize"),
                       value.name="estimate", variable.name = "metric"))
  
  #  subset for age groups
  met_age <- met_long[met_long$data2=="age",]
  met_age$Method = factor(met_age$data,
                          levels=c("mets_age1", "mets_age2"),
                          labels=c("Route 1", "Route 2"))
  table(met_age$age_grp)
  #  subset for education groups
  met_edu <- met_long[met_long$data2=="edu",]
  met_edu$Method = factor(met_edu$data,
                          levels=c("mets_edu1", "mets_edu2"),
                          labels=c("Route 1", "Route 2"))
  # use only educational level data from the full data
  met_edu <- met_edu[met_edu$age_grp==32,] # any length of age group can be used since the educational level is invariant with the length of the age groups
  met_age$age_group <- factor(met_age$age_grp)
  
  ##----- RMSE
  dim(rmse1 <- met_age[met_age$metric=="RMSE",]) # for age
  dim(rmse2 <- met_edu[met_edu$metric=="RMSE",]) # for edu
  
  
  ##----- MAE
  dim(mae1 <- met_age[met_age$metric=="MAE",]) # for age
  dim(mae2 <- met_edu[met_edu$metric=="MAE",]) # for edu
  
  
  ##----- Correlation
  dim(corr1 <- met_age[met_age$metric=="corr",]) # for age
  dim(corr2 <- met_edu[met_edu$metric=="corr",]) # for edu
  
  
  # dot plots for RMSE
  # for age groups
  library(ggpubr)
  
  
  #  - for age and sex groups
  # rmse
  rmse1$age_group <- factor(rmse1$age_grp)
  rmse1$sample_size <- factor(rmse1$ssize)
  rmse_all <- ggline(rmse1, x = "sample_size", y = "estimate",
                     error.plot = "estimate",
                     facet.by = "Method",
                     #panel.labs= list(method=c("Scaled", "Unscaled")),
                     color = "age_group",
                     palette = "aaas",
                     point.size=0.7,
                     #orientation = "horiz",
                     #color = "steelblue",
                     #label = TRUE, lab.pos = "in", lab.col = "white",
                     # position = position_dodge(),
                     #linetype = "pop_bias",
                     size=1)+
    theme_minimal()
  
  rrmse_all <- ggpar(rmse_all,
                     main = "",
                     legend = "top", legend.title=element_text("Number of \n age groups"),
                     font.legend=c(18),
                     # palette = c("aaa"),
                     font.label = list(size = 18, face = "bold", color ="red"),
                     font.x = c(18),
                     font.y = c(18),
                     font.main=c(18),
                     font.xtickslab =c(16),
                     font.ytickslab =c(18),
                     xtickslab.rt = 45, ytickslab.rt = 45,
                     
                     ggtheme = theme(strip.text.x = element_text(size = 20)),
                     xlab = "Sample size",
                     ylab = "RMSE")
  rrmse_all
  
  
  
  
  
  
  
  # mae
  mae1$age_group <- factor(mae1$age_grp)
  mae1$sample_size <- factor(mae1$ssize)
  mae_all <- ggline(mae1, x = "sample_size", y = "estimate",
                    error.plot = "estimate",
                    facet.by = "Method",
                    #panel.labs= list(method=c("Scaled", "Unscaled")),
                    color = "age_group",
                    palette = "aaas",
                    point.size=0.7,
                    #orientation = "horiz",
                    #color = "steelblue",
                    #label = TRUE, lab.pos = "in", lab.col = "white",
                    # position = position_dodge(),
                    #linetype = "pop_bias",
                    size=1)+
    theme_minimal()
  
  rmae_all <- ggpar(mae_all,
                    main = "",
                    legend = "top", legend.title=element_text("Number of \n age groups"),
                    font.legend=c(18),
                    # palette = c("aaa"),
                    font.label = list(size = 18, face = "bold", color ="red"),
                    font.x = c(18),
                    font.y = c(18),
                    font.main=c(18),
                    font.xtickslab =c(16),
                    font.ytickslab =c(18),
                    xtickslab.rt = 45, ytickslab.rt = 45,
                    
                    ggtheme = theme(strip.text.x = element_text(size = 20)),
                    xlab = "Sample size",
                    ylab = "MAE")
  rmae_all
  
  
  
  
  
  
  # corr
  corr1$age_group <- factor(corr1$age_grp)
  corr1$sample_size <- factor(corr1$ssize)
  corr_all <- ggline(corr1, x = "sample_size", y = "estimate",
                     error.plot = "estimate",
                     facet.by = "Method",
                     #panel.labs= list(method=c("Scaled", "Unscaled")),
                     color = "age_group",
                     palette = "aaas",
                     point.size=0.7,
                     #orientation = "horiz",
                     #color = "steelblue",
                     #label = TRUE, lab.pos = "in", lab.col = "white",
                     # position = position_dodge(),
                     #linetype = "pop_bias",
                     size=1)+
    theme_minimal()
  
  rcorr_all <- ggpar(corr_all,
                     main = "",
                     legend = "top", legend.title=element_text("Number of \n age groups"),
                     font.legend=c(18),
                     # palette = c("aaa"),
                     font.label = list(size = 18, face = "bold", color ="red"),
                     font.x = c(18),
                     font.y = c(18),
                     font.main=c(18),
                     font.xtickslab =c(16),
                     font.ytickslab =c(18),
                     xtickslab.rt = 45, ytickslab.rt = 45,
                     ggtheme = theme(strip.text.x = element_text(size = 20)),
                     xlab = "Sample size",
                     ylab = "CC")
  rcorr_all
  
  
  
  ggarrange(rrmse_all, rmae_all,rcorr_all,
            nrow=3, ncol=1)
  