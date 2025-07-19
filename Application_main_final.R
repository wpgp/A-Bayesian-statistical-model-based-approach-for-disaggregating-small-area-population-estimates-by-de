rm(list=ls())
#devtools::install_github("wpgp/jollofR")
#remove.packages("jollofR")
library(INLA)
#library(jollofR)
library(tidyr) # use install.packages("tidyr") to install, if not available
library(ggplot2) # use install.packages("ggplot2") to install, if not available
library(dplyr)   # use install.packages("dplyr") to install, if not available
library(terra)   # use install.packages("terra") to install, if not available
library(raster)  # use install.packages("raster") to install, if not available
library(sf)
#-------------------------------------------------
# Real data example
#-----------------------------------------------
path <- "C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper"
data_path <- paste0(path, "/data/")
output_path <- paste0(path, "/output/")

admin_data <- read.csv(paste0(data_path, "Arrondisement_data.csv")) # admin data
grid_data <- readRDS(paste0(data_path, "CMR_grid_data.rds")) # grid data
shp3 <- st_read(paste0(data_path, "arrondissement_shp.shp")) # read in the subdivisional shapefile

# st_crs(shp3) # check the projection
# Reproject to WGS84 (longitude/latitude)
shp3_lnlt <- st_transform(shp3, crs = 4326) # EPSG:4326 is WGS84

                                                                                           
plot(shp3["pop2005"])
names(grid_data)

grid_data$grd_id <- grid_data$grd
names(admin_data)
table(admin_data$set_typ <- factor(admin_data$set_typ)) # settlement type


# extract class names for age
age_dat <-  admin_data %>% dplyr::select(starts_with("age_"))



# define output matrices for age disaggregated counts and proportions
pred_dt <- prop_dt <- age_dat # mean
pred_dtL <- prop_dtL <- age_dat # lower bound
pred_dtU <- prop_dtU <- age_dat # upper bound

# obtain rowsum of age group classes
age_classes <-  names(age_dat)
age_dat$total <- apply(age_dat[,age_classes], 1, sum, na.rm=T)# rowsum of age group counts
age_dat$total[age_dat$total==0] = NA # set all zero total counts to NA


# extract geospatial covariates from the admin data
cov_names <- c("x3", "x17", "x21", "x32", "x34", "x40", "x42") # best fit covariates
covs2 <- admin_data %>% dplyr::select(cov_names)

# standardize the selected covariates to make them comparable
stdize <- function(x)
{
  stdz <- (x - mean(x, na.rm=T))/sd(x, na.rm=T)
  return(stdz)
}

covs <- data.frame(apply(covs2, 2, stdize))
age_dat$admin_id <- admin_data$admin_id
age_dat$pop <- admin_data$total # add the total count to be disaggregated
age_dat$set_typ <- factor(admin_data$set_typ) # add settlement type
age_dat <- cbind(age_dat,covs) # join the age data and the geospatial covariates

# extract female age data, the classes and the output matrices
f_dat = admin_data %>% dplyr::select(starts_with("fage_"))
f_age_classes <-  names(f_dat)
 # output matrices for female age-sex disaggregated data
f.pred_dt <- f.prop_dt <-m.pred_dt <- m.prop_dt <- f_dat
f.pred_dtL <- f.prop_dtL <-m.pred_dtL <- m.prop_dtL <- f_dat
f.pred_dtU <- f.prop_dtU <-m.pred_dtU <- m.prop_dtU <- f_dat

f_dat$admin_id <- admin_data$admin_id
# extract educational group data
edu_dat <-  admin_data %>% dplyr::select(starts_with("edu_"))

  # output matrices for educational level disaggregated data
pred_edu_dt <- prop_edu_dt <- edu_dat
pred_edu_dtL <- prop_edu_dtL <- edu_dat
pred_edu_dtU <- prop_edu_dtU <- edu_dat
edu_classes <-  names(edu_dat)

edu_dat$total <- apply(edu_dat[,edu_classes], 1, sum, na.rm=T) # row sum across educational levels
edu_dat$total[edu_dat$total==0] = NA # set zero total counts to NA
edu_dat$admin_id <- admin_data$admin_id
edu_dat$pop <- admin_data$total # Add total counts to be disaggregated
edu_dat <- cbind(edu_dat,covs) # join educational level data to the geospatial covariates


# Create grid cell level (100m by 100m) outputs matrices
 # age proportions
male_grd <- fem_grd <- prop_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(age_classes))
male_grdL <- fem_grdL <- prop_grdL <- male_grd
male_grdU<- fem_grdU <- prop_grdU <- male_grd
prop_grd2 <- prop_grd3 <- prop_gr4   <- prop_grd

  # age counts
male_grdt <- fem_grdt <- pop_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(age_classes))
male_grdtL <- fem_grdtL <- pop_grdL <- male_grd
male_grdtU<- fem_grdtU <- pop_grdU <- male_grd
pop_grd2 <- pop_grd3 <- pop_gr4   <- pop_grd

# For education
prop_edu_grd <- pop_edu_grd <- matrix(0, nrow=nrow(grid_data), ncol=length(edu_classes))
prop_edu_grdL <- prop_edu_grdU <- prop_edu_grd
pop_edu_grdL <- pop_edu_grdU <- pop_edu_grd
prop_edu_grd2 <- prop__edu_grd3 <- prop_edu_gr4   <- prop_edu_grd
pop_edu_grd2 <- pop__edu_grd3 <- pop_edu_gr4   <- pop_edu_grd


# Build mesh. SPDE object and projection matrix for both estimation and prediction

  # admin - level estimation
plot(admin_data$lon, admin_data$lat)

# Build the mesh, Projection matix and spde at admin level
coords <- cbind(admin_data$lon, admin_data$lat)
plot(coords)

# non-convex hull mesh
non_convex_bdry <- inla.nonconvex.hull(coords, -0.03, -0.05, resolution = c(100, 100))
plot(non_convex_bdry$loc)

mesh <- inla.mesh.2d(boundary = non_convex_bdry, max.edge=c(0.05,0.7),
                      offset = c(0.1, 0.5),
                      cutoff = 0.1)
par(mfrow=c(1,1))
plot(mesh)
plot(coords, add=T, lwd=1.5, col=2, pch=16)
plot(mesh, add=T)

mesh$n # 3214 mesh nodes


#  Create the SPDE
spde <- inla.spde2.matern(mesh, alpha=2)

#  Build projector matrix A
A <- inla.spde.make.A(mesh=mesh, loc=coords);dim(A)

#  specify the observation indices for estimation
iset <- inla.spde.make.index(name = "s", spde$n.spde)

## fit age models
for(i in 1:length(age_classes))
{
  # i=2
  # Disaggregate Each age group count by sex
  prior.prec <- list(prec = list(prior = "pc.prec",
                                 param = c(1, 0.01))) # using PC prior
  print(paste(paste0("(",i,")"),paste0(age_classes[i], " model is running")))

  age_dat[,colnames(age_dat)[i]] <- round(age_dat[,i])

  covars <- age_dat[,c(cov_names,"set_typ", "admin_id")]; dim(covars) ##---Population density

  # stack for age
  stk <- inla.stack(data=list(y=age_dat[,i], n=age_dat$total), #the response

                    A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)

                    effects=list(c(list(Intercept=1), #the Intercept
                                   iset),  #the spatial index
                                 #the covariates
                                 list(covars)
                    ),
                    #this is a quick name so you can call upon easily
                    tag='est')


  print(form_age <- as.formula(paste0("y", " ~ ",
                                      paste(c("-1","Intercept", cov_names), collapse = " + "),
                                      " +   f(admin_id, model = 'iid', hyper = prior.prec)+
                                      f(set_typ, model = 'iid', hyper = prior.prec) +  f(s, model = spde)"))) # This is how you include  a typical random effect.
  mod_age <-inla(form_age, #the formula
                 data=inla.stack.data(stk,spde=spde),  #the data stack
                 family= 'binomial', Ntrials = n,  #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE)

  summary(mod_age)

  
  # Save fixed and random effects estimates for each group
  parameter_dir <- paste0(output_path, "fixed_and_random_effects_age/",age_classes[i])
  if (!dir.exists(parameter_dir)) {
    dir.create(parameter_dir, recursive = TRUE)
    message(paste("Directory", parameter_dir, "created successfully."))
  }else{
    message(paste("Directory", parameter_dir, "already exists."))
  }
  capture.output(summary(mod_age), file = paste0(parameter_dir, "/posterior_estimates.txt"))



  indx <-inla.stack.index(stk, "est")$data #---extract the data location indices
  prop_dt[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"mean"]),6)
  prop_dtL[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.025quant"]),6)
  prop_dtU[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.975quant"]),6)

  pred_dt[,i] = round(prop_dt[,i]*age_dat$pop)
  pred_dtL[,i] = round(prop_dtL[,i]*age_dat$pop)
  pred_dtU[,i] = round(prop_dtU[,i]*age_dat$pop)

  f_dat[,colnames(f_dat)[i]] <- round(f_dat[,i]) # ensure counts are integers


# sex disaggregation
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


  f.pred_dt[,i] = round(f.prop_dt[,i]*pred_dt[,i])# female counts
  f.pred_dtL[,i] = round(f.prop_dt[,i]*pred_dtL[,i])# female counts - lower
  f.pred_dtU[,i] = round(f.prop_dt[,i]*pred_dtU[,i])# female counts  upper

  m.pred_dt[,i] = pred_dt[,i] - f.pred_dt[,i] # male cunts
  m.pred_dtL[,i] = pred_dtL[,i] - f.pred_dtL[,i] # male cunts
  m.pred_dtU[,i] = pred_dtU[,i] - f.pred_dtU[,i] # male cunts

}



# disaggregate grid pop counts based on the grid in admin props
for (i in admin_data$admin_id) {

  dim(grid_df <- grid_data[grid_data$admin_id == i, ])
  ids <- which(grid_data$admin_id == i)


  for (j in 1:length(age_classes)) {

    print(paste(paste0("age class ", j, " of admin ",
                       i, " is running")))

    # type 1 - sprinkle
    # admin disaggregated
    prop_grd2[ids, j] <- prop_dt[i, j]
    pop_grd2[ids, j] <- round(prop_grd2[ids, j] * grid_df$total,5)


    fem_grd[ids, j] <- f.prop_dt[i,j]
    fem_grdt[ids, j] <- round(fem_grd[ids, j] * pop_grd2[ids, j], 5)
    male_grdt[ids, j] <- pop_grd2[ids, j] - fem_grdt[ids, j]


    # Lower
    prop_grdL[ids, j] <- prop_dtL[i, j]
    pop_grdL[ids, j] <- round(prop_grdL[ids, j] * grid_df$total,2)

    #Upper
    prop_grdU[ids, j] <- prop_dtU[i, j]
    pop_grdU[ids, j] <- round(prop_grdU[ids, j] * grid_df$total,2)
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

                    A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)

                    effects=list(c(list(Intercept=1), #the Intercept
                                   iset),  #the spatial index
                                 #the covariates
                                 list(covars)
                    ),
                    #this is a quick name so you can call upon easily
                    tag='est')


  print(form_edu <- as.formula(paste0("y", " ~ ",
                                      paste(c("-1","Intercept", cov_names), collapse = " + "),
                                      " +  f(admin_id, model = 'iid', hyper = prior.prec)+
                                f(s, model = spde)"))) # This is how you include  a typical random effect.

  mod_edu <-inla(form_edu, #the formula
                 data=inla.stack.data(stk,spde=spde),  #the data stack
                 family= 'binomial', Ntrials = n,  #which family the data comes from
                 control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                 control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                 verbose = FALSE)

  summary(mod_edu)

  # Save fixed and random effects estimates for each group
  parameter_dir <- paste0(output_path, "fixed_and_random_effects_edu/",edu_classes[i])
  if (!dir.exists(parameter_dir)) {
    dir.create(parameter_dir, recursive = TRUE)
    message(paste("Directory", parameter_dir, "created successfully."))
  }
  else {
    message(paste("Directory", parameter_dir, "already exists."))
  }
  capture.output(summary(mod_edu), file = paste0(parameter_dir, "/posterior_estimates.txt"))

  prop_edu_dt[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"mean"]),6)
  prop_edu_dtL[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.025quant"]),6)
  prop_edu_dtU[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.975quant"]),6)

  pred_edu_dt[,i] = round(prop_edu_dt[,i]*edu_dat$pop)
  pred_edu_dtL[,i] = round(prop_edu_dtL[,i]*edu_dat$pop)
  pred_edu_dtU[,i] = round(prop_edu_dtU[,i]*edu_dat$pop)

  summary(mod_edu)

}





## route 2
# disaggregate grid pop counts based on the grid props
for (i in admin_data$admin_id) {

  dim(grid_df <- grid_data[grid_data$admin_id == i, ])
  ids <- which(grid_data$admin_id == i)


  for (j in 1:length(edu_classes)) {

    print(paste(paste0("educational level class ", j, " of admin ",
                       i, " is running")))
    # admin disaggregated

    # Mean
    prop_edu_grd2[ids, j] <- prop_edu_dt[i, j]
    pop_edu_grd2[ids, j] <- round(prop_edu_grd2[ids, j] * grid_df$total,2)

    # Lower
    prop_edu_grdL[ids, j] <- prop_edu_dtL[i, j]
    pop_edu_grdL[ids, j] <- round(prop_edu_grdL[ids, j] * grid_df$total,2)

    #Upper
    prop_edu_grdU[ids, j] <- prop_edu_dtU[i, j]
    pop_edu_grdU[ids, j] <- round(prop_edu_grdU[ids, j] * grid_df$total,2)
  }
}


#### Extract grid cells posterior estimates
# for age
#pop_grd2 <- apply(pop_grd2, 2, round)
pop_grd2 <- data.frame(pop_grd2)
names(pop_grd2) <- age_classes
grid_data2 <- cbind(grid_data, pop_grd2)
grid_data2$pred_age_tot <- apply(pop_grd2, 1, sum)
grid_data2$pred_age_totL <- apply(pop_grd2, 1, quantile, probs=c(0.025), na.rm=T)
grid_data2$pred_age_totU <- apply(pop_grd2, 1, quantile, probs=c(0.975), na.rm=T)


# for education
#pop_edu_grd2 <- apply(pop_edu_grd2, 2, round)
pop_edu_grd2 <- data.frame(pop_edu_grd2)
colnames(pop_edu_grd2) <- edu_classes
grid_data2 <- cbind(grid_data2, pop_edu_grd2)
grid_data2$pred_edu_tot <- apply(pop_edu_grd2, 1, sum)
grid_data2$pred_edu_totL <- apply(pop_edu_grd2, 1, quantile, probs=c(0.025), na.rm=T)
grid_data2$pred_edu_totU <- apply(pop_edu_grd2, 1, quantile, probs=c(0.975), na.rm=T)

# for age - female
fem_grdt <- data.frame(fem_grdt)
names(fem_grdt) <- f_age_classes
fem_grdt <- round(fem_grdt,2)
grid_data2 <- cbind(grid_data2, fem_grdt)


# for age - male
male_grdt <- data.frame(male_grdt)
names(male_grdt) <- paste0("mage_",age_classes)
male_grdt <- round(male_grdt,2)
grid_data2 <- cbind(grid_data2, male_grdt)


#### Extract admin levels posterior estimates
 # for age-sex
names(admin_data)
admin_data2 <- cbind(admin_data[,c("admin_id", "total")], pred_dt)
admin_data2$predicted <- apply(pred_dt, 1, sum)
admin_data2$predictedL <- apply(pred_dtL, 1, sum)
admin_data2$predictedU <- apply(pred_dtU, 1, sum)


# for age-sex
names(admin_data)
admin_data3 <- cbind(admin_data[,c("admin_id", "total")], pred_edu_dt)
admin_data3$predicted_edu <- apply(pred_edu_dt, 1, sum)
admin_data3$predicted_eduL <- apply(pred_edu_dtL, 1, sum)
admin_data3$predicted_eduU <- apply(pred_edu_dtU, 1, sum)


saveRDS(grid_data2,paste0(output_path, "CMR_post_grid_data.rds"))
# Visualise scatter plots
par(mfrow=c(1,2))
plot(grid_data2$total, grid_data2$pred_age_tot)
plot(grid_data2$total, grid_data2$pred_edu_tot)

# run model fit checks
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
mets_age2  <- unlist(model_metrics(grid_data2$total, grid_data2$pred_age_tot)) # all ages - no covs
mets_edu2 <- unlist(model_metrics(grid_data2$total, grid_data2$pred_edu_tot)) # all educational levels combined - no covs

# save model fit metrics
mets <- data.frame(rbind(mets_age2, mets_edu2))
print(mets)
write.csv(mets, file=paste0(output_path,"/fit_metrics.csv"), row.names=F) # save all predicted by age group


# ----Make pyramid plot----------------

  ### female data
f_mat <- f.pred_dt
age_classes <- names(f_mat)
f_mat$id <- 1:nrow(f_mat)
(female_pop <- reshape2::melt(f_mat, id=c("id"), value.name="Population",
                              variable.name="Age"))
female_pop$Age <- factor(female_pop$Age)
levels(female_pop$Age) <- gsub("pp_fage_", "", age_classes)
levels(female_pop$Age) <- gsub("fage_", "", levels(female_pop$Age))
levels(female_pop$Age) <- gsub("_", "-", levels(female_pop$Age))
levels(female_pop$Age)[18] <- "85+"
female_pop$Gender <- rep("Female", nrow(female_pop))# create gender variable


### male data
# mean
m_mat <- m.pred_dt #
age_classes <- names(m_mat)
m_mat$id <- 1:nrow(m_mat)
(male_pop <- reshape2::melt(m_mat, id=c("id"), value.name="Population",
                            variable.name="Age"))
male_pop$Age <- factor(male_pop$Age)
levels(male_pop$Age) <- gsub("pp_fage_", "", age_classes)
levels(male_pop$Age) <- gsub("fage_", "", levels(male_pop$Age))
levels(male_pop$Age) <- gsub("_", "-", levels(male_pop$Age))
levels(male_pop$Age)[18] <- "85+"
male_pop$Gender <- rep("Male", nrow(male_pop))# create gender variable


# combine both datasets
dim(pop_pyramid <- rbind(female_pop, male_pop)) # mean

pop_pyramid$Population2 <- pop_pyramid$Population/1000
# Create the age-sex population pyramid
# national - mean
population_pyramid1  <-  ggplot(
  pop_pyramid,
  aes(
    x = Age,
    fill = Gender,
    y = ifelse(
      test = Gender == "Male",
      yes = -Population2,
      no = Population2
    )
  )
) +
  geom_bar(stat = "identity") +
  scale_y_continuous(
    labels = abs
  ) +
  labs(
    x = "Age",
    y = "Population",
    fill = "Gender"#,
  ) +
  theme_minimal() +
  coord_flip()

pyramid1 <- ggpubr::ggpar(population_pyramid1, ylab="Population Count ('000)", xlab="Age (years)",
                          legend = "right", legend.title = "Gender",size=22,
                          font.legend=c(16),
                          palette = "lancet",
                          font.label = list(size = 15, face = "bold", color ="red"),
                          font.x = c(16),
                          font.y = c(16),
                          font.main=c(14),
                          font.xtickslab =c(14),
                          font.ytickslab =c(16)#,
                          # orientation = "reverse",
                          #xtickslab.rt = 45, ytickslab.rt = 45
)
print(pyramid1)


# Scatter plots with 95% credible interval lines

# age -sex
names(admin_data2)
varsRequired2 <- c("total","predicted", "predictedL", "predictedU")
age_pred <- admin_data2[, varsRequired2] # extract key variables for age -sex grid data
names(age_pred) <- c("Total","Predicted", "Lower", "Upper")
age_pred$data <- rep("AgeSex", nrow(age_pred))
age_pred$Total <- round(age_pred$Total)
age_pred$age_0_4 <- admin_data$age_0_4


age_pred[order(age_pred$Total),]
# educational level 
varsRequired3 <- c("total","predicted_edu", "predicted_eduL", "predicted_eduU")
edu_pred <- admin_data3[, varsRequired3] # extract key variables for educational level grid data
names(edu_pred) <- c("Total","Predicted", "Lower", "Upper")
edu_pred$data <- rep("Education", nrow(edu_pred))
edu_pred$Total <- round(edu_pred$Total)
edu_pred$edu_0_4 <- admin_data$edu_0_4
dat_lng <- rbind(age_pred, edu_pred)


dat_lng$data <- factor(dat_lng$data)


# Create scatter plot with vertical confidence intervals
scat <- ggplot(dat_lng, aes(x = Total,y = Predicted)) +
  geom_point(size = 3, color = "black") +  # Scatter points
  geom_line(data = dat_lng, aes(x = Total,y = Predicted), color = "darkblue",size=1.5) +
  geom_errorbar(aes(ymin = Lower, ymax = Upper), width =2,size=0.6, color = "red") +  # CI lines
  #geom_pointrange(aes(ymin = Lower, ymax = Upper), size=0.6, color = "black") +  # CI lines
  theme_minimal()+
facet_wrap(~data)

pscat <- ggpubr::ggpar(scat , xlab="Total input population", ylab="Predicted total population",
                           legend = "none", legend.title = "",size=22,
                           font.legend=c(18),
                           palette = "lancet",
                           font.label = list(size = 18, face = "bold", color ="red"),
                           font.x = c(18),
                           font.y = c(18),
                           font.main=c(18),
                           font.xtickslab =c(18),
                           font.ytickslab =c(18),
                           xtickslab.rt = 45
)
pscat

#---- box and line plots ------------------------------------
  # combined
classes <- names(pred_dt)
dmat_lng <- tidyr::gather(as.data.frame(pred_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("age_", "", classes)
levels(dmat_lng$key) <- gsub("_", "-", levels(dmat_lng$key))
levels(dmat_lng$key)[18] <- "85+"

dmat_lng$value2 <- dmat_lng$value/1000
plt_box1 <- ggplot(dmat_lng, aes(y=value2, x=key, col = key)) +
  geom_boxplot(fill="grey")+
  theme_minimal()


page_box1 <- ggpubr::ggpar(plt_box1 , xlab="Age group (years)", ylab="Population count ('000)",
                          legend = "none", legend.title = "",size=22,
                          font.legend=c(18),
                          palette = "",
                          font.label = list(size = 14, face = "bold", color ="red"),
                          font.x = c(16),
                          font.y = c(14),
                          font.main=c(14),
                          font.xtickslab =c(14),
                          font.ytickslab =c(14),
                          xtickslab.rt = 45
)
page_box1

# females
classes <- names(f.pred_dt)
dmat_lng <- tidyr::gather(as.data.frame(f.pred_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("fage_", "", classes)
levels(dmat_lng$key) <- gsub("_", "-", levels(dmat_lng$key))
levels(dmat_lng$key)[18] <- "85+"

dmat_lng$value2 <- dmat_lng$value/1000
plt_box2 <- ggplot(dmat_lng, aes(y=value2, x=key, col = key)) +
  geom_boxplot(fill="grey")+
  theme_minimal()


page_box2 <- ggpubr::ggpar(plt_box2 , xlab="Age group (years)", ylab="Population count ('000)",
                           legend = "none", legend.title = "",size=22,
                           font.legend=c(18),
                           palette = "",
                           font.label = list(size = 18, face = "bold", color ="red"),
                           font.x = c(18),
                           font.y = c(18),
                           font.main=c(18),
                           font.xtickslab =c(16),
                           font.ytickslab =c(16),
                           xtickslab.rt = 90
)
page_box2

# males
classes <- names(m.pred_dt)
dmat_lng <- tidyr::gather(as.data.frame(m.pred_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("fage_", "", classes)
levels(dmat_lng$key) <- gsub("_", "-", levels(dmat_lng$key))
levels(dmat_lng$key)[18] <- "85+"

dmat_lng$value2 <- dmat_lng$value/1000
plt_box3 <- ggplot(dmat_lng, aes(y=value2, x=key, col = key)) +
  geom_boxplot(fill="grey")+
  theme_minimal()


page_box3 <- ggpubr::ggpar(plt_box3 , xlab="Age group (years)", ylab="Population count ('000)",
                           legend = "none", legend.title = "",size=22,
                           font.legend=c(18),
                           palette = "",
                           font.label = list(size = 18, face = "bold", color ="red"),
                           font.x = c(18),
                           font.y = c(18),
                           font.main=c(18),
                           font.xtickslab =c(16),
                           font.ytickslab =c(16),
                           xtickslab.rt = 90
)
page_box3

# Make line graph
 # all

classes <- names(pred_dt)
total1 <- apply(pred_dt, 2, sum)
tot_dat1 <- data.frame(Population = total1,
                      Key = 1:length(classes),
                      label = levels(dmat_lng$key))

tot_dat1$Population2 <- tot_dat1$Population/1000
plt_line1 <- ggplot(tot_dat1, aes(y=Population2, x=Key)) +
  geom_point(size=3)+
  geom_line(size=1.5, color = "darkblue")+
  scale_x_continuous(breaks = seq(1, length(classes),
                                  by =1), labels =levels(dmat_lng$key))+
  theme_minimal()
page_line1 <- ggpubr::ggpar(plt_line1 , xlab="Age group (years)", ylab="Predicted total population ('000)",
                           legend = "none", legend.title = "",size=22,
                           font.legend=c(18),
                           palette = "pnj",
                           font.label = list(size = 18, face = "bold", color ="red"),
                           font.x = c(18),
                           font.y = c(18),
                           font.main=c(18),
                           font.xtickslab =c(16),
                           font.ytickslab =c(16),
                           xtickslab.rt = 90
)
page_line1


# females
total2 <- apply(f.pred_dt, 2, sum)
tot_dat2 <- data.frame(Population = total2,
                       Key = 1:length(classes),
                       label = levels(dmat_lng$key))

tot_dat2$Population2 <- tot_dat2$Population/1000
plt_line2 <- ggplot(tot_dat2, aes(y=Population2, x=Key)) +
  geom_line(size=1)+
  scale_x_continuous(breaks = seq(1, length(classes),
                                  by =1), labels =levels(dmat_lng$key))+
  theme_bw()
page_line2 <- ggpubr::ggpar(plt_line2 , xlab="Age group (years)", ylab="Predicted total population ('000)",
                            legend = "none", legend.title = "",size=22,
                            font.legend=c(18),
                            palette = "pnj",
                            font.label = list(size = 16, face = "bold", color ="red"),
                            font.x = c(16),
                            font.y = c(16),
                            font.main=c(16),
                            font.xtickslab =c(16),
                            font.ytickslab =c(16),
                            xtickslab.rt = 90
)
page_line2



# males
total3 <- apply(m.pred_dt, 2, sum)
tot_dat3 <- data.frame(Population = total3,
                       Key = 1:length(classes),
                       label = levels(dmat_lng$key))

tot_dat3$Population2 <- tot_dat3$Population/1000
plt_line3 <- ggplot(tot_dat3, aes(y=Population2, x=Key)) +
  geom_line(size=1)+
  scale_x_continuous(breaks = seq(1, length(classes),
                                  by =1), labels =levels(dmat_lng$key))+
  theme_bw()
page_line3 <- ggpubr::ggpar(plt_line3 , xlab="Age group (years)", ylab="Predicted total population ('000)",
                            legend = "none", legend.title = "",size=22,
                            font.legend=c(16),
                            palette = "pnj",
                            font.label = list(size = 16, face = "bold", color ="red"),
                            font.x = c(16),
                            font.y = c(16),
                            font.main=c(16),
                            font.xtickslab =c(16),
                            font.ytickslab =c(16),
                            xtickslab.rt = 90
)
page_line3


print(ggpubr::ggarrange(#page_box1,page_line1,
page_line2,page_line3,
  #labels = c("(a)","(b)",
  #           "(c)","(d)",
  #            "(e)","(f)"),
  nrow=1, ncol=2))



print(ggpubr::ggarrange(#page_box1,page_line1,
                        page_box2, page_line2,
                        page_box3, page_line3,
                        #labels = c("(a)","(b)",
                        #           "(c)","(d)",
                       #            "(e)","(f)"),
                        nrow=2, ncol=2))




#----histograms------------------------------------------------

# females
classes <- names(f.pred_dt)
dmat_lng <- tidyr::gather(as.data.frame(f.pred_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("fage_", "", classes)
levels(dmat_lng$key) <- gsub("_", "-", levels(dmat_lng$key))
levels(dmat_lng$key)[18] = "85+"

dens_dist2 <- ggplot(dmat_lng, aes(x=value, color=key)) +
  geom_histogram(alpha=0.3,size=0.8)+
  theme_bw()+
  theme(strip.text = element_text(size=18))+
  facet_wrap(~key)
pdens_dist2 <- ggpubr::ggpar(dens_dist2, ylab="Frequency", xlab="Population Count",
                            legend = "none", legend.title = "",size=18,
                            font.legend=c(18),
                            palette = "pnj",
                            font.label = list(size = 15, face = "bold", color ="red"),
                            font.x = c(18),
                            font.y = c(18),
                            font.main=c(18),
                            font.xtickslab =c(18),
                            font.ytickslab =c(18),
                            xtickslab.rt = 45
)
print(pdens_dist2)



# males
classes <- names(m.pred_dt)
dmat_lng <- tidyr::gather(as.data.frame(m.pred_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("fage_", "", classes)
levels(dmat_lng$key) <- gsub("_", "-", levels(dmat_lng$key))
levels(dmat_lng$key)[18] = "85+"

dens_dist3 <- ggplot(dmat_lng, aes(x=value, color=key)) +
  geom_histogram(alpha=0.3,size=0.8)+
  theme_bw()+
  theme(strip.text = element_text(size=18))+
  facet_wrap(~key)
pdens_dist3 <- ggpubr::ggpar(dens_dist3, ylab="Frequency", xlab="Population Count",
                             legend = "none", legend.title = "",size=18,
                             font.legend=c(18),
                             palette = "pnj",
                             font.label = list(size = 15, face = "bold", color ="red"),
                             font.x = c(18),
                             font.y = c(18),
                             font.main=c(18),
                             font.xtickslab =c(18),
                             font.ytickslab =c(18),
                             xtickslab.rt = 45
)
print(pdens_dist3)



# ---------------------------------------------------------
# Make box and line plots and histogram for educational groups
#-------------------------------------------------------------
#---- box and line plots ------------------------------------
# combined
classes <- names(pred_edu_dt)
dmat_lng <- tidyr::gather(as.data.frame(pred_edu_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("edu_", "", classes)

dmat_lng$value2 <- dmat_lng$value/1000
plt_box1b <- ggplot(dmat_lng, aes(y=value2, x=key, col = key)) +
  geom_boxplot(fill="grey")+
  theme_minimal()


page_box1b <- ggpubr::ggpar(plt_box1b , xlab="Educational level", ylab="Predicted total population ('000)",
                           legend = "none", legend.title = "",size=22,
                           font.legend=c(18),
                           palette = "",
                           font.label = list(size = 18, face = "bold", color ="red"),
                           font.x = c(18),
                           font.y = c(16),
                           font.main=c(18),
                           font.xtickslab =c(18),
                           font.ytickslab =c(18),
                           xtickslab.rt = 45
)
page_box1b



# Make line graph
# all
total1b <- apply(pred_edu_dt, 2, sum, na.rm=T)
tot_dat1b <- data.frame(Population = total1b,
                       Key = 1:length(classes),
                       label = levels(dmat_lng$key))

tot_dat1b$Population2 <- tot_dat1b$Population/1000
plt_line1b <- ggplot(tot_dat1b, aes(y=Population2, x=Key)) +
  geom_point(size=3)+
  geom_line(size=1.5, color = "darkblue")+
  scale_x_continuous(breaks = seq(1, length(classes),
                                  by =1), labels =levels(dmat_lng$key))+
  theme_minimal()
page_line1b <- ggpubr::ggpar(plt_line1b , xlab="Educational level", ylab="Predicted total population ('000)",
                            legend = "none", legend.title = "",size=22,
                            font.legend=c(18),
                            palette = "lancet",
                            font.label = list(size = 18, face = "bold", color ="red"),
                            font.x = c(18),
                            font.y = c(16),
                            font.main=c(18),
                            font.xtickslab =c(18),
                            font.ytickslab =c(18),
                            xtickslab.rt = 45
)
page_line1b


print(ggpubr::ggarrange(#page_box1,page_line1,
  page_line1, page_line1b,
  nrow=1, ncol=2))


print(ggpubr::ggarrange(#page_box1,page_line1,
  page_box1b, page_line1b,
  nrow=1, ncol=2))

#----histograms for educational levels
classes <- names(pred_edu_dt)
dmat_lng <- tidyr::gather(as.data.frame(pred_edu_dt))
dmat_lng$key <- factor(dmat_lng$key,
                       levels=classes,
                       labels=classes)
levels(dmat_lng$key) <- gsub("edu_", "", classes)

dmat_lng$value2 <- dmat_lng$value/1000
dens_dist2b <- ggplot(dmat_lng, aes(x=value2, color=key)) +
  geom_histogram(alpha=0.3,size=0.8)+
  theme_bw()+
  theme(strip.text = element_text(size=18))+
  facet_wrap(~key)
pdens_dist2b <- ggpubr::ggpar(dens_dist2b, ylab="Frequency", xlab="Predicted total population ('000)",
                             legend = "none", legend.title = "",size=18,
                             font.legend=c(18),
                             palette = "lancet",
                             font.label = list(size = 15, face = "bold", color ="red"),
                             font.x = c(18),
                             font.y = c(18),
                             font.main=c(18),
                             font.xtickslab =c(18),
                             font.ytickslab =c(18),
                             xtickslab.rt = 45
)
print(pdens_dist2b)



#-------------------------------------------------------------------------------------
#---raster plots -------------------------------------------
#------------------------------------------------------------------------------------

# Write the raster files
   # age and sex rasters


# specify the ereference coordinates
ref_coords <- cbind(grid_data2$lon, grid_data2$lat)
xx <- as.matrix(ref_coords)
# specify the names of the raster files to save
rclass <- paste0("CMR_population_v1_0_",paste0("age_",
                                               paste0(seq(0,85, by=5), "_",c(seq(4,84, by=5), "above"))))
rclassL <- paste0(rclass, "L") # lower bound
rclassU <- paste0(rclass, "U") # upper bound
frclass <- gsub("_age", "_agesex_f", rclass) # female age classes
mrclass <- gsub("_age", "_agesex_m", rclass) # male age classes

output_dir <- output_path
#pop_edu_grd2
pop_grd2
for (k in 1:length(rclass)) {

   print(paste0(output_dir, "/pop_", rclass[k], ".tif"))
    z1b <- as.matrix(pop_grd2[, k])
    h1b <- rasterFromXYZ(cbind(xx, z1b))
    writeRaster(h1b, filename = paste0(output_dir, "/pop_",
                                       rclass[k], ".tif"), overwrite = TRUE, options = c(COMPRESS = "LZW"))


    print(paste0(output_dir, "/pop_", rclass[k], "_lower",
                 ".tif"))
    z2b <- as.matrix(pop_grdL[, k])
    h2b <- rasterFromXYZ(cbind(xx, z2b))
    writeRaster(h2b, filename = paste0(output_dir, "/pop_",
                                       rclass[k], "_lower", ".tif"), overwrite = TRUE,
                options = c(COMPRESS = "LZW"))


    print(paste0(output_dir, "/pop_", rclass[k], "_upper",
                 ".tif"))
    z3b <- as.matrix(pop_grdU[, k])
    h3b <- rasterFromXYZ(cbind(xx, z3b))
    writeRaster(h3b, filename = paste0(output_dir, "/pop_",
                                       rclass[k], "_upper", ".tif"), overwrite = TRUE,
                options = c(COMPRESS = "LZW"))


    print(paste0(output_dir, "/pop_", frclass[k], ".tif"))
    z4b <- as.matrix(fem_grdt[, k])
    h4b <- rasterFromXYZ(cbind(xx, z4b))
    writeRaster(h4b, filename = paste0(output_dir, "/pop_",
                                       frclass[k], ".tif"), overwrite = TRUE, options = c(COMPRESS = "LZW"))


    print(paste0(output_dir, "/pop_", mrclass[k], ".tif"))
    z5b <- as.matrix(male_grdt[, k])
    h5b <- rasterFromXYZ(cbind(xx, z5b))
    writeRaster(h5b, filename = paste0(output_dir, "/pop_",
                                       mrclass[k], ".tif"), overwrite = TRUE, options = c(COMPRESS = "LZW"))
    print(paste0(output_dir, "/pop_", mrclass[k], "_lower",
                 ".tif"))
  }


# Visualise rasters

# AGE 0-4 years - combined
ras1<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output//pop_CMR_population_v1_0_age_0_4.tif")
ras1L<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output//pop_CMR_population_v1_0_age_0_4_lower.tif")
ras1U<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output//pop_CMR_population_v1_0_age_0_4_upper.tif")

par(mfrow=c(3,3), mar=c(3,3,2,1))
plot(ras1, cex.lab=1.5,
     cex.axis=2) # visulize raster
plot(ras1L)
plot(ras1U)


#---------------------------------------------------------------------------------
## WRITE RASTERS FOR EDUCATION LEVELS
# specify the ereference coordinates
ref_coords <- cbind(grid_data2$lon, grid_data2$lat)
xx <- as.matrix(ref_coords)
# specify the names of the raster files to save
rclass <- paste0("CMR_population_v1_0_",paste0("edu_",c("None", "Primary", "Secondary", "Higher", "Unknown")))
rclassL <- paste0(rclass, "L") # lower bound
rclassU <- paste0(rclass, "U") # upper bound


output_dir <- paste0(output_path, "Education")

for (k in 1:length(rclass)) {
  # k= 2
  print(paste0(output_dir, "/pop_", rclass[k], ".tif"))
  z1b <- as.matrix(pop_edu_grd2[, k])
  h1b <- rasterFromXYZ(cbind(xx, z1b))
  writeRaster(h1b, filename = paste0(output_dir, "/pop_",
                                     rclass[k], ".tif"), overwrite = TRUE, options = c(COMPRESS = "LZW"))


  print(paste0(output_dir, "/pop_", rclass[k], "_lower",
               ".tif"))
  z2b <- as.matrix(pop_edu_grdL[, k])
  h2b <- rasterFromXYZ(cbind(xx, z2b))
  writeRaster(h2b, filename = paste0(output_dir, "/pop_",
                                     rclass[k], "_lower", ".tif"), overwrite = TRUE,
              options = c(COMPRESS = "LZW"))


  print(paste0(output_dir, "/pop_", rclass[k], "_upper",
               ".tif"))
  z3b <- as.matrix(pop_edu_grdU[, k])
  h3b <- rasterFromXYZ(cbind(xx, z3b))
  writeRaster(h3b, filename = paste0(output_dir, "/pop_",
                                     rclass[k], "_upper", ".tif"), overwrite = TRUE,
              options = c(COMPRESS = "LZW"))
}

ras1<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output/Education/pop_CMR_population_v1_0_edu_None.tif")
ras2<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output/Education/pop_CMR_population_v1_0_edu_Primary.tif")
ras3<- rast("C:/Users/ccn1r22/OneDrive - University of Southampton/Documents/packages/main/jollofR_scripts/paper/output/Education/pop_CMR_population_v1_0_edu_Secondary.tif")

par(mfrow=c(3,3), mar=c(3,3,2,1))
plot(ras1, cex.lab=1.5,
     cex.axis=2) # visulize raster
plot(ras2)
plot(ras3)

#--------------------------------------------------------------------
##### Cross validation
#-------------------------------------------------------------------

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

  output <- list(#MAPE = MAPE,
    MAE  = MAE ,
    #BIAS = abs(BIAS),
    RMSE = RMSE,
    corr = corr)
  return(output)
}



seed = 1325
set.seed(seed)


dtt = admin_data
N <- nrow(dtt)

n.folds <- 10
#  k_fold


######

table(ind_train <- factor(sample(x = rep(1:n.folds,c(rep(35,9), 44)),  # Sample IDs for training data
                                 size = N))) # to allow 10 folds

table(as.numeric(ind_train))
dtt$k_fold <- as.numeric(ind_train)



k_uniq <-sort(unique(dtt$k_fold))

#---------------------------------------------------------------
#                   in-sample
#---------------------------------------------------------------

grid_list_in <- list()
for(u in 1:length(k_uniq))
{
  #u =1

  print(paste0("in-sample cross-validation using fold ", u, sep=""))
  test_ind <- which(dtt$k_fold==k_uniq[u])

  dim(test <- dtt[test_ind, ]) #---test set for fold i
  age_train <- age_dat

  for(i in 1:length(age_classes))
  {
    # i=2
    # Disaggregate Each age group count by sex
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.01))) # using PC prior
    print(paste(paste0("(",i,")"),paste0(age_classes[i], " model is running")))

    age_train[,colnames(age_train)[i]] <- round(age_train[,i])

    covars <- age_train[,c(cov_names,"set_typ", "admin_id")]; dim(covars) ##---Population density

    # stack for age
    stk <- inla.stack(data=list(y=age_train[,i], n=age_train$total), #the response

                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)

                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')


    print(form_age <- as.formula(paste0("y", " ~ ",
                                        paste(c("-1","Intercept", cov_names), collapse = " + "),
                                        " +   f(admin_id, model = 'iid', hyper = prior.prec)+
                                      f(set_typ, model = 'iid', hyper = prior.prec) +  f(s, model = spde)"))) # This is how you include  a typical random effect.
    mod_age <-inla(form_age, #the formula
                   data=inla.stack.data(stk,spde=spde),  #the data stack
                   family= 'binomial', Ntrials = n,  #which family the data comes from
                   control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                   control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                   verbose = FALSE)

    summary(mod_age)


    indx <-inla.stack.index(stk, "est")$data #---extract the data location indices
    prop_dt[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"mean"]),6)
    prop_dtL[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.025quant"]),6)
    prop_dtU[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.975quant"]),6)

    pred_dt[,i] = round(prop_dt[,i]*age_train$pop)
    pred_dtL[,i] = round(prop_dtL[,i]*age_train$pop)
    pred_dtU[,i] = round(prop_dtU[,i]*age_train$pop)

    f_dat[,colnames(f_dat)[i]] <- round(f_dat[,i]) # ensure counts are integers

    # sex disaggregation
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


    f.pred_dt[,i] = round(f.prop_dt[,i]*pred_dt[,i])# female counts
    f.pred_dtL[,i] = round(f.prop_dt[,i]*pred_dtL[,i])# female counts - lower
    f.pred_dtU[,i] = round(f.prop_dt[,i]*pred_dtU[,i])# female counts  upper

    m.pred_dt[,i] = pred_dt[,i] - f.pred_dt[,i] # male cunts
    m.pred_dtL[,i] = pred_dtL[,i] - f.pred_dtL[,i] # male cunts
    m.pred_dtU[,i] = pred_dtU[,i] - f.pred_dtU[,i] # male cunts

  }


  # disaggregate grid pop counts based on the grid props
  for (v in test_ind) {

    dim(grid_df <- grid_data[grid_data$admin_id == v, ])
    ids <- which(grid_data$admin_id == v)


    for (j in 1:length(age_classes)) {

      print(paste(paste0("age class ", j, " of admin ",
                         v, " is running")))

      # type 1 - sprinkle
      # admin disaggregated
      prop_grd2[ids, j] <- prop_dt[v, j]
      pop_grd2[ids, j] <- round(prop_grd2[ids, j] * grid_df$total,5)
      # Lower
      prop_grdL[ids, j] <- prop_dtL[v, j]
      pop_grdL[ids, j] <- round(prop_grdL[ids, j] * grid_df$total,2)

      #Upper
      prop_grdU[ids, j] <- prop_dtU[v, j]
      pop_grdU[ids, j] <- round(prop_grdU[ids, j] * grid_df$total,2)
    }
  }

  # Amin unit checks
  admin_cv <- admin_data
  admin_cv$tot2 <- apply(pred_dt, 1, sum, na.rm=T)
 
  # grid checks
  grid_cv <- grid_data
  grid_cv$tot2 <- apply(pop_grd2, 1, sum, na.rm=T)

  grid_list_in[[u]] <- unlist(model_metrics(grid_cv$total[grid_cv$admin_id%in%test_ind], 
                                            grid_cv$tot2[grid_cv$admin_id%in%test_ind]))
}
grid_list_in_dat <- do.call(rbind,grid_list_in)
metrics_in_grid <- apply(grid_list_in_dat, 2, mean,na.rm=T)



#  #---------------------------------------------------------------
#                   out-of-sample
#---------------------------------------------------------------

met_list_out <- list()
grid_list_out <- list()
for(u in 1:length(k_uniq))
{
  #u =1

  print(paste0("out-of-sample cross-validation using fold ", u, sep=""))
  test_ind <- which(dtt$k_fold==k_uniq[u])

  dim(test <- dtt[test_ind, ]) #---test set for fold i
  age_train <- age_dat

  age_train[test_ind,age_classes] =NA
  age_train$total <- round(age_dat$total)

  age_classes <- names(age_train %>% dplyr::select(starts_with("age_")))
  for(i in 1:length(age_classes))
  {
    # i=1
    # Disaggregate Each age group count by sex
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.01))) # using PC prior
    print(paste(paste0("(",i,")"),paste0(age_classes[i], " model is running")))

    age_train[,colnames(age_train)[i]] <- round(age_train[,i])

    covars <- age_train[,c(cov_names,"set_typ", "admin_id")]; dim(covars) ##---Population density

    # stack for age
    stk <- inla.stack(data=list(y=age_train[,i], n=age_train$total), #the response

                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)

                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')


    print(form_age <- as.formula(paste0("y", " ~ ",
                                        paste(c("-1","Intercept", cov_names), collapse = " + "),
                                        " +   f(admin_id, model = 'iid', hyper = prior.prec)+
                                      f(set_typ, model = 'iid', hyper = prior.prec) +  f(s, model = spde)"))) # This is how you include  a typical random effect.
    mod_age <-inla(form_age, #the formula
                   data=inla.stack.data(stk,spde=spde),  #the data stack
                   family= 'binomial', Ntrials = n,  #which family the data comes from
                   control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                   control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                   verbose = FALSE)

    summary(mod_age)


    indx <-inla.stack.index(stk, "est")$data #---extract the data location indices
    prop_dt[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"mean"]),6)
    prop_dtL[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.025quant"]),6)
    prop_dtU[,i] = round(plogis(mod_age$summary.linear.predictor[indx,"0.975quant"]),6)

    pred_dt[,i] = round(prop_dt[,i]*age_train$pop)
    pred_dtL[,i] = round(prop_dtL[,i]*age_train$pop)
    pred_dtU[,i] = round(prop_dtU[,i]*age_train$pop)

    f_dat[,colnames(f_dat)[i]] <- round(f_dat[,i]) # ensure counts are integers

    # sex disaggregation
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


    f.pred_dt[,i] = round(f.prop_dt[,i]*pred_dt[,i])# female counts
    f.pred_dtL[,i] = round(f.prop_dt[,i]*pred_dtL[,i])# female counts - lower
    f.pred_dtU[,i] = round(f.prop_dt[,i]*pred_dtU[,i])# female counts  upper

    m.pred_dt[,i] = pred_dt[,i] - f.pred_dt[,i] # male cunts
    m.pred_dtL[,i] = pred_dtL[,i] - f.pred_dtL[,i] # male cunts
    m.pred_dtU[,i] = pred_dtU[,i] - f.pred_dtU[,i] # male cunts

  }


  # disaggregate grid pop counts based on the grid props
  for (v in test_ind) {

    dim(grid_df <- grid_data[grid_data$admin_id == v, ])
    ids <- which(grid_data$admin_id == v)


    for (j in 1:length(age_classes)) {

      print(paste(paste0("age class ", j, " of admin ",
                         v, " is running")))

      # type 1 - sprinkle
      # admin disaggregated
      prop_grd2[ids, j] <- prop_dt[v, j]
      pop_grd2[ids, j] <- round(prop_grd2[ids, j] * grid_df$total,5)
      # Lower
      prop_grdL[ids, j] <- prop_dtL[v, j]
      pop_grdL[ids, j] <- round(prop_grdL[ids, j] * grid_df$total,2)

      #Upper
      prop_grdU[ids, j] <- prop_dtU[v, j]
      pop_grdU[ids, j] <- round(prop_grdU[ids, j] * grid_df$total,2)
    }
  }

  # Amin unit checks
  admin_cv <- admin_data
  admin_cv$tot2 <- apply(pred_dt, 1, sum, na.rm=T)


  # grid checks
  grid_cv <- grid_data
  grid_cv$tot2 <- apply(pop_grd2, 1, sum, na.rm=T)

  grid_list_out[[u]] <- unlist(model_metrics(grid_cv$total[grid_cv$admin_id%in%test_ind], 
                                             grid_cv$tot2[grid_cv$admin_id%in%test_ind]))
}
grid_list_out_dat <- do.call(rbind,grid_list_out)
metrics_out_grid <- apply(grid_list_out_dat, 2, mean,na.rm=T)

mets_cv <- rbind(metrics_in_grid, metrics_out_grid)

write.csv(mets_cv, paste0(output_path, "/cross_validated_metrics.csv"), row.names=F)


#---------------------------------------------------------------------
### Cross-validation for educational level disaggregation
#-------------------------------------------------------------------------
dtt = admin_data
N <- nrow(dtt)

n.folds <- 10
#  k_fold


######

table(ind_train <- factor(sample(x = rep(1:n.folds,c(rep(35,9), 44)),  # Sample IDs for training data
                                 size = N))) # to allow 10 folds

table(as.numeric(ind_train))
dtt$k_fold <- as.numeric(ind_train)



k_uniq <-sort(unique(dtt$k_fold))

#---------------------------------------------------------------
#                   in-sample
#---------------------------------------------------------------

grid_edu_list_in <- list()

for(u in 1:length(k_uniq))
{
  #u =1
  
  print(paste0("in-sample cross-validation using fold ", u, sep=""))
  test_ind <- which(dtt$k_fold==k_uniq[u])
  
  dim(test <- dtt[test_ind, ]) #---test set for fold i
  edu_train <- edu_dat
  
  for(i in 1:length(edu_classes))
  {
    # i=2
    # Disaggregate Each age group count by sex
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.01))) # using PC prior
    print(paste(paste0("(",i,")"),paste0(edu_classes[i], " model is running")))
    
    edu_train[,colnames(edu_train)[i]] <- round(edu_train[,i])
    
    covars <- edu_dat[,c(cov_names, "admin_id")]; dim(covars) ##---Population density
    
    # stack for age
    stk <- inla.stack(data=list(y=edu_train[,i], n=edu_train$total), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')
    
    
    
    print(form_edu <- as.formula(paste0("y", " ~ ",
                                        paste(c("-1","Intercept", cov_names), collapse = " + "),
                                        " +  f(admin_id, model = 'iid', hyper = prior.prec)+
                                f(s, model = spde)"))) # This is how you include  a typical random effect.
    
    mod_edu <-inla(form_edu, #the formula
                   data=inla.stack.data(stk,spde=spde),  #the data stack
                   family= 'binomial', Ntrials = n,  #which family the data comes from
                   control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                   control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                   verbose = FALSE)
    
    summary(mod_edu)
    
    
    prop_edu_dt[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"mean"]),6)
    prop_edu_dtL[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.025quant"]),6)
    prop_edu_dtU[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.975quant"]),6)
    
    pred_edu_dt[,i] = round(prop_edu_dt[,i]*edu_dat$pop)
    pred_edu_dtL[,i] = round(prop_edu_dtL[,i]*edu_dat$pop)
    pred_edu_dtU[,i] = round(prop_edu_dtU[,i]*edu_dat$pop)
    
    summary(mod_edu)
  }
  
  
  # disaggregate grid pop counts based on the grid props
  for (v in test_ind) {
    
    dim(grid_df <- grid_data[grid_data$admin_id == v, ])
    ids <- which(grid_data$admin_id == v)
    
    
    for (j in 1:length(edu_classes)) {
      
      print(paste(paste0("educational level class ", j, " of admin ",
                         v, " is running")))
      # admin disaggregated
      
      # Mean
      prop_edu_grd2[ids, j] <- prop_edu_dt[v, j]
      pop_edu_grd2[ids, j] <- round(prop_edu_grd2[ids, j] * grid_df$total,2)
      
      # Lower
      prop_edu_grdL[ids, j] <- prop_edu_dtL[v, j]
      pop_edu_grdL[ids, j] <- round(prop_edu_grdL[ids, j] * grid_df$total,2)
      
      #Upper
      prop_edu_grdU[ids, j] <- prop_edu_dtU[v, j]
      pop_edu_grdU[ids, j] <- round(prop_edu_grdU[ids, j] * grid_df$total,2)
    }
  }
  
  
  
  # grid checks
  grid_edu_cv <- grid_data
  grid_edu_cv$tot2 <- apply(pop_edu_grd2, 1, sum, na.rm=T)
  
  grid_edu_list_in[[u]] <- unlist(model_metrics(grid_edu_cv$total[grid_edu_cv$admin_id%in%test_ind], 
                                                grid_edu_cv$tot2[grid_edu_cv$admin_id%in%test_ind]))
}
grid_edu_list_in_dat <- do.call(rbind,grid_edu_list_in)
metrics_edu_in_grid <- apply(grid_edu_list_in_dat, 2, mean,na.rm=T)



#  #---------------------------------------------------------------
#                   out-of-sample
#---------------------------------------------------------------

grid_edu_list_out <- list()
for(u in 1:length(k_uniq))
{
  #u =1
  
  print(paste0("out-of-sample cross-validation using fold ", u, sep=""))
  test_ind <- which(dtt$k_fold==k_uniq[u])
  
  dim(test <- dtt[test_ind, ]) #---test set for fold i
  edu_train <- edu_dat
  
  edu_train[test_ind,edu_classes] =NA
  edu_train$total <- round(edu_dat$total)
  
  edu_classes <- names(edu_train %>% dplyr::select(starts_with("edu_")))
  for(i in 1:length(edu_classes))
  {
    # i=2
    # Disaggregate Each age group count by sex
    prior.prec <- list(prec = list(prior = "pc.prec",
                                   param = c(1, 0.01))) # using PC prior
    print(paste(paste0("(",i,")"),paste0(edu_classes[i], " model is running")))
    
    edu_train[,colnames(edu_train)[i]] <- round(edu_train[,i])
    
    covars <- edu_dat[,c(cov_names, "admin_id")]; dim(covars) ##---Population density
    
    # stack for age
    stk <- inla.stack(data=list(y=edu_train[,i], n=edu_train$total), #the response
                      
                      A=list(A,1),  #the A matrix; the 1 is included to make the list(covariates)
                      
                      effects=list(c(list(Intercept=1), #the Intercept
                                     iset),  #the spatial index
                                   #the covariates
                                   list(covars)
                      ),
                      #this is a quick name so you can call upon easily
                      tag='est')
    
    
    
    print(form_edu <- as.formula(paste0("y", " ~ ",
                                        paste(c("-1","Intercept", cov_names), collapse = " + "),
                                        " +  f(admin_id, model = 'iid', hyper = prior.prec)+
                                f(s, model = spde)"))) # This is how you include  a typical random effect.
    
    mod_edu <-inla(form_edu, #the formula
                   data=inla.stack.data(stk,spde=spde),  #the data stack
                   family= 'binomial', Ntrials = n,  #which family the data comes from
                   control.predictor=list(A=inla.stack.A(stk),compute=TRUE),  #compute gives you the marginals of the linear predictor
                   control.compute = list(dic = TRUE, waic = TRUE, cpo=TRUE,config = TRUE), #model diagnostics and config = TRUE gives you the GMRF
                   verbose = FALSE)
    
    summary(mod_edu)
    
    
    prop_edu_dt[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"mean"]),6)
    prop_edu_dtL[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.025quant"]),6)
    prop_edu_dtU[,i] = round(plogis(mod_edu$summary.linear.predictor[indx,"0.975quant"]),6)
    
    pred_edu_dt[,i] = round(prop_edu_dt[,i]*edu_dat$pop)
    pred_edu_dtL[,i] = round(prop_edu_dtL[,i]*edu_dat$pop)
    pred_edu_dtU[,i] = round(prop_edu_dtU[,i]*edu_dat$pop)
    
    summary(mod_edu)
  }
  
  
  # disaggregate grid pop counts based on the grid props
  for (v in test_ind) {
    
    dim(grid_df <- grid_data[grid_data$admin_id == v, ])
    ids <- which(grid_data$admin_id == v)
    
    
    for (j in 1:length(edu_classes)) {
      
      print(paste(paste0("educational level class ", j, " of admin ",
                         v, " is running")))
      # admin disaggregated
      
      # Mean
      prop_edu_grd2[ids, j] <- prop_edu_dt[v, j]
      pop_edu_grd2[ids, j] <- round(prop_edu_grd2[ids, j] * grid_df$total,2)
      
      # Lower
      prop_edu_grdL[ids, j] <- prop_edu_dtL[v, j]
      pop_edu_grdL[ids, j] <- round(prop_edu_grdL[ids, j] * grid_df$total,2)
      
      #Upper
      prop_edu_grdU[ids, j] <- prop_edu_dtU[v, j]
      pop_edu_grdU[ids, j] <- round(prop_edu_grdU[ids, j] * grid_df$total,2)
    }
  }
  
  
  # grid checks
  grid_edu_cv <- grid_data
  grid_edu_cv$tot2 <- apply(pop_grd2, 1, sum, na.rm=T)
  
  grid_edu_list_out[[u]] <- unlist(model_metrics(grid_edu_cv$total[grid_edu_cv$admin_id%in%test_ind], 
                                                 grid_edu_cv$tot2[grid_edu_cv$admin_id%in%test_ind]))
}
grid_edu_list_out_dat <- do.call(rbind,grid_edu_list_out)
metrics_edu_out_grid <- apply(grid_edu_list_out_dat, 2, mean,na.rm=T)

mets_edu_cv <- rbind(metrics_edu_in_grid, metrics_edu_out_grid)
#                       MAE      RMSE      corr
#metrics_edu_in_grid  0.1428564 0.5144165 0.9998678
#metrics_edu_out_grid 0.1269514 0.4381988 0.9999062
save.image(paste0(output_path, "CMR_application_output.RData"))
