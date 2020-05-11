##Statistical analysis to accompany KPSC DME Epidemiology Paper, 2008-2018 EHR data
##Code and analyses by Joan Casey
##Updated May 11, 2020

#load libraries
library(tidyverse)
library(sp)
library(data.table)
library(viridis)
library(ggcorrplot)
library(sf)
library(gee)
library(broom)
library(RColorBrewer)
library(spdep)
library(hglm)
library(scales)
library(ggspatial)
library(ggmap)
library(ggsn)
library(egg)
library(table1)

#Scale by SD
scale_this <- function(x){
  (x) / sd(x, na.rm=TRUE)
}

#Read in ehr data
dat <- read_csv("uni_mrn_dme_demo_ses3_medical_bp.csv")
head(dat)
colnames(dat)
table(dat$dme_yes)
n_distinct(dat$ZIP) #1260 unique ZIP codes
n_distinct(dat$studyID) #249297

#Link in DME identifiers
type_dme <- read_csv("dme_categories_year_id.csv")
dat <- left_join(dat, type_dme, by=c("studyID"="studyid", "year"="year"))

#Dropping if age is missing because it means person wasn't yet born in that year
table(dat$agecat)
dat <- dat %>% filter(agecat!="Missing")
n_distinct(dat$studyID) #249297
n_distinct(dat$ZIP) #1260
start <- 249297

#Selection criteria:  1. Has a ZIP code; 2. in S. Cali ZIP; 3. Total population >30; 4. >30 KPSC patients in ZIP
#Read in denominators by ZIP code

#1 Dropping data without a ZIP code
summary(dat$ZIP)
dat <- dat %>% drop_na("ZIP")
n_distinct(dat$ZIP) #1259 unique ZIP codes
n_distinct(dat$studyID) #244563

#2 located in S. Cali
zip_kpsc <- read_csv("zips_kpsc.csv")
dat <- right_join(dat, zip_kpsc, by = c("ZIP"="zip_code"))

#3 remove men with breast pumps (n = 1) and women outside reproductive age (15-49, n = 0)
dat <- dat %>% mutate(drop_male_bp = ifelse(bp_overall==1 & GENDER=="M",1,0))
table(dat$drop_male_bp)
dat <- dplyr::filter(dat, drop_male_bp!=1)

#4 women not repro age getting breast pumps
dat <- dat %>% mutate(drop_f_age = ifelse(bp_overall==1 & (age<15 | age>49),1,0))
dat <- dplyr::filter(dat, drop_f_age!=1)

#5 dropping ZIPs with <30 KPSC members
denom <- read_csv("dme_denom_zip_kpsc_all_final.csv")
dat <- left_join(dat, denom, by = c("ZIP" = "zip", "year"="year"))

#6 dropping when census variables are missing
dat <- dat %>% drop_na(age_gt65_p  ,   asian_p    ,    black_p     ,  med_inc,
                       edu_lt_hs_p,hispan_p,       ling_iso_p   ,    pov_p    ,     
                       rent_p  ,      rural_p  ,      snap_p     ,   total_pop  ,    unemploy_p   ,  white_p)

#Number of people using breast pumps vs. other DME
dat %>% dplyr::filter(bp_overall==1) %>% summarise(count = n_distinct(studyID)) #96534
dat %>% dplyr::filter(bp_overall==1) %>% summarise(count = n_distinct(ZIP)) #648 ZIP
dat %>% dplyr::filter(bp_overall==1) %>% summarise(count = n_distinct(studyID))/start7*100
dat %>% dplyr::filter(notbp_overall==1) %>% summarise(count = n_distinct(studyID)) #147295
dat %>% dplyr::filter(notbp_overall==1) %>% summarise(count = n_distinct(ZIP)) #673 ZIP
dat %>% dplyr::filter(notbp_overall==1) %>% summarise(count = n_distinct(studyID)) /start7*100

#Link distance to nearest healthcare facility#
health_dist <- read_csv("facility_zip_dist_m.csv")
dat <- left_join(dat, health_dist, by = c("ZIP" = "ZIP"))

###Add type of DME###
type_dme <- read_csv("type_of_dme_kpsc.csv")
dat <- left_join(dat, type_dme, by=c("studyID"="studyid", "year"="year"))
dat <- dat %>% mutate(age_cat3 = case_when(age < 18 ~ 1, age >=18 & age <65 ~2, age >=65 ~3)) #create three age bins
dat <- dat %>% mutate(age_cat6 = case_when(age < 18 & GENDER=="F" ~ 1, 
                                           age < 18 & GENDER=="M" ~ 2,
                                           (age >=18 & age <65) & GENDER=="F" ~3,
                                           (age >=18 & age <65) & GENDER=="M" ~4,
                                           age >=65 & GENDER=="F" ~5,
                                           age >=65 & GENDER=="M" ~6)) #create 6 sex / age bins

#Create racial/ethnic categories
dat <- dat %>% mutate(black = ifelse(race_eth=="Black",1,0))
dat <- dat %>% mutate(asian = ifelse(race_eth=="Asian",1,0))
dat <- dat %>% mutate(hispan = ifelse(race_eth=="Hispanic",1,0))
dat <- dat %>% mutate(white = ifelse(race_eth=="White",1,0))
dat <- dat %>% mutate(other = ifelse(race_eth=="Other",1,0))

#Variable indicating use of 2+ categories of DME
dat <- dat %>% mutate(dme_2plus = ifelse(dme_tot_types>1, 1, 0))
table(dat$dme_2plus)

########################################
#Check residual spatial autocorrelation#
########################################
##Breast pumps##
dat_bp <- filter(dat, bp_overall==1)
n_distinct(dat_bp$studyID) # 96798

#Create demographics dat_bpaframe
census_dem <- dplyr::select(dat_bp, "ZIP"  ,   "asian_p"    ,    "black_p"     ,  "med_inc",
                            "edu_lt_hs_p" ,   "gisjoin"   ,"hispan_p",       "ling_iso_p" ,    "med_inc"    ,    "pov_p"    ,     
                            "rent_p"   ,      "rural_p"  ,      "snap_p"      ,   "total_pop"  ,    "unemploy_p"   ,  "white_p")
census_dem <- distinct(census_dem)

#Summarize dme to year/zip
dat_year_zip <- dat %>% group_by(ZIP, year) %>%
  dplyr::summarize(tot_dme_bp = sum(bp_yes))

#add denominator
dat_year_zip <- left_join(dat_year_zip, denom, by = c("ZIP" = "zip", "year"="year"))

#Scale so they're per 10% change in Medicaid prevalence
dat_year_zip$kpsc_f_p <- dat_year_zip$kpsc_f/dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_medicaid_p <- dat_year_zip$kpsc_medicaid / dat_year_zip$kpsc_n *10
dat_year_zip$kpsc_f_p <- dat_year_zip$kpsc_f /dat_year_zip$kpsc_f * 10
dat_year_zip$kpsc_black_p <- dat_year_zip$kpsc_black / dat_year_zip$kpsc_n * 10
dat_year_zip$kpsc_white_p <- dat_year_zip$kpsc_white /dat_year_zip$kpsc_n * 10
dat_year_zip$kpsc_asian_p <- dat_year_zip$kpsc_asian /dat_year_zip$kpsc_n * 10
dat_year_zip$kpsc_hispan_p <- dat_year_zip$kpsc_hispan/dat_year_zip$kpsc_n * 10
dat_year_zip$kpsc_other_p <- dat_year_zip$kpsc_other / dat_year_zip$kpsc_n *10

#Add demographics 
dat_year_zip <- left_join(dat_year_zip, census_dem, by = c("ZIP" = "ZIP"))
#Scale all demographics, so we interpret coefficient as "for a 10% increase in X"
divide10 <- function(census_var) {
  census_var/10
}

dat_year_zip$med_inc <-scale_this(dat_year_zip$med_inc ) 
dat_year_zip$edu_lt_hs_p <-scale_this(dat_year_zip$edu_lt_hs_p ) 
dat_year_zip$ling_iso_p <-scale_this(dat_year_zip$ling_iso_p ) 
dat_year_zip$pov_p <-scale_this(dat_year_zip$pov_p ) 
dat_year_zip$rent_p <-scale_this(dat_year_zip$rent_p ) 
dat_year_zip$rural_p <-scale_this(dat_year_zip$rural_p ) 
dat_year_zip$snap_p <-scale_this(dat_year_zip$snap_p ) 
dat_year_zip$total_pop <-scale_this(dat_year_zip$total_pop ) 
dat_year_zip$rent_p <- scale_this(dat_year_zip$rent_p)
dat_year_zip$unemploy_p <-scale_this(dat_year_zip$unemploy_p )

#Scale for 10% increase 
dat_year_zip$black_p <-dat_year_zip$black_p  / 10
dat_year_zip$hispan_p <-dat_year_zip$hispan_p /10
dat_year_zip$white_p <-dat_year_zip$white_p  /10
dat_year_zip$asian_p <-dat_year_zip$asian_p /10
summary(dat_year_zip$hispan_p)

ZIP <- st_read("usps_zips_kpsc_v2.shp")
head(ZIP)
dat_year_zip$ZIP <- factor(dat_year_zip$ZIP)
ZIP_2 <- left_join(ZIP, dat_year_zip, by=c("ZIP_CODE"="ZIP"))
table(ZIP_2$year)
ZIP_2 <- dplyr::select(ZIP_2, year, kpsc_mean_age , kpsc_black_p  ,kpsc_hispan_p ,
                       kpsc_medicaid_p   , med_inc   ,rural_p , 
                       ling_iso_p , pov_p  , edu_lt_hs_p ,rural_p ,
                       snap_p , unemploy_p, kpsc_f1549, tot_dme_bp, ZIP_CODE)
ZIP_2 <- ZIP_2 %>% drop_na()

#Show final model for breast pumps
model0 <- hglm(fixed = tot_dme_bp ~  as.factor(year) + kpsc_mean_age + kpsc_black_p  + kpsc_hispan_p +
                 kpsc_medicaid_p  + pov_p  + rent_p +
                 snap_p, random = ~1|ZIP_CODE, offset = log(kpsc_f1549),
               data = ZIP_2,
               family = poisson(link=log)) 
summary(model0)

#Spatial pattern of residuals
ZIP_2$residuals <- model0$res

ZIP_2 %>% 
  ggplot() +
  geom_sf(aes(fill = residuals),  lwd = 0) +
  scale_fill_viridis("Residuals", na.value = "grey")+
  theme_grey(base_size = 14) +theme(legend.position="bottom", legend.box = "horizontal")
#Random!

#Moran's I for residuals
#Create neighbor matrix
nb <- poly2nb(ZIP_2, queen=FALSE)

#Assigning weights
nb_lw <- nb2listw(nb, style="B")

moran.test(ZIP_2$residuals, nb_lw) #p-value = 1 indicating no support non-randomly distributed residuals

##############################
##Non-breast pump area-level##
##############################
dat_bp <- filter(dat, notbp_overall==1)

#Create demographics dat_bpaframe
census_dem <- dplyr::select(dat_bp, "ZIP"  ,   "asian_p"    ,    "black_p"     ,  "med_inc",
                            "edu_lt_hs_p" ,   "gisjoin"   ,"hispan_p",       "ling_iso_p" ,    "med_inc"    ,    "pov_p"    ,     
                            "rent_p"   ,      "rural_p"  ,      "snap_p"      ,   "total_pop"  ,    "unemploy_p"   ,  "white_p", "dist_m")
census_dem <- distinct(census_dem)

#Summarize dme to year/zip
dat_year_zip <- dat %>% group_by(ZIP, year) %>%
  dplyr::summarize(tot_dme_bp = sum(notbp_yes))

#add denominator
dat_year_zip <- left_join(dat_year_zip, denom, by = c("ZIP" = "zip", "year"="year"))

#Scale so they're per 10% change in Medicaid prevalence
dat_year_zip$kpsc_f_p <- dat_year_zip$kpsc_f/dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_medicaid_p <- dat_year_zip$kpsc_medicaid / dat_year_zip$kpsc_n *100
dat_year_zip$kpsc_f_p <- dat_year_zip$kpsc_f /dat_year_zip$kpsc_f * 100
dat_year_zip$kpsc_black_p <- dat_year_zip$kpsc_black / dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_white_p <- dat_year_zip$kpsc_white /dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_asian_p <- dat_year_zip$kpsc_asian /dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_hispan_p <- dat_year_zip$kpsc_hispan/dat_year_zip$kpsc_n * 100
dat_year_zip$kpsc_other_p <- dat_year_zip$kpsc_other / dat_year_zip$kpsc_n *100
dat_year_zip$dist_km <-dat_year_zip$dist_m/1000

#Add demographics 
dat_year_zip <- left_join(dat_year_zip, census_dem, by = c("ZIP" = "ZIP"))
#Scale all demographics, so we interpret coefficient as "for a 10% increase in X"
divide10 <- function(census_var) {
  census_var/10
}

dat_year_zip$med_inc <-scale_this(dat_year_zip$med_inc ) 
dat_year_zip$edu_lt_hs_p <-scale_this(dat_year_zip$edu_lt_hs_p ) 
dat_year_zip$ling_iso_p <-scale_this(dat_year_zip$ling_iso_p ) 
dat_year_zip$pov_p <-scale_this(dat_year_zip$pov_p ) 
dat_year_zip$rent_p <-scale_this(dat_year_zip$rent_p ) 
dat_year_zip$rural_p <-scale_this(dat_year_zip$rural_p ) 
dat_year_zip$snap_p <-scale_this(dat_year_zip$snap_p ) 
dat_year_zip$total_pop <-scale_this(dat_year_zip$total_pop ) 
dat_year_zip$rent_p <- scale_this(dat_year_zip$rent_p)
dat_year_zip$unemploy_p <-scale_this(dat_year_zip$unemploy_p )

#Scale for 10% increase 
dat_year_zip$black_p <-dat_year_zip$black_p  / 10
dat_year_zip$hispan_p <-dat_year_zip$hispan_p /10
dat_year_zip$white_p <-dat_year_zip$white_p  /10
dat_year_zip$asian_p <-dat_year_zip$asian_p /10
summary(dat_year_zip$hispan_p)

ZIP <- st_read("usps_zips_kpsc_v2.shp")
head(ZIP)
dat_year_zip$ZIP <- factor(dat_year_zip$ZIP)
ZIP_2 <- left_join(ZIP, dat_year_zip, by=c("ZIP_CODE"="ZIP"))
table(ZIP_2$year)
ZIP_2 <- dplyr::select(ZIP_2, year, kpsc_mean_age , kpsc_black_p  ,kpsc_hispan_p ,
                       kpsc_medicaid_p   , med_inc   , rent_p , 
                       ling_iso_p , pov_p  , edu_lt_hs_p ,rural_p ,
                       snap_p , unemploy_p, kpsc_n, tot_dme_bp, ZIP_CODE)
ZIP_2 <- ZIP_2 %>% drop_na()

#final model for non-breast pumps
model0 <- hglm(fixed = tot_dme_bp ~  as.factor(year) + kpsc_mean_age + kpsc_black_p  + kpsc_hispan_p +
                 kpsc_medicaid_p  + ling_iso_p + unemploy_p + edu_lt_hs_p + rent_p,
                 random = ~1|ZIP_CODE, offset = log(kpsc_n),
               data = ZIP_2,
               family = poisson(link=log)) 
summary(model0)

#Spatial pattern of residuals
ZIP_2$residuals <- model0$res

ZIP_2 %>% 
  ggplot() +
  geom_sf(aes(fill = residuals),  lwd = 0) +
  scale_fill_viridis("Residuals", na.value = "grey", breaks = c(-150, -50, 0, 50, 100))+
  theme_grey(base_size = 14) +theme(legend.position="bottom", legend.box = "horizontal")
#Random!

#Moran's I for residuals
#Create neighbor matrix
nb <- poly2nb(ZIP_2, queen=FALSE)

#Assigning weights
nb_lw <- nb2listw(nb, style="B")

moran.test(ZIP_2$residuals, nb_lw) #p-value = 1 indicating no support non-randomly distributed residuals

###################
#TABLE 2 MODELS   #
###################
####BREAST PUMP####
#Year, medicaid, age
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + 
                 kpsc_medicaid_p + med_inc + kpsc_black_p + kpsc_hispan_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#Black race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  +
                 kpsc_medicaid_p + med_inc + kpsc_black_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#Asian race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  +
                 kpsc_medicaid_p + med_inc + kpsc_asian_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#Hispanic race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  +
                 kpsc_medicaid_p + med_inc + kpsc_hispan_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#White race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  +
                 kpsc_medicaid_p + med_inc + kpsc_white_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#Other race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  +
                 kpsc_medicaid_p + med_inc + kpsc_other_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

######################
#AREA-LEVEL VARIABLES#
######################
#Poisson regression for each of the 8 SES factor individually 
####BREAST PUMP####
#1 Median household income
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#2 Poverty
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + pov_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

#3 SNAP
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + snap_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
#4 ling iso
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + ling_iso_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
#5 renters
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + rent_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
#6 education
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + edu_lt_hs_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
#7 rural
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + rural_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
#8 unemployment
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + unemploy_p, random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))
##############################################################
#Fully adjusted model picking important predictors with LASSO
# Split the data into training and test set BREAST PUMPS
##############################################################
dat_year_zip_lass <- dat_year_zip %>% dplyr::select(year, kpsc_mean_age, kpsc_mean_age ,kpsc_black_p , kpsc_hispan_p,
                                                    kpsc_medicaid_p,  med_inc,
                                                    edu_lt_hs_p, rural_p, ling_iso_p, pov_p, rent_p, rural_p,
                                                    snap_p,unemploy_p,  tot_dme_bp, kpsc_f1549)

dat_year_zip_lass$year <- factor(dat_year_zip_lass$year, levels=c(2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018))
dat_year_zip_lass <- dat_year_zip_lass[,2:16]

dat_year_zip_lass <- na.omit(dat_year_zip_lass)

set.seed(317)
library(glmnet)
library(caret)

# Offset 
offset <- log(dat_year_zip_lass$kpsc_f1549)
dat_year_zip_lass <- dat_year_zip_lass[,-15]

# Predictor variables
x <- model.matrix(tot_dme_bp ~., dat_year_zip_lass)

# Outcome variable
y <- dat_year_zip_lass$tot_dme_bp

train <- sample(1:nrow(x), nrow(x)*.8)
test <- (-train)
y.test <- y[test]
dim (x[train,-1])
length(y[train])

#Identifying best lamba
cv <- cv.glmnet(x[train,-1], y[train],alpha=1, family = "poisson", offset=offset[train])
plot(cv)
opt.lam = c(cv$lambda.min, cv$lambda.1se)
coef(cv, s = opt.lam)

# Display the best lambda value
cv$lambda.1se

#LASSO model
# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, lambda = cv$lambda.1se, family = "poisson", offset=offset)
# Display regression coefficients
coef(model) # omitted med_inc, education, rural, living iso, unemployment

model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p +  pov_p + snap_p +
                 rent_p , 
               random = ~1|ZIP, offset = log(kpsc_f1549), 
               dat=dat_year_zip, family = poisson(link = log))

###################
#SOCIODEMOGRAPHICS#
###################
####NONBREAST PUMP#
#Year, medicaid, age
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p + 
                 kpsc_medicaid_p + med_inc + kpsc_black_p + kpsc_hispan_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))

#Black race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p +
                 kpsc_medicaid_p + med_inc + kpsc_black_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#Asian race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p +
                 kpsc_medicaid_p + med_inc + kpsc_asian_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#Hispanic race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p +
                 kpsc_medicaid_p + med_inc + kpsc_hispan_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#White race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p +
                 kpsc_medicaid_p + med_inc + kpsc_white_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#Other race KPSC
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)  + kpsc_mean_age  + kpsc_f_p +
                 kpsc_medicaid_p + med_inc + kpsc_other_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))

######################
#AREA-LEVEL VARIABLES#
######################
####NONBREAST PUMP#
#Poisson regression for each of the 8 SES factor individually 
#1 Median household income
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#2 Poverty
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + pov_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#3 SNAP
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + snap_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#4 ling iso
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + ling_iso_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#5 renters
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year)+ kpsc_f_p  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + rent_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#6 education
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p  + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + edu_lt_hs_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#7 rural
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + rural_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#8 unemployment
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + unemploy_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#9 >=65
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + age_gt65_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))

#################
#AREA-LEVEL RACE#
#################
####NONBREAST PUMP#
#1 Black race/ethnicity
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc + black_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#2 Asian race/ethnicity
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc + asian_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#3 Hispanic race/ethnicity
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc + hispan_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))
#4 Hispanic race/ethnicity
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + kpsc_hispan_p +
                 kpsc_medicaid_p + med_inc + white_p, random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))

##############################################################
#Fully adjusted model picking important predictors with LASSO
# Split the data into training and test set NON BREAST PUMPS
##############################################################
dat_year_zip_lass <- dat_year_zip %>% dplyr::select(year, kpsc_mean_age, kpsc_mean_age ,kpsc_black_p , kpsc_hispan_p,
                                                    kpsc_medicaid_p,  med_inc, kpsc_f_p, dist_km, age_gt65_p,
                                                    edu_lt_hs_p, rural_p, ling_iso_p, pov_p, rent_p, rural_p,
                                                    snap_p,unemploy_p,  tot_dme_bp, kpsc_n)

dat_year_zip_lass$year <- factor(dat_year_zip_lass$year, levels=c(2008,2009,2010,2011,2012,2013,2014,2015,2016,2017,2018))
dat_year_zip_lass <- dat_year_zip_lass[,2:19]

dat_year_zip_lass <- na.omit(dat_year_zip_lass)

set.seed(317)
library(glmnet)
library(caret)

# Offset 
offset <- log(dat_year_zip_lass$kpsc_n)
dat_year_zip_lass <- dat_year_zip_lass[,-18]

# Predictor variables
x <- model.matrix(tot_dme_bp ~., dat_year_zip_lass)

# Outcome variable
y <- dat_year_zip_lass$tot_dme_bp

train <- sample(1:nrow(x), nrow(x)*.8)
test <- (-train)
y.test <- y[test]
dim (x[train,-1])
length(y[train])

#Identifying best lamba
cv <- cv.glmnet(x[train,-1], y[train],alpha=1, family = "poisson", offset=offset[train])
plot(cv)
opt.lam = c(cv$lambda.min, cv$lambda.1se)
coef(cv, s = opt.lam)

#LASSO model
# Fit the final model on the training data
model <- glmnet(x, y, alpha = 1, lambda = cv$lambda.1se, family = "poisson", offset=offset)
# Display regression coefficients
coef(model) ##Keeping ling_iso, renters, unemployment, < high school education, age >=65

#1-se lamba
model0 <- hglm(fixed = tot_dme_bp ~ as.factor(year) + kpsc_f_p + kpsc_mean_age + kpsc_black_p + 
                 kpsc_hispan_p + kpsc_medicaid_p + edu_lt_hs_p + 
                 ling_iso_p + rent_p + unemploy_p + age_gt65_p  , 
               random = ~1|ZIP, offset = log(kpsc_n), 
               dat=dat_year_zip, family = poisson(link = log))

###################
#TABLE 3 MODELS   #
###################
####BREAST PUMP####
####################################################################################
#####Days overall during the study period re: individual and area-level factors#####
####################################################################################
model0 <- hglm(fixed = days_overall ~ age + medicaid + race_f + med_inc, random = ~1|ZIP, 
               dat=dat_indiv_bp, family = gaussian(link = identity))

#####################
# NON BP ANALYSES ###
#####################
model0 <- hglm(fixed = days_overall ~ age + medicaid + sex + race_f + med_inc, random = ~1|ZIP, 
               dat=dat_indiv_2, family = gaussian(link = identity))


#########Logisitic regression for using 2+ DME categories in 1 year##########
model0 <- hglm(fixed = dme_2plus ~ age + medicaid + sex + race_f + med_inc, random = ~1|ZIP, 
               dat=dat_indiv_2, family = binomial(link = "logit"))

#################################
#INDIVIDUAL-LEVEL LASSO##########
#################################
####BREAST PUMP####
dat_year_zip_lass <- dat_indiv_bp %>% drop_na(days_overall) %>% 
  dplyr::select(days_overall, studyID, medicaid, age, race_f, med_inc,
                edu_lt_hs_p, rural_p, ling_iso_p, pov_p, rent_p, rural_p, dist_km,
                snap_p,unemploy_p)
dat_year_zip_lass <- dat_year_zip_lass[,-2]
set.seed(317)
library(glmnet)
library(caret)

# Predictor variables
x <- model.matrix(days_overall ~., dat_year_zip_lass)

# Outcome variable
y <- dat_year_zip_lass$days_overall

train <- sample(1:nrow(x), nrow(x)*.8)
test <- (-train)
y.test <- y[test]
dim (x[train,-1])
length(y[train])

#Identifying best lamba
cv <- cv.glmnet(x[train,-1], y[train],alpha=1, family = "gaussian")
plot(cv)
opt.lam = c(cv$lambda.min, cv$lambda.1se)
coef(cv, s = exp(0.5))

# Display the best lambda value
cv$lambda.1se #NO VARIABLES RETAINED!

####FOR NON-BREAST PUMP####
dat_year_zip_lass <- dat_indiv_2 %>% drop_na(days_overall) %>% 
  dplyr::select(days_overall, studyID, medicaid, age, race_f, med_inc, sex,
                edu_lt_hs_p, rural_p, ling_iso_p, pov_p, rent_p, rural_p, dist_km,
                snap_p,unemploy_p, age_gt65_p)
dat_year_zip_lass <- dat_year_zip_lass[,-2]
set.seed(317)
library(glmnet)
library(caret)

# Predictor variables
x <- model.matrix(days_overall ~., dat_year_zip_lass)

# Outcome variable
y <- dat_year_zip_lass$days_overall

train <- sample(1:nrow(x), nrow(x)*.8)
test <- (-train)
y.test <- y[test]
dim (x[train,-1])
length(y[train])

#Identifying best lamba
cv <- cv.glmnet(x[train,-1], y[train],alpha=1, family = "gaussian")
plot(cv)
opt.lam = c(cv$lambda.min, cv$lambda.1se)
coef(cv, s = opt.lam) # retains only unemployment

model0 <- hglm(fixed = days_overall ~ age + medicaid + race_f +
                 unemploy_p, 
               random = ~1|ZIP,  
               dat=dat_indiv_2, family = gaussian(link = "identity"))

