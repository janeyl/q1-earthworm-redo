#_______________________________________________________
#Project Info----
#_______________________________________________________
# Started Feb 24, 2024
# Janey R. Lienau
# Qualifying Exams: Question 1
# Analysis of Phillips et al. 2019 "Global distribution of earthworm diversity" DOI: 10.1126/science.aax4851
#_______________________________________________________
# 1) Contents----
#_______________________________________________________
# 1) Contents
# 2) Load libraries
# 3) Import Data
# 4) Cleaning Data
# 5) Richness Specification Curve Analysis
# 6) Abundance Specification Curve Analysis
# 7) Biomass Specification Curve Analysis
# 8) Regional Map
# 9) Random Forest Model
# 10) Applying predictors from real data to Soilsgrid data for fun
#_______________________________________________________
# 2) Load libraries----
#_______________________________________________________
#figure out how to make library(groundhog) work...
library(pacman)
p_load(RColorBrewer, 
       tidyverse, 
       cowplot,
       corrplot,
       car,
       specr,
       lme4,
       ggplot2,
       sf,
       rnaturalearth,
       rnaturalearthdata,
       randomForest)
#_______________________________________________________
# 3) Load data----
#_______________________________________________________
# this data is downloaded from public repository from this 
# Phillips, H. R. P. et al.(2019). Global distribution of earthworm diversity (Version 1.1) . 
# Dataset available at iDiv Data Repository. http://idata.idiv.de/ddm/Data/ShowData/1818
# the following csv files were downloaded on Monday Feb 24, 2025 at 12:45 pm and included in this repository
# 1) 1818_2_sWormModelData_meta-data_correction_2020April_complete.csv
# 2) 1818_2_sWormModelData_correction_2020April.csv
# 3) SiteData_sWorm_2021-02-18.csv was added from the dataset https://idata.idiv.de/ddm/Data/ShowData/1880?version=0 DOI: doi.org/10.1038/s41597-021-00912-z

#set the working directory to this repository to this repository
setwd("/Users/janeylienau/Desktop/GitHubRepository/q1-earthworm-redo")
#load in data
metadata.df <- read.csv("raw-data/1818_15_1818_2_sWormModelData_meta-data_correction_2020April_complete.csv") 
modeldata.df <- read.csv("raw-data/1818_15_1818_2_sWormModelData_correction_2020April.csv")
coords <- read.csv("raw-data/SiteData_sWorm_2021-02-18.csv")
#_______________________________________________________
#4) Clean Data----
#_______________________________________________________
#CHELSA Climatologies data based on coordinates notes for what these mean
#bio10_1 (annual mean temperature)
#bio10_4 (temperature seasonality)
#bio10_7 (temperature annual range) - Phillips et al used for richness
#bio10_12 (annual precipitation) 
#bio10_15 (precipitation seasonality)

#look at data
names(modeldata.df)
str(modeldata.df)
head(modeldata.df)

origData <- modeldata.df

#the first step of our analysis is to remove all the SoilGrid data and ESA because it's an intermediate
modeldata.df <- modeldata.df%>%
  select(-ESA,-PHIHOX,-CLYPPT,-SLTPPT,-SNDPPT,-CECSOL,-ORCDRC,-phFinal,-ClayFinal,-SandFinal,-SiltFinal,-OCFinal)
#remove useless columns
modeldata.df <- modeldata.df%>%
  select(-rowID,-file,-PH,-PH_Collection_Method, -PH_mean,-OC_mean,-CN_mean,-sand_silt_clay_mean,-SpeciesRichnessUnit,-Site_WetBiomass,-Site_WetBiomassUnits,-Site_AbundanceUnits, -Site_Abundance)
#clean up some weird things with cleaning function
clean_string <- function(x) {
  x <- iconv(x, to = "ASCII", sub = "")  # Convert to ASCII
  x <- gsub("[^\x01-\x7F]", "", x)   # Remove any remaining non-ASCII characters
  return(x)
}
modeldata.df$Study_Name <- sapply(modeldata.df$Study_Name, clean_string)
modeldata.df$Site_Name <- sapply(modeldata.df$Site_Name, clean_string)
#_______________________________________________________#_______________________________________________________
# 5) RICHNESS----
#_______________________________________________________#_______________________________________________________

# make Species Richness df and remove na's
richness <- modeldata.df%>%
  select(-Site_Biomassm2,-Sites_Abundancem2,-logBiomass,-logAbundance)%>%
  drop_na() #only leaves 245 observations because soil data from SoilsGrid filled in a lot
hist(richness$SpeciesRichness) #left skewed
names(richness)
str(richness)

#_______________________________________________________
#Correlation Plot
#_______________________________________________________
#draw out possible relationship with all variables
names(richness)
#correlation plot (not including categorical variables: ESA, SnowMonths_cat)
cor.df <- richness %>%
  select(SpeciesRichness,
         Organic_Carbon__percent,
         C.N_ratio,
         Sand__percent,
         Silt__percent,
         Clay__percent,
         ph_new,
         bio10_1,
         bio10_4,
         bio10_7,
         bio10_12,
         bio10_15,
         elevation,
         Aridity,
         PETyr,
         PET_SD
  )

cor_matrix <- cor(cor.df, use = "pairwise.complete.obs")

corrplot.mixed(cor_matrix, lower.col = "black", upper = "ellipse", tl.col = "black", 
               number.cex = 0.7, tl.pos = "d", tl.cex = .55, 
               p.mat = cor.mtest(cor.df, conf.level = 0.95, use = "pairwise.complete.obs")$p, 
               sig.level = 0.05)

cor_long <- as.data.frame(as.table(cor_matrix))
colnames(cor_long) <- c("var1", "var2", "correlation")
cor_long <- cor_long[upper.tri(cor_matrix), ]
cor_long_filtered <- cor_long %>%
  filter(abs(correlation) < 0.7)

cor_long_filtered <- cor_long_filtered %>%
  mutate(var1 = as.character(var1), var2 = as.character(var2)) %>%
  mutate(pair_id = pmin(var1, var2)) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

cortable <- as.data.frame(cor(cor.df)) #or df table to look at numbers

#now iterative remove high VIF scores until everything is under the 3 threshold
remove_high_vif <- function(data, threshold = 3, response_var = "SpeciesRichness") {
  predictors <- setdiff(names(data), response_var)
  while(TRUE) {
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    model <- lm(formula, data = data)
    vif_values <- vif(model)
    if(all(vif_values < threshold)) {
      break  # Exit the loop if all VIFs are below threshold
    }
    
    max_vif_var <- names(which.max(vif_values))     # Find the variable with the highest VIF
    predictors <- setdiff(predictors, max_vif_var)  # Remove the variable with the highest VIF
    
    cat("Removed variable:", max_vif_var, "with VIF:", max(vif_values), "\n")
  }
  
  return(data[, c(response_var, predictors)]) # Return the filtered dataset
}
str(richness)
filtered_data <- richness%>% #need to remove any categories and non-relevant data
  select(-Study_Name, -Site_Name, - SnowMonths_cat, -SnowMonths)
#run the function
filtered_data <- remove_high_vif(filtered_data, threshold = 3)
print(colnames(filtered_data)) #these are the variables to stay

#VIF final variables from Phillips et al # bio10_7,bio10_15,CECSOL,elevation,Aridity,PETyr,phFinal,ClayFinal,SiltFinal,OCFinal
#VIF final variables from me # Organic_Carbon__percent, C.N_ratio, Silt__percent, Clay__percent, ph_new, bio10_12, elevation

#scale all continuous variables
vars_to_scale <- c("Organic_Carbon__percent", "C.N_ratio", "Silt__percent", "Clay__percent", "ph_new", "bio10_12", "elevation")
richness <- richness %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(richness))
#_______________________________________________________
#Tets models on final variables
#_______________________________________________________
#test model before curve
names(richness)
model <- glmer(SpeciesRichness ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                 Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                 elevation_scaled + (1|Study_Name),
             data = richness, family = poisson)
vif(model) # all low VIF
#now add back in categorical models and check VIF again
model <- glmer(SpeciesRichness ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                 Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                 elevation_scaled + SnowMonths_cat + (1|Study_Name),
             data = richness, family = poisson)
vif(model) #all still under 3

#_______________________________________________________
#Rich. RealDat Curve
#_______________________________________________________
str(richness)
richness$SnowMonths_cat <- as.factor(richness$SnowMonths_cat)

#Richness should be poisson because of left skewed count but it doesn't work because it's havin trouble with cat.
lm_poisson <- function(formula, data) {
  glm(formula = formula, 
      data = data, 
      family = poisson(link = "identity"))
}

#test organic carbon
vars <- c("Organic_Carbon__percent_scaled") # consistently positive
covariates <- c("C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")

specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_OC_Rich <- plot(results, type = "curve")+labs(x = "Org. C")
summary(results, type = "curve") #extract the median 
#plot(results)

#test C:N
vars <- c("C.N_ratio_scaled") # mostly not significant but some positive
covariates <- c("Organic_Carbon__percent_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")

specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)

R_C.N_Rich <- plot(results, type = "curve")+labs(x = "C:N Ratio")
summary(results, type = "curve") #extract the median 

#test silt 
vars <- c("Silt__percent_scaled") # overwhelmingly positive
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")

specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_Silt_Rich <- plot(results, type = "curve")+labs(x = "% Silt")
summary(results, type = "curve") #extract the median 


#test clay
vars <- c("Clay__percent_scaled") # overwhelmingly negative
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled","Silt__percent_scaled",
                 "ph_new_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")

specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_Clay_Rich <- plot(results, type = "curve")+labs(x = "% Clay")
summary(results, type = "curve") #extract the median 

#test ph
vars <- c("ph_new_scaled")# highly variable on covariates (many are not significant)
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled","bio10_12_scaled",
                "Clay__percent_scaled", "elevation_scaled", "SnowMonths_cat")
specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_ph_Rich <- plot(results, type = "curve")+labs(x = "pH")
summary(results, type = "curve") #extract the median 


#test elevation 
vars <- c("elevation_scaled")# elevation consistently neg
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled","bio10_12_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat")
specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_elevation_Rich <- plot(results, type = "curve")+labs(x = "Elevation")
summary(results, type = "curve") #extract the median 

#test bio 12 
vars <- c("bio10_12_scaled") #consistently positive
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled","elevation_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat")
specs <- setup(data = richness,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
R_bio12_Rich <- plot(results, type = "curve")+labs(x = "Precip")
summary(results, type = "curve") #extract the median 

#Rich. RealDat Fig
pRichReal<-plot_grid(R_OC_Rich,R_Silt_Rich,R_Clay_Rich,R_ph_Rich,R_elevation_Rich,R_bio12_Rich,R_C.N_Rich,
                    align = "v",
                    axis = "rbl",
                    ncol = 1)

if(TRUE){
  pdf("Richness-realDat.pdf", width = 3, height = 9)
  plot(pRichReal)
  dev.off()}


#_______________________________________________________#_______________________________________________________
# 6) ABUNDANCE----
#_______________________________________________________#_______________________________________________________


# make Species abundance df and remove na's
abundance <- modeldata.df%>%
  select(-Site_Biomassm2,-Sites_Abundancem2,-logBiomass, -SpeciesRichness)%>%
  drop_na() #only leaves 289 observations because soil data from SoilsGrid filled in a lot
hist(abundance$logAbundance) #left skewed
names(abundance)
str(abundance)

#_______________________________________________________
#Correlation Plot
#_______________________________________________________
#draw out possible relationship with all variables
names(abundance)
#correlation plot (not including categorical variables: ESA, SnowMonths_cat)
cor.df <- abundance %>%
  select(logAbundance,
         Organic_Carbon__percent,
         C.N_ratio,
         Sand__percent,
         Silt__percent,
         Clay__percent,
         ph_new,
         bio10_1,
         bio10_4,
         bio10_7,
         bio10_12,
         bio10_15,
         elevation,
         Aridity,
         PETyr,
         PET_SD
  )

cor_matrix <- cor(cor.df, use = "pairwise.complete.obs")

corrplot.mixed(cor_matrix, lower.col = "black", upper = "ellipse", tl.col = "black", 
               number.cex = 0.7, tl.pos = "d", tl.cex = .55, 
               p.mat = cor.mtest(cor.df, conf.level = 0.95, use = "pairwise.complete.obs")$p, 
               sig.level = 0.05)

cor_long <- as.data.frame(as.table(cor_matrix))
colnames(cor_long) <- c("var1", "var2", "correlation")
cor_long <- cor_long[upper.tri(cor_matrix), ]
cor_long_filtered <- cor_long %>%
  filter(abs(correlation) < 0.7)

cor_long_filtered <- cor_long_filtered %>%
  mutate(var1 = as.character(var1), var2 = as.character(var2)) %>%
  mutate(pair_id = pmin(var1, var2)) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

cortable <- as.data.frame(cor(cor.df)) #or df table

#now iterative remove high VIF scores until everything is under the 3 threshold
remove_high_vif <- function(data, threshold = 3, response_var = "logAbundance") {
  predictors <- setdiff(names(data), response_var)
  while(TRUE) {
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    model <- lm(formula, data = data)
    vif_values <- vif(model)
    if(all(vif_values < threshold)) {
      break  # Exit the loop if all VIFs are below threshold
    }
    
    max_vif_var <- names(which.max(vif_values))     # Find the variable with the highest VIF
    predictors <- setdiff(predictors, max_vif_var)  # Remove the variable with the highest VIF
    
    cat("Removed variable:", max_vif_var, "with VIF:", max(vif_values), "\n")
  }
  
  return(data[, c(response_var, predictors)]) # Return the filtered dataset
}
str(abundance)
filtered_data <- abundance%>% #need to remove any categories and non-relevant data
  select(-Study_Name, -Site_Name, - SnowMonths_cat, -SnowMonths)
#run the function
filtered_data <- remove_high_vif(filtered_data, threshold = 3)
print(colnames(filtered_data)) #these are the variables to stay

#VIF final variables from Phillips et al # bio10_7,bio10_15,CECSOL,elevation,Aridity,PETyr,phFinal,ClayFinal,SiltFinal,OCFinal
#VIF final variables from me # Organic_Carbon__percent, C.N_ratio, Silt__percent, Clay__percent, ph_new, bio10_12, elevation

#scale all continuous variables
vars_to_scale <- c("Organic_Carbon__percent", "C.N_ratio", "Silt__percent", "Clay__percent", "ph_new", "bio10_12", "elevation")
abundance <- abundance %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(abundance))
#_______________________________________________________
#Tets models on final variables
#_______________________________________________________
#test model before curve
sum(is.na(abundance))
names(abundance)
model <- glmer(logAbundance ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                 Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                 elevation_scaled + (1|Study_Name),
               data = abundance, family = gaussian)
vif(model) # all low VIF
#now add back in categorical models and check VIF again
model <- lmer(logAbundance ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                 Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                 elevation_scaled  + SnowMonths_cat + (1|Study_Name),
               data = abundance)
vif(model) #nothing to remove
#_______________________________________________________
#Abun. RealDat Curve----
#_______________________________________________________
str(abundance)
abundance$SnowMonths_cat <- as.factor(abundance$SnowMonths_cat)

#abundance should be poisson because of left skewed count but it doesn't work
lm_gaussian <- function(formula, data) {
  glm(formula = formula, 
      data = data, 
      family = gaussian(link = "identity"))
}

#test each variable as the predictor

#test organic carbon
vars <- c("Organic_Carbon__percent_scaled") # consistently neg
covariates <- c("C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_OC_Abun <- plot(results, type = "curve")+labs(x = "Org. C")
summary(results, type = "curve") #extract the median 

#test C:N
vars <- c("C.N_ratio_scaled") # not significant but a few neg
covariates <- c("Organic_Carbon__percent_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_C.N_Abun <- plot(results, type = "curve") + labs(x = "C:N Ratio")
summary(results, type = "curve") #extract the median 

#test silt 
vars <- c("Silt__percent_scaled") # mostly positive
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "elevation_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_Silt_Abun <- plot(results, type = "curve")+labs(x = "% Silt")
summary(results, type = "curve") #extract the median 

#test clay
vars <- c("Clay__percent_scaled") # mostly negative
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled","Silt__percent_scaled",
                "ph_new_scaled", "elevation_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_Clay_Abun <- plot(results, type = "curve")+labs(x = "% Clay")
summary(results, type = "curve") #extract the median 

#test ph
vars <- c("ph_new_scaled")# mostly negative (many are no significant)
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "elevation_scaled", "SnowMonths_cat", "bio10_12_scaled")
specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_ph_Abun <- plot(results, type = "curve")+labs(x = "pH")
summary(results, type = "curve") #extract the median 

#test elevation 
vars <- c("elevation_scaled")# elevation consistently neg
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat", "bio10_12_scaled")
specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_elevation_Abun <- plot(results, type = "curve")+labs(x = "Elevation")
summary(results, type = "curve") #extract the median 

#test precip 
vars <- c("bio10_12_scaled")# elevation consistently neg
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat","elevation_scaled")
specs <- setup(data = abundance,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_bio12_Abun <- plot(results, type = "curve")+labs(x = "Precip")
summary(results, type = "curve") #extract the median 

#Abun. RealDat Fig
pAbunReal<-plot_grid(R_OC_Abun,R_Silt_Abun,R_Clay_Abun,R_ph_Abun,R_elevation_Abun,R_bio12_Abun,R_C.N_Rich,
                     align = "v",
                     axis = "rbl",
                     ncol = 1)

if(TRUE){
  pdf("Abundance-realDat.pdf", width = 3, height = 9)
  plot(pAbunReal)
  dev.off()}



#_______________________________________________________#_______________________________________________________
#7) BIOMASS----
#_______________________________________________________#_______________________________________________________

# make Biomass df and remove na's
biomass <- modeldata.df%>%
  select(-Site_Biomassm2,-Sites_Abundancem2,-logAbundance, -SpeciesRichness)%>%
  drop_na() #only leaves 208 observations because soil data from SoilsGrid filled in a lot
hist(biomass$logBiomass) #left skewed
names(biomass)
str(biomass)

#_______________________________________________________
#Correlation Plot
#_______________________________________________________
#draw out possible relationship with all variables
names(biomass)
#correlation plot (not including categorical variables: ESA, SnowMonths_cat)
cor.df <- biomass %>%
  select(logBiomass,
         Organic_Carbon__percent,
         C.N_ratio,
         Sand__percent,
         Silt__percent,
         Clay__percent,
         ph_new,
         bio10_1,
         bio10_4,
         bio10_7,
         bio10_12,
         bio10_15,
         elevation,
         Aridity,
         PETyr,
         PET_SD
  )

cor_matrix <- cor(cor.df, use = "pairwise.complete.obs")

corrplot.mixed(cor_matrix, lower.col = "black", upper = "ellipse", tl.col = "black", 
               number.cex = 0.7, tl.pos = "d", tl.cex = .55, 
               p.mat = cor.mtest(cor.df, conf.level = 0.95, use = "pairwise.complete.obs")$p, 
               sig.level = 0.05)

cor_long <- as.data.frame(as.table(cor_matrix))
colnames(cor_long) <- c("var1", "var2", "correlation")
cor_long <- cor_long[upper.tri(cor_matrix), ]
cor_long_filtered <- cor_long %>%
  filter(abs(correlation) < 0.7)

cor_long_filtered <- cor_long_filtered %>%
  mutate(var1 = as.character(var1), var2 = as.character(var2)) %>%
  mutate(pair_id = pmin(var1, var2)) %>%
  distinct(pair_id, .keep_all = TRUE) %>%
  select(-pair_id)

cortable <- as.data.frame(cor(cor.df)) #or df table

#now iterative remove high VIF scores until everything is under the 3 threshold
remove_high_vif <- function(data, threshold = 3, response_var = "logBiomass") {
  predictors <- setdiff(names(data), response_var)
  while(TRUE) {
    formula <- as.formula(paste(response_var, "~", paste(predictors, collapse = " + ")))
    model <- lm(formula, data = data)
    vif_values <- vif(model)
    if(all(vif_values < threshold)) {
      break  # Exit the loop if all VIFs are below threshold
    }
    
    max_vif_var <- names(which.max(vif_values))     # Find the variable with the highest VIF
    predictors <- setdiff(predictors, max_vif_var)  # Remove the variable with the highest VIF
    
    cat("Removed variable:", max_vif_var, "with VIF:", max(vif_values), "\n")
  }
  
  return(data[, c(response_var, predictors)]) # Return the filtered dataset
}
str(biomass)
filtered_data <- biomass%>% #need to remove any categories and non-relevant data
  select(-Study_Name, -Site_Name, - SnowMonths_cat, -SnowMonths)
#run the function
filtered_data <- remove_high_vif(filtered_data, threshold = 3)
print(colnames(filtered_data)) #these are the variables to stay

#VIF final variables from Phillips et al # bio10_7,bio10_15,CECSOL,elevation,Aridity,PETyr,phFinal,ClayFinal,SiltFinal,OCFinal
#VIF final variables from me # Organic_Carbon__percent, C.N_ratio, Silt__percent, Clay__percent, ph_new, bio10_12

#scale all continuous variables
vars_to_scale <- c("Organic_Carbon__percent", "C.N_ratio", "Silt__percent", "Clay__percent", "ph_new", "bio10_12")
biomass <- biomass %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(biomass))
#_______________________________________________________
#Tets models on final variables
#_______________________________________________________
#test model before curve
sum(is.na(biomass))
names(biomass)
model <- glmer(logBiomass ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                 Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                 (1|Study_Name),
               data = biomass, family = gaussian)
vif(model) # all low VIF
#now add back in categorical models and check VIF again
model <- lmer(logBiomass ~ Organic_Carbon__percent_scaled + C.N_ratio_scaled + 
                Silt__percent_scaled + Clay__percent_scaled + ph_new_scaled + bio10_12_scaled + 
                SnowMonths_cat + (1|Study_Name),
              data = biomass)
vif(model) #nothing to remove
#_______________________________________________________
#Bio. RealDat Curve
#_______________________________________________________
str(biomass)
biomass$SnowMonths_cat <- as.factor(biomass$SnowMonths_cat)

#test organic carbon
vars <- c("Organic_Carbon__percent_scaled") # mostly neg, some ns 
covariates <- c("C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_OC_Bio <- plot(results, type = "curve")+labs(x = "Org. C")
summary(results, type = "curve") #extract the median 

#test C:N
vars <- c("C.N_ratio_scaled") # not significant 
covariates <- c("Organic_Carbon__percent_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_C.N_Bio <- plot(results, type = "curve")+labs(x = "C:N Ratio")
summary(results, type = "curve") #extract the median 

#test silt 
vars <- c("Silt__percent_scaled") # overwhelmingly positive
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled",
                "Clay__percent_scaled", "ph_new_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_Silt_Bio <- plot(results, type = "curve")+labs(x = "% Silt")
summary(results, type = "curve") #extract the median 

#test clay
vars <- c("Clay__percent_scaled") # not sig
covariates <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled","Silt__percent_scaled",
                "ph_new_scaled", "SnowMonths_cat", "bio10_12_scaled")

specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_Clay_Bio <- plot(results, type = "curve")+labs(x = "% Clay")
summary(results, type = "curve") #extract the median 

#test ph
vars <- c("ph_new_scaled")# not significant
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "SnowMonths_cat", "bio10_12_scaled")
specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_ph_Bio <- plot(results, type = "curve")+labs(x = "pH")
summary(results, type = "curve") #extract the median 

#test ph
vars <- c("bio10_12_scaled")# not significant
covariates <- c("Organic_Carbon__percent_scaled","C.N_ratio_scaled", "Silt__percent_scaled",
                "Clay__percent_scaled", "SnowMonths_cat", "ph_new_scaled")
specs <- setup(data = biomass,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm_gaussian"))

results <- specr(specs)
R_bio12_Bio <- plot(results, type = "curve")+labs(x = "Precip")
summary(results, type = "curve") #extract the median 

# Bio. RealDat Fig.
pBioReal<-plot_grid(R_OC_Bio, R_Silt_Bio, R_Clay_Bio, R_ph_Bio, R_bio12_Bio, R_C.N_Bio, 
             align = "v",
             axis = "rbl",
             ncol = 1)

if(TRUE){
  pdf("Biomass-realDat.pdf", width = 3, height = 9)
  plot(pBioReal)
  dev.off()}



#make fig all real dat

realdatFIG <- plot_grid(
  plot_grid(R_OC_Rich, R_C.N_Rich, R_Silt_Rich, R_Clay_Rich, R_ph_Rich,R_bio12_Rich, R_elevation_Rich, 
            ncol = 1, labels = c("a", "b", "c", "d", "e", "f", "g")),
  plot_grid(R_OC_Abun, R_C.N_Rich, R_Silt_Abun, R_Clay_Abun, R_ph_Abun, R_bio12_Abun,  R_elevation_Abun, 
            ncol = 1, labels = c("h", "i", "j", "k", "l", "m", "n")),
  plot_grid(R_OC_Bio,R_C.N_Bio, R_Silt_Bio, R_Clay_Bio, R_ph_Bio,  R_bio12_Bio, NULL,
            ncol = 1, labels = c("o", "p", "q", "r", "", "s", "")), 
  ncol = 3,
  align = "hv",
  axis = "tblr",
  label_size = 8
)

realdatFIG

if(TRUE){
  pdf("Fig-2.pdf", width = 6, height = 10)
  plot(realdatFIG)
  dev.off()}


#8) Regional Test----
coords%>%select(Country)%>%distinct()
coords <- coords%>%
  select(Site_Name, Study_Name, Latitude_decimal_degrees, Longitude_decimal_degrees, Country)%>%
  mutate(Study_Name_Site_Name = paste(Study_Name, Site_Name, sep = "_"))


# Add 'Region' column based on 'Country'
coords <- coords %>%
  mutate(Region = case_when(
    Country %in% c("United States", "Canada", "Mexico") ~ "North America",
    Country %in% c("Argentina", "Brazil", "Bolivia", "Peru", "Colombia", "Ecuador", "Uruguay", "Venezuela") ~ "South America",
    Country %in% c("Russia") ~ "Eastern Europe & Asia",
    Country %in% c("United Kingdom", "Ireland", "Germany", "France", "Poland", "Romania", "Hungary", "Austria", 
                   "Netherlands", "Belgium", "Sweden", "Finland", "Denmark", "Switzerland", "Estonia", "Lithuania",
                   "Slovakia", "Czech Republic", "Slovenia", "Croatia") ~ "Europe",
    Country %in% c("Spain", "Portugal", "Italy", "Greece") ~ "Southern Europe",
    Country %in% c("Australia", "New Zealand") ~ "Oceania",
    Country %in% c("China", "India", "Iran", "Japan", "Malaysia", "Philippines") ~ "Asia",
    Country %in% c("Kenya", "Nigeria", "Ghana", "Benin", "Niger", "Cote d'Ivoire", "Cameroon", "Burkina Faso", 
                   "Malawi", "Libya", "Congo", "Madagascar") ~ "Africa",
    Country %in% c("Hawaii", "Puerto Rico") ~ "Islands",
    Country == "Border" ~ "Unknown",
    TRUE ~ "Other"  # For any country not listed explicitly
  ))


richness_region <- richness%>%
  mutate(Study_Name_Site_Name = paste(Study_Name, Site_Name, sep = "_"))%>%
  left_join(coords, by = "Study_Name_Site_Name")


abundance_region <- abundance%>%
  mutate(Study_Name_Site_Name = paste(Study_Name, Site_Name, sep = "_"))%>%
  left_join(coords, by = "Study_Name_Site_Name")

biomass_region <- biomass%>%
  mutate(Study_Name_Site_Name = paste(Study_Name, Site_Name, sep = "_"))%>%
  left_join(coords, by = "Study_Name_Site_Name")

biomass_region%>%select(Country)%>%distinct()

#Map
# Get world map data
world <- ne_countries(scale = "medium", returnclass = "sf")
map <- ggplot() +
  geom_sf(data = world, fill = "#bae4bc", color = "lightgrey") +
  geom_point(data = abundance_region, aes(x = Longitude_decimal_degrees, y = Latitude_decimal_degrees),
             color = "darkgrey", size = 2) +
  coord_sf(xlim = c(-90, 60), ylim = c(0, 70), expand = TRUE) + # Adjust limits as needed
  labs(title = "",
       x = "Longitude", y = "Latitude") +
  theme_minimal()

map

if(TRUE){
  pdf("Fig-1.pdf", width = 4, height = 3)
  plot(map)
  dev.off()}

#9) Random Forest Test----

# Function to calculate variable share (relative importance)
var.share <- function(rf.obj, members) {
  count <- table(rf.obj$forest$bestvar)[-1]  # Count the number of times each variable is used
  names(count) <- names(rf.obj$forest$ncat)  # Assign names to the count table
  share <- count[members] / sum(count[members])  # Calculate relative share (importance)
  return(share)
}

# Function to calculate group importance
group.importance <- function(rf.obj, groups) {
  var.imp <- as.matrix(sapply(groups, function(g) {
    # Sum the individual variable importance for each group, weighted by var.share
    sum(importance(rf.obj, 2)[g, ] * var.share(rf.obj, g))
  }))
  colnames(var.imp) <- "MeanDecreaseGini"  # Label the output as 'MeanDecreaseGini'
  return(var.imp)
}

# Function to order the group importances (from highest to lowest)
OrderImportance <- function(groupedImportances){
  order <- order(groupedImportances, decreasing = TRUE)  # Sort in descending order
  return(order)
}

#Richness Random Forest
# Define your main effects variables for richness model
RichmainEffects <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
                      "Clay__percent_scaled", "ph_new_scaled", "bio10_12_scaled", 
                      "elevation_scaled", "SnowMonths_cat")

# Fit the Random Forest model for richness
richness_rf <- randomForest(y = richness$SpeciesRichness, 
                             x = richness[,names(richness) %in% RichmainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Elevation <- "elevation_scaled"
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
          "Clay__percent_scaled", "ph_new_scaled")

# Create the groups list
groups <- list(
  Elevation = Elevation, 
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
richness_import_split <- group.importance(richness_rf, groups)

# Print the grouped importance
print(richness_import_split)

# Optionally, order the importance from highest to lowest
ordered_importance <- OrderImportance(richness_import_split)
print(ordered_importance)

#Abundance RandomForest
# Define your main effects variables for abundance model
AbundmainEffects <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
                      "Clay__percent_scaled", "ph_new_scaled", "bio10_12_scaled", 
                      "elevation_scaled", "SnowMonths_cat")

# Fit the Random Forest model for abundance
abundance_rf <- randomForest(y = abundance$logAbundance, 
                             x = abundance[,names(abundance) %in% AbundmainEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Elevation <- "elevation_scaled"
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
          "Clay__percent_scaled", "ph_new_scaled")

# Create the groups list
groups <- list(
  Elevation = Elevation, 
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
abundance_import_split <- group.importance(abundance_rf, groups)

# Print the grouped importance
print(abundance_import_split)

# Optionally, order the importance from highest to lowest
ordered_importance <- OrderImportance(abundance_import_split)
print(ordered_importance)


#Biomass Random Forest
# Define your main effects variables for biomass model
BiomainEffects <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
                     "Clay__percent_scaled", "ph_new_scaled", "bio10_12_scaled", 
                      "SnowMonths_cat")
names(biomass)
# Fit the Random Forest model for biomass
biomass_rf <- randomForest(y = biomass$logBiomass, 
                            x = biomass[,names(biomass) %in% BiomainEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("Organic_Carbon__percent_scaled", "C.N_ratio_scaled", "Silt__percent_scaled", 
          "Clay__percent_scaled", "ph_new_scaled")

# Create the groups list
groups <- list(
  Elevation = NA,
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
biomass_import_split <- group.importance(biomass_rf, groups)

# Print the grouped importance
print(biomass_import_split)

ordered_importance <- OrderImportance(biomass_import_split)
print(ordered_importance)

#Fig Random Forest

# Combine the results into a table
group_importance_table <- data.frame(
  Model = rep(c("Richness", "Abundance", "Biomass"), each = length(groups)),
  Group = rep(names(groups), times = 3),
  MeanDecreaseGini = c(richness_import_split, abundance_import_split, biomass_import_split)
)

# Create the bubble plot
rf <- ggplot(group_importance_table, aes(x = reorder(Group, MeanDecreaseGini), y = factor(Model, levels = c("Biomass","Abundance", "Richness")), size = MeanDecreaseGini)) +
  geom_point(alpha = 0.7) +  # Set transparency to make bubbles more visible
  scale_size_continuous(range = c(3, 15)) +  # Adjust the range of bubble sizes
  labs(title = "",
       x = "", 
       y = "Model", 
       size = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")

if(TRUE){
  pdf("Fig-3.pdf", width = 5, height = 3)
  plot(rf)
  dev.off()}

#_______________________________________________________
#Appendix SoilsGrid----
#_______________________________________________________
#CHecking if the same predictors, when applied to siols grid data, mean anything??
#don't have time to include

# a little cleaning
richnessGrid <- origData%>%
  select(Study_Name, Site_Name, OCFinal, SiltFinal, ClayFinal, phFinal, elevation, SnowMonths_cat, bio10_12, SpeciesRichness)

vars_to_scale <- c("OCFinal", "SiltFinal", "ClayFinal", "phFinal", "elevation","bio10_12")

richnessGrid <- richnessGrid %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(richnessGrid))

richnessGrid$SnowMonths_cat <- as.factor(richnessGrid$SnowMonths_cat)

sum(is.na(richnessGrid))
columns_to_check <- c("SpeciesRichness", "bio10_12", "SiltFinal_scaled", "OCFinal_scaled", 
                      "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat")

richnessGrid <- richnessGrid %>%
  filter(if_all(all_of(columns_to_check), ~!is.na(.)))
sum(is.na(richnessGrid))

#organic carbon
vars <- c("OCFinal_scaled") # soil -OC 1/2 positive, 1/2 ns
covariates <- c("SiltFinal_scaled","bio10_12",
                "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat")
specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_OC_Rich <- plot(results, type = "curve")+labs(x = "Org. C")

#SiltFinal_scaled
vars <- c("SiltFinal_scaled") # silt 1/2 positive 1/2 neg
covariates <- c("bio10_12","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat")
specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Silt_Rich <- plot(results, type = "curve")+labs(x = "% Silt")

#ClayFinal_scaled
vars <- c("ClayFinal_scaled") # clay mostly negative, some positive
covariates <- c("SiltFinal_scaled","bio10_12","OCFinal_scaled",
                "phFinal_scaled", "elevation_scaled", "SnowMonths_cat")
specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Clay_Rich <- plot(results, type = "curve")+labs(x = "% Clay")


#ph
vars <- c("phFinal_scaled") # consistently positive
covariates <- c("SiltFinal_scaled","bio10_12","OCFinal_scaled",
                "ClayFinal_scaled",  "SnowMonths_cat","elevation_scaled")

specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_ph_Rich <- plot(results, type = "curve")+labs(x = "pH")

#elevation_scaled
vars <- c("elevation_scaled")# highly negative 
covariates <- c("SiltFinal_scaled","bio10_12","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat")
specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_elevation_Rich <- plot(results, type = "curve")+labs(x = "Elevation")

#bio 12
vars <- c("bio10_12")# highly negative 
covariates <- c("SiltFinal_scaled","elevation_scaled","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat")
specs <- setup(data = richnessGrid,
               x = vars,
               y = c("SpeciesRichness"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_bio12_Rich <- plot(results, type = "curve")+labs(x = "Precip")

#Rich. Grid Fig---
pRichGrid<-plot_grid(G_OC_Rich,G_Silt_Rich,G_Clay_Rich,G_ph_Rich,G_elevation_Rich,G_bio12_Rich,
                     align = "v",
                     axis = "rbl",
                     ncol = 1)

if(TRUE){
  pdf("Richness-Grid.pdf", width = 3, height = 9)
  plot(pRichGrid)
  dev.off()}

richnessFIG <- plot_grid(
  plot_grid(R_OC_Rich, R_Silt_Rich, R_Clay_Rich, R_ph_Rich, R_elevation_Rich, R_bio12_Rich, R_C.N_Rich, ncol = 1),
  plot_grid(G_OC_Rich, G_Silt_Rich, G_Clay_Rich, G_ph_Rich, G_elevation_Rich, G_bio12_Rich, NULL, ncol = 1),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1.2),
  labels = c("a","b"),
  vjust = 1,
  label_size = 12
)


if(TRUE){
  pdf("Richness-Both.pdf", width = 4.5, height = 9)
  plot(richnessFIG)
  dev.off()}

#_______________________________________________________
#Abun. Grid Curve---
#_______________________________________________________
# a little cleaning
abundanceGrid <- origData%>%
  select(Study_Name, Site_Name, OCFinal, SiltFinal, ClayFinal, phFinal, elevation, SnowMonths_cat, bio10_12, logAbundance)

vars_to_scale <- c("OCFinal", "SiltFinal", "ClayFinal", "phFinal", "elevation", "bio10_12")
abundanceGrid <- abundanceGrid %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(abundanceGrid))

abundanceGrid$SnowMonths_cat <- as.factor(abundanceGrid$SnowMonths_cat)

sum(is.na(abundanceGrid))
columns_to_check <- c("logAbundance", "SiltFinal_scaled", "OCFinal_scaled", 
                      "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")

abundanceGrid <- abundanceGrid %>%
  filter(if_all(all_of(columns_to_check), ~!is.na(.)))
sum(is.na(abundanceGrid))# zero zero!

#check each predictor

#organic carbon
vars <- c("OCFinal_scaled") # soil -OC some negative mostly NS
covariates <- c("SiltFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_OC_Abun <- plot(results, type = "curve")+labs(x = "Org. C")

#SiltFinal_scaled
vars <- c("SiltFinal_scaled") # overwgelming neg
covariates <- c("OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Silt_Abun <- plot(results, type = "curve")+labs(x = "% Silt")

#ClayFinal_scaled
vars <- c("ClayFinal_scaled") # clay all negative 
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "phFinal_scaled", "elevation_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Clay_Abun <- plot(results, type = "curve")+labs(x = "% Clay")

#ph
vars <- c("phFinal_scaled") # mostly negative
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "ClayFinal_scaled",  "SnowMonths_cat","elevation_scaled","bio10_12_scaled")

specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_ph_Abun <- plot(results, type = "curve")+labs(x = "ph")

#elevation_scaled
vars <- c("elevation_scaled")# highly negative 
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_elevation_Abun <- plot(results, type = "curve")+labs(x = "Elevation")

#bio 12
vars <- c("bio10_12_scaled")# mostly positive 
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled","elevation_scaled")
specs <- setup(data = abundanceGrid,
               x = vars,
               y = c("logAbundance"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_bio12_Abun <- plot(results, type = "curve")+labs(x = "Precip")
#-----------------------------------------#
#Abun. Fig---
#-----------------------------------------#
abundanceFIG <- plot_grid(
  plot_grid(R_OC_Abun, R_Silt_Abun, R_Clay_Abun, R_ph_Abun, R_elevation_Abun, R_bio12_Abun, R_C.N_Abun, ncol = 1),
  plot_grid(G_OC_Abun, G_Silt_Abun, G_Clay_Abun, G_ph_Abun, G_elevation_Abun, G_bio12_Abun, NULL, ncol = 1),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1.2),
  labels = c("a","b"),
  vjust = 1,
  label_size = 12
)


if(TRUE){
  pdf("abundanceFIG-Both.pdf", width = 4.5, height = 9)
  plot(abundanceFIG)
  dev.off()}

#_______________________________________________________
#Bio. Grid Curve---
#_______________________________________________________
# a little cleaning
biomassGrid <- origData%>%
  select(Study_Name, Site_Name, OCFinal, SiltFinal, ClayFinal, phFinal, SnowMonths_cat, bio10_12, logBiomass)

vars_to_scale <- c("OCFinal", "SiltFinal", "ClayFinal", "phFinal", "bio10_12")
biomassGrid <- biomassGrid %>%
  mutate(across(all_of(vars_to_scale), ~ as.numeric(scale(.)), .names = "{.col}_scaled"))
print(colnames(biomassGrid))


biomassGrid$SnowMonths_cat <- as.factor(biomassGrid$SnowMonths_cat)

sum(is.na(biomassGrid))
columns_to_check <- c("logBiomass", "SiltFinal_scaled", "OCFinal_scaled", 
                      "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")

biomassGrid <- biomassGrid %>%
  filter(if_all(all_of(columns_to_check), ~!is.na(.)))
sum(is.na(biomassGrid))# zero zero!

#check each predictor

#organic carbon
vars <- c("OCFinal_scaled") # soil -OC positive 
covariates <- c("SiltFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = biomassGrid,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm"))


results <- specr(specs)
G_OC_Bio <- plot(results, type = "curve")+ labs(x = "Org. C")

#SiltFinal_scaled
vars <- c("SiltFinal_scaled") # overwhelmingly neg
covariates <- c("OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = biomassGrid,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Silt_Bio <- plot(results, type = "curve")+ labs(x = "% Silt")

#ClayFinal_scaled
vars <- c("ClayFinal_scaled") # clay 1/4 negative , 3/4 ns
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = biomassGrid,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_Clay_Bio <- plot(results, type = "curve")+ labs(x = "% Clay")

#ph
vars <- c("phFinal_scaled") # half positive half ns
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "ClayFinal_scaled",  "SnowMonths_cat","bio10_12_scaled")

specs <- setup(data = biomassGrid,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_ph_Bio <- plot(results, type = "curve")+ labs(x = "pH")

#bio 12
vars <- c("bio10_12_scaled")# NOT SIG
covariates <- c("SiltFinal_scaled","OCFinal_scaled",
                "ClayFinal_scaled", "phFinal_scaled", "SnowMonths_cat","bio10_12_scaled")
specs <- setup(data = biomassGrid,
               x = vars,
               y = c("logBiomass"),
               controls = covariates,
               model = c("lm"))

results <- specr(specs)
G_bio12_Bio <- plot(results, type = "curve")+ labs(x = "Precip")




#Bio. Grid Fig---
pBioGrid<-plot_grid(G_OC_Bio, G_Silt_Bio, G_Clay_Bio, G_ph_Bio, G_bio12_Bio,  
                    align = "v",
                    axis = "rbl",
                    ncol = 1)

if(TRUE){
  pdf("Biomass-SoilsGrid.pdf", width = 3, height = 9)
  plot(pBioGrid)
  dev.off()}

#Bio. Fig---

biomassFIG <- plot_grid(
  plot_grid(R_OC_Bio, R_Silt_Bio, R_Clay_Bio, R_ph_Bio, R_bio12_Bio, ncol = 1),
  plot_grid(G_OC_Bio, G_Silt_Bio, G_Clay_Bio, G_ph_Bio, G_bio12_Bio, ncol = 1),
  ncol = 2,
  align = "hv",
  axis = "tblr",
  rel_widths = c(1, 1.2),
  labels = c("a","b"),
  vjust = 1,
  label_size = 12
)


if(TRUE){
  pdf("biomassFIG-Both.pdf", width = 4.5, height = 9)
  plot(biomassFIG)
  dev.off()}


#make fig all soilsgrid
GriddatFIG <- plot_grid(
  plot_grid(G_OC_Rich,G_Silt_Rich,G_Clay_Rich,G_ph_Rich,G_elevation_Rich,G_bio12_Rich, ncol = 1),
  plot_grid(G_OC_Abun, G_Silt_Abun, G_Clay_Abun, G_ph_Abun, G_elevation_Abun, G_bio12_Abun, ncol = 1),
  plot_grid(G_OC_Bio, G_Silt_Bio, G_Clay_Bio, G_ph_Bio, NULL,G_bio12_Bio, ncol = 1),
  ncol = 3,
  align = "hv",
  axis = "tblr",
  #  rel_widths = c(1, 1.2),
  labels = c("a","b","c"),
  vjust = 1,
  label_size = 12
)
GriddatFIG
if(TRUE){
  pdf("GriddatFIG.pdf", width = 5, height = 8)
  plot(GriddatFIG)
  dev.off()}



#SoilGrid Random Forest---

#Abundance Random Forest
# Define your main effects variables for abundanceGrid model
names(abundanceGrid)
AbundmainGEffects <- c("OCFinal_scaled", "SiltFinal_scaled", "ClayFinal_scaled", 
                       "phFinal_scaled", "elevation_scaled", "bio10_12_scaled", 
                       "SnowMonths_cat")

# Fit the Random Forest model for abundanceGrid
abundanceGrid_rf <- randomForest(y = abundanceGrid$logAbundance, 
                                 x = abundanceGrid[,names(abundanceGrid) %in% AbundmainGEffects], 
                                 ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Elevation <- "elevation_scaled"
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("OCFinal_scaled", "SiltFinal_scaled", "ClayFinal_scaled", 
          "phFinal_scaled")

# Create the groups list
groups <- list(
  Elevation = Elevation, 
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
abundance_import_split <- group.importance(abundanceGrid_rf, groups)

# Print the grouped importance
print(abundance_import_split)

# Optionally, order the importance from highest to lowest
ordered_importance <- OrderImportance(abundance_import_split)
print(ordered_importance)

#Richness Random Forest
names(richnessGrid)
# Define your main effects variables for richness model
RichmainGEffects <- c("OCFinal_scaled",  "SiltFinal_scaled", 
                      "ClayFinal_scaled", "phFinal_scaled", "bio10_12_scaled", 
                      "elevation_scaled", "SnowMonths_cat")

# Fit the Random Forest model for richness
richnessG_rf <- randomForest(y = richnessGrid$SpeciesRichness, 
                             x = richnessGrid[,names(richnessGrid) %in% RichmainGEffects], 
                             ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Elevation <- "elevation_scaled"
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("OCFinal_scaled",  "SiltFinal_scaled", 
          "ClayFinal_scaled", "phFinal_scaled")

# Create the groups list
groups <- list(
  Elevation = Elevation, 
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
richness_import_split <- group.importance(richnessG_rf, groups)

# Print the grouped importance
print(richness_import_split)

# Optionally, order the importance from highest to lowest
ordered_importance <- OrderImportance(richness_import_split)
print(ordered_importance)


#Biomass Random Forest
# Define your main effects variables for biomass model
names(biomassGrid)
BiomainGEffects <- c("OCFinal_scaled", "SiltFinal_scaled", "ClayFinal_scaled", 
                     "phFinal_scaled", "bio10_12_scaled", 
                     "SnowMonths_cat")

# Fit the Random Forest model for biomass
biomassG_rf <- randomForest(y = biomassGrid$logBiomass, 
                            x = biomassGrid[,names(biomassGrid) %in% BiomainGEffects], 
                            ntree=501, importance=TRUE, proximity = TRUE)

# Define your groups
Precip <- c("bio10_12_scaled", "SnowMonths_cat")
Soil <- c("OCFinal_scaled", "SiltFinal_scaled", "ClayFinal_scaled")

# Create the groups list
groups <- list(
  Elevation = NA,
  Precip = Precip,
  Soil = Soil
)

# Calculate the group importance for the abundance Random Forest model
biomass_import_split <- group.importance(biomassG_rf, groups)

# Print the grouped importance
print(biomass_import_split)

ordered_importance <- OrderImportance(biomass_import_split)
print(ordered_importance)

#Random Forest Fig

# Combine the results into a table
group_importance_table <- data.frame(
  Model = rep(c("Richness", "Abundance", "Biomass"), each = length(groups)),
  Group = rep(names(groups), times = 3),
  MeanDecreaseGini = c(richness_import_split, abundance_import_split, biomass_import_split)
)

# Create the bubble plot
rf <- ggplot(group_importance_table, aes(x = reorder(Group, MeanDecreaseGini), y = factor(Model, levels = c("Biomass","Abundance", "Richness")), size = MeanDecreaseGini)) +
  geom_point(alpha = 0.7) +  # Set transparency to make bubbles more visible
  scale_size_continuous(range = c(3, 15)) +  # Adjust the range of bubble sizes
  labs(title = "",
       x = "", 
       y = "Model", 
       size = "") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  theme(legend.position = "none")
rf
#end---


#-------------------END---------------------------------------

