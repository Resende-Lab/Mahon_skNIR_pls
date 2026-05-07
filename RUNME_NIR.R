############################################################################################
###        Prediction of phytoglycogen content in sweet corn using 
###             Single-kernel Near-Infrared spectroscopy 
###
###   Mahon et al. 2024
###   doi: 
### 
############################################################################################

  rm(list = ls())


###-----------------------------------
####---- 1. Packages
###-----------------------------------

library(dplyr)
library(pls)


###-----------------------------------
####---- 2. PLS_stats Function
###-----------------------------------

PLS_stats <- function(trait, ncomp, k_fold = 10, seg = 10) {
  cat(trait, "\n")
  
  #----------- 1. Functions
  # SNV function
  SNV <- function(nir) { 
    nir <- as.matrix(nir)
    NIRsnv <- t(scale(t(nir), center = TRUE, scale = TRUE))
    return(NIRsnv)
  }   
 
   # RPD function
  RPD <- function(measured, RMSEP) {
    sd(measured)/RMSEP
  }
  
  #----------- 2. Dataset loading
  
  # dataset
  if (trait == "weight_mg") {
    dataset = "w"
  } else if (trait == "PG_PCT") {
    dataset = "pg"
  } else if (trait %in% c("starch_PCT", "starch_mg", "total_sugars_mg", "glucose_PCT", "sucrose_mg", "total_carbohydrates_mg")) {
    dataset = "s"
  }
  
  
  #----------- 3. Spectra pre-treatment
  
  if (trait %in% c("starch_PCT","glucose_PCT", "PG_PCT")) {
    spec = "snv"
  } else {
    spec = "raw"
  }
  
  #----------- 4 read calibration dataset 
  # weight
  if (dataset == "w") {
    d <- read.csv("weight_WL.NIR.csv") 
    d = d[d$cv_weight < 0.10, ]
    
    # Phytoglycogen
  } else if (dataset == "pg") {
    d <- read.csv("PG_WL.NIR.csv") 
    
    # sugars, starch
  } else if (dataset == "s") {
    d <- read.csv("sugars_WL.NIR.csv")
    if (trait != "total_carbohydrates_mg") {
      cv <- paste0("cv_", gsub("_(mg|PCT)$", "", trait))
      d <- d %>%  filter(d[, cv] < 0.11)
    }
  }
  d <- d[!(is.na(d[,trait])), ]
  
  #----------- 5 data.frame with genotype + measured + NIR  
  
  NIR <- d[, c(paste0("X", 940:1640))] %>%  as.matrix()
  data <- data.frame(d[,c("genotype",trait)], I(NIR)) 

  #----------- 6 Create objects to hold results
  
  output = list()
  model_stats = c()                   
  
  #----------- 7. Model
  
  # 10 iterations
  for (iter in 1:10) {
    cat("iter", iter, "\n")
    set.seed(iter)
    dat <- data %>%  mutate(set = sample(1:k_fold, nrow(data), replace = TRUE))
    
    # create an object for predicted values.
    Res = c()
    
    for (fold in 1:k_fold) {
      # split dataset into test and train
      test = dat %>% filter(set == fold)
      train = dat %>% filter(set != fold)
      
      # build train.pls model for each fold
      if (spec == "raw") { # raw NIR
        train.pls <- plsr(get(trait) ~ NIR, ncomp = 25, data = train, validation = "CV", segments = seg)
      } else if (spec == "snv") { # SNV(NIR)
        train.pls <- plsr(get(trait) ~ SNV(NIR), ncomp = 25, data = train, validation = "CV", segments = seg)
      }
      
      # predict test data (outer-validation) with train.pls
      t <- test %>% 
        select(genotype, dplyr::all_of(trait)) %>% 
        mutate(predicted = predict(train.pls, ncomp = ncomp, newdata = test) )
      # data
      Res <- rbind(Res, t)   #cols = c("genotype", trait(measured), "pred")                                      
    }
    
    # grouped by genotype
    Res.geno <- Res %>% group_by(genotype) %>% 
      summarise(across(.col = 1:2, mean)) %>%
      as.data.frame()  
    
    ##----------- 8. Results
    
    # Individual: correlation, RMSEP, RPD   
    X = Res[,trait]
    Y = Res[,"predicted"]
    
    r = cor(X, Y) %>% round(2)
    rmse = sqrt(mean((X - Y)^2)) %>% round(2)
    rpd = RPD(X, rmse) %>% round(2)
    individual <- c(r, rmse, rpd)
    
    # Genotype : correlation, RMSEP, RPD
    X = Res.geno[,trait]
    Y = Res.geno[,"predicted"]
    
    r.g = cor(X, Y) %>% round(2)
    rmse.g = sqrt(mean((X - Y)^2)) %>% round(2)
    rpd.g  = RPD(X, rmse.g) %>% round(2)
    genotype <- c(r.g, rmse.g, rpd.g)
    
    # Grouping
    t <- data.frame(iter = c(iter,iter),
                    type = c("individual", "genotype"), 
                    r = c(r, r.g), 
                    RMSEP = c(rmse, rmse.g), 
                    RPD = c(rpd, rpd.g))
    model_stats <- rbind(model_stats, t)
  }
  print(model_stats)
  model_stats_mean <- model_stats %>%
    group_by(type) %>%
    summarise(across(
      .cols = c("r","RMSEP", "RPD"),
      .fns = mean))
  
  print(model_stats_mean)
  output["model_stats"] <- list(model_stats)
  output["model_stats_mean"] <- list(model_stats_mean)
  return(output)
}


###------------------
####--- RUN
###------------------



#----------- 1. Running
result <- list()
result[[trait]] <- PLS_stats(trait = trait, ncomp = ncomp)


#----------- 2. Observation
# Obs.: After carefully analyzing the data, we found that every trait has one
# number of components that better fit the data.
# For the list of traits, a different component will be selected.
# Specify the trait and ncomp below it

# # (KW)
trait = "weight_mg"
ncomp = 11

# # (ST - %)
trait = "starch_PCT"
ncomp = 9

# # (ST - mg)
trait = "starch_mg"
ncomp = 9

# # (SUC)
trait = "sucrose_mg"
ncomp = 5

# # (GC)
trait = "glucose_PCT"
ncomp = 8

# # (PG)
trait = "PG_PCT"
ncomp = 8

# # (TS)
trait = "total_sugars_mg"
ncomp = 9

# # (TC)
trait = "total_carbohydrates_mg"
ncomp = 10



  
  
