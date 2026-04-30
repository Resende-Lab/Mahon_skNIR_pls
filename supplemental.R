
######################################################################
#  get result from RUNME_NIR0428.R
######################################################################
# for Table 1 
#data: result, list = (sample_stats, model_stats, model_stats_mean)

library(dplyr)

# sample_stats retrieved from result list
traits <- names(result)
stats <- c()
for (i in 1:length(traits)) {
  t <- result[[i]][[1]] |>  as.data.frame()
  t$trait <- traits[i]
  stats <- rbind(stats, t)
}
# remove "weight_mg genotype" row
stats <- stats[-2,]
sample_stats <- stats[c(7, 1:6)]

# model_stats_mean retrieved from result list
stats <- c()
for (i in 1:length(traits)) {
  t <- result[[i]][[3]] |>  as.data.frame()
  t$trait <- traits[i]
  stats <- rbind(stats, t)
}

# remove weight_mg, genotype row
stats <- stats[-1,]
model_stats_mean <- stats[c(5, 1:4)]
model_stats_mean  <-  model_stats_mean |> 
  mutate(across(where(is.numeric), \(x) round(x, 2)))

#save(result, sample_stats, model_stats_mean, file = "my_data.RData")
#names <- c("KW(mg)","ST(%)", "ST(mg)","SUC(mg)","GC(PCT)","PG(PCT)","TS(mg)","TC(mg)")



########################
# Supplemental figure
########################

#============================================================================
#  SD plot 
#     -compare sd of PG (Wet chemistry) and sd of NIR_PG (prediction)
#============================================================================
# calculate Coefficient of Variation as a decimal and sd
# plot SD comparison


library(ggplot2)

SNV <- function(x) {
  nir <- as.matrix(x)
  nir_snv <- t(scale(t(nir), center = TRUE, scale = TRUE))
  return(nir_snv)
}
cv <- function(x) {
  sd(x, na.rm = TRUE)/mean(x, na.rm = TRUE)
}

k_fold = 10

# data prep
dt <- read.csv("PG_WL.NIR.csv") |> 
  rename(PG_i = PG_PCT) |> 
  filter(!is.na(PG_i))
nir <- dt |>  
  select(X940:X1640) 
WLi_NIRi <- data.frame(dt[,c("genotype","PG_i")], NIR = I(SNV(nir))) 
set.seed(711)
dat <- WLi_NIRi %>%  mutate(set = sample(1:k_fold, nrow(WLi_NIRi), replace = TRUE))

# Create a dataframe to record predicted values
Res = data.frame(genotype = as.character(),
                 PG_i = as.numeric(),
                 pred_i = as.numeric())

#------- Model

for (fold in 1:k_fold) {
  #split dataset into test and train
  test = dat %>% filter(set == fold)
  train = dat %>% filter(set != fold)
  
  # build train.pls model for each fold
  train.pls <- plsr(PG_i ~ NIR, ncomp = 25, data = train, validation = "CV")
  
  # predict test data (outer-validation) with train.pls
  t <- test %>% 
    select(genotype, PG_i) %>% 
    mutate(pred_i = predict(train.pls, ncomp = 8, newdata = test) )
  # data
  Res <- rbind(Res, t) 
}
Res <- Res[order(Res$genotype),]

# calculate mean, variance, CV by genotype for PG and NIR_PG
t1 <- Res[c("genotype", "PG_i", "pred_i")] |> 
  na.omit() |> 
  group_by(genotype) |> 
  summarise(
    across(.cols = c("PG_i", "pred_i"),
           .fns = list(mean = mean,
                       sd = sd, 
                       cv = cv)))

# average 
t1_avg <- t1[-1] |> 
  summarise_all(mean)
t1_avg

# SD plot
t1 |> 
  ggplot(aes(x = pred_i_sd, y = PG_i_sd)) +
  geom_point() +
  geom_abline(slope = 1) +
  xlab("SD of NIR_PG") +
  ylab("SD of PG") +
  xlim(0,9) +
  ylim(0,10) +
  ggtitle("SD comparison")

#==========================================================================
#  1k training model using avg PG 
#==========================================================================

# New Figure: build a new model using average PG as response values and
#                                     1k NIR(averaged spec) as predictor

# data prep (WLg_NIRi)
df <- read.csv("PG_WL.NIR.csv") |> 
  rename(PG_i = PG_PCT) |> 
  filter(!is.na(PG_i))
NIR <- select(df, X940:X1640) |> SNV()
WLi <- df[,c("genotype","PG_i")] 
WLg <- WLi |> 
  group_by(genotype) |> 
  summarise(PG_g = mean(PG_i, na.rm = TRUE))
WLg_NIRi <- data.frame(left_join(WLi,WLg, by = "genotype"), I(NIR)) 
dt <- WLg_NIRi 

genos <- unique(dt$genotype)
cor <- data.frame(iter = as.integer(),
                  r = as.numeric())
RES <- list()
set.seed(3092)
iters = 5

# create a new dataset by sampling 1 kernel for each genotype
data <- c()
for (i in genos) {
  d <- dt[dt$genotype == i, ]
  t <- d[sample(nrow(d), 1),]
  data <- rbind(data, t)
}

# build PLS model
cat("iter: ")
for (j in 1:iters) {
  set.seed(j * 5 + 1000)
  cat(j,",")
  # conditions
  k_fold = 10
  ncomp = 4 
  set.seed(j*4)
  # Create a dataframe to record predicted values
  Res = data.frame(genotype = as.character(),
                   PG_g = as.numeric(),
                   pred_i = as.numeric())
  dat <- data %>%  mutate(set = sample(1:k_fold, nrow(data), replace = TRUE))
  
  #------- Model
  for (fold in 1:k_fold) {
    #split dataset into test and train
    test = dat %>% filter(set == fold)
    train = dat %>% filter(set != fold)
    
    # build train.pls model for each fold
    train.pls <- plsr(PG_g ~ NIR, ncomp = 25, data = train, validation = "CV")
    # predict test data (outer-validation) with train.pls
    t <- test %>% 
      select(genotype, PG_g) %>% 
      mutate(pred_i = predict(train.pls, ncomp = ncomp, newdata = test) )
    Res <- rbind(Res, t) 
  }
  Res <- Res[order(Res$genotype),]
  RES[[paste0("iter", j)]] <- Res
  cor[j,"iter"] <- j
  cor[j,"r"] <- cor(Res$PG_g, Res$pred_i)
} 
cor
cor_mean <- mean(cor$r) |> round(2)  


# scatter plot for 1k-training model Res(iter 5)
label = paste0("r = ", cor_mean)  
Res |> 
  ggplot(aes(x = pred_i, y = PG_g)) +
  geom_point() +
  xlab("NIR_PG (indi)") +
  ylab("PG (geno)") +
  # xlim(0,30) +
  # ylim(0,30) +
  geom_smooth(method = "lm", se = FALSE) +
  annotate("text", x = 5, y = 20, label = label, parse = FALSE) +
  ggtitle("1k-training model with PG(geno)" )



