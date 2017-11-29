load("teddy_metadata.RData")
library(tidyverse)
library(stringr)
library(randomForest)

data <- read.table("teddy_metacyc_pathways.tsv", header=T, sep="\t")
data <- data.frame(t(data))

rownames(data) <- substr(rownames(data), 2, 100)
metadata <- metadata[ metadata$sample_id %in% rownames(data) , ]
metadata$sample_id <- as.character(metadata$sample_id)
data <- data[ metadata$sample_id , ]
all(rownames(data) == metadata$sample_id)

data <- data / (metadata$num_wgs_reads / 1e6)

metadata_t1d <- metadata %>% 
  filter(!(is.na(T1D_casecontrol_ind)))

# compute age wrt. T1D diagnosis
for (i in metadata_t1d$T1D_casecontrol_ind) {
  case_ages <- metadata_t1d$age_t1d[metadata_t1d$T1D_casecontrol_ind == i & metadata_t1d$T1D_casecontrol_outcome == 1]
  stopifnot(length(unique(case_ages)) == 1)
  seroconvertion_age <- unique(case_ages)
  metadata_t1d$age_wrt_t1d[ metadata_t1d$T1D_casecontrol_ind == i ] <- as.integer(metadata_t1d$age_at_collection[ metadata_t1d$T1D_casecontrol_ind == i ] - seroconvertion_age)
}

metadata_t1d <- metadata_t1d %>% filter(age_wrt_t1d <= 30)

table(metadata_t1d$T1D_casecontrol_outcome) 
data_t1d <- data[ metadata_t1d$sample_id , ]

# leave one case-control pair out cross-validation
set.seed(1523)
mtree <- 2000
errors <- c()
forests <- list()
i <- 0
for (leave_out in unique(metadata_t1d$T1D_casecontrol_ind)) {
  training <- metadata_t1d %>% filter(!(T1D_casecontrol_ind %in% leave_out))
  test <- metadata_t1d %>% filter(T1D_casecontrol_ind %in% leave_out)
  if (length(unique(test$T1D_casecontrol_outcome)) < 2) { next }
  i <- i + 1
  
  training_data <- bind_cols(data_t1d[ training$sample_id , ], training[ , c("cc","age_at_collection")]) 
  test_data <- bind_cols(data_t1d[ test$sample_id , ], test[ , c("cc","age_at_collection")])

  nstrat <- min(table(training$T1D_casecontrol_outcome))
  rf <- randomForest(training_data, factor(as.integer(training$T1D_casecontrol_outcome)), 
                     xtest = test_data, ytest = factor(as.integer(test$T1D_casecontrol_outcome)),
                     importance = T, ntree=mtree, strata = training$T1D_casecontrol_outcome,
                     sampsize = rep(nstrat, 2))
  forests[[i]] <- rf 
  errors <- c(errors, rf$test$err.rate[mtree,1])

}

save(forests, errors, file="forests_T1D_pathways.RData")