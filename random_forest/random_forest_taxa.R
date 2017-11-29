load("teddy_metadata.RData")
library(tidyverse)
library(randomForest)
library(microbiomics)

species <- read_metaphlan_table("teddy_metaphlan_profiles.tsv", lvl = 7)
species <- species[ , apply(species,2,function(x) sum(x>0)) > 30 , ]

rownames(species) <- substr(rownames(species),2,100)
colnames(species) <- sapply(colnames(species), function(x) strsplit(x,".",fixed = T)[[1]][7])

metadata$sample_id <- as.character(metadata$sample_id)
metadata <- metadata[ metadata$sample_id %in% rownames(species) , ]
species <- species[ as.character(metadata$sample_id) , ]
all(rownames(species) == as.character(metadata$sample_id))

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
species_t1d <- species[ metadata_t1d$sample_id , ]

# leave one case-control pair out cross-validation
set.seed(1523)
mtree <- 2000
errors <- c()
forests <- list()
i <- 0
for (leave_out in unique(metadata_t1d$T1D_casecontrol_ind)) {
  training <- metadata_t1d %>% filter(!(T1D_casecontrol_ind %in% leave_out))
  test <- metadata_t1d %>% filter(T1D_casecontrol_ind %in% leave_out)
  i <- i + 1
  
  training_data <- bind_cols(species_t1d[ training$sample_id , ], training[ , c("cc","age_at_collection")]) 
  test_data <- bind_cols(species_t1d[ test$sample_id , ], test[ , c("cc","age_at_collection")])
  
  nstrat <- min(table(training$T1D_casecontrol_outcome))
  rf <- randomForest(training_data, factor(as.integer(training$T1D_casecontrol_outcome)), 
                     xtest = test_data, ytest = factor(as.integer(test$T1D_casecontrol_outcome)),
                     importance = F, ntree=mtree, strata = training$T1D_casecontrol_outcome,
                     sampsize = rep(nstrat, 2))
  forests[[i]] <- rf 
  errors <- c(errors, rf$test$err.rate[mtree,1])

}

save(forests, errors, file="forests_T1D_taxa.RData")