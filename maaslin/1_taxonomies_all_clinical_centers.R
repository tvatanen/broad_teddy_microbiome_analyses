require(microbiomics)
require(Maaslin)
require(tidyverse)

species <- read_metaphlan_table("teddy_metaphlan_profiles.tsv", lvl = 7)
species <- species[ , apply(species,2,function(x) sum(x>0)) > 30 , ]
colnames(species) <- sapply(colnames(species), function(x) strsplit(x,".",fixed = T)[[1]][7])

load("teddy_metadata.RData")

metadata$sample_id <- paste("X", metadata$sample_id, sep = "")
# remove samples where breastfeeding data is not available
metadata <- metadata[ !(is.na(metadata$breastfeeding)) , ]
metadata <- metadata[ metadata$sample_id %in% rownames(species) , ]
species <- species[ as.character(metadata$sample_id) , ]
all(rownames(species) == as.character(metadata$sample_id))

metadata_ia <- metadata %>% 
  filter(!(is.na(IA_casecontrol_ind))) %>%
  group_by(subject_id) %>% 
  filter(age_at_collection == min(age_at_collection)) %>%
  ungroup() %>% 
  group_by(IA_casecontrol_ind) %>% 
  summarize(num_subjects=n()) %>%
  filter(num_subjects == 2) %>% 
  left_join(metadata)

## generate age wrt. first AAB
for (i in metadata_ia$IA_casecontrol_ind) {
  case_ages <- metadata_ia$age_first_pos[metadata_ia$IA_casecontrol_ind == i & metadata_ia$IA_casecontrol_outcome == 1]
  stopifnot(length(unique(case_ages)) == 1)
  seroconvertion_age <- unique(case_ages)
  metadata_ia$age_wrt_first_pos[ metadata_ia$IA_casecontrol_ind == i ] <- as.integer(metadata_ia$age_at_collection[ metadata_ia$IA_casecontrol_ind == i ] - seroconvertion_age)
}
metadata_ia %>% filter(age_wrt_first_pos >=0, age_wrt_first_pos < 30) %>% dplyr::select(age_wrt_first_pos)
# many samples where age_wrt_firstpos between [0,30]

metadata_ia <- metadata_ia %>% filter(age_wrt_first_pos <= 30)
species_ia <- species[ metadata_ia$sample_id , ]

# IA casecontrol outcome
variables <- c("age_at_collection","delivery_simple","cc","breastfeeding","solid_food","IA_casecontrol_outcome","num_wgs_reads")
#select_IA_cohort <- which(!(is.na(metadata$IA_casecontrol_ind)))

res <- Maaslin.wrapper(species_ia, 
                       metadata_ia, 
                       variables = variables, 
                       dMinSamp = 0.1, 
                       fOmitLogFile = T, 
                       strForcedPredictors = c("IA_casecontrol_outcome"), 
                       strVerbosity = "ERROR", 
                       strRandomCovariates = c("subject_id"))

# write results:
write.table(res$IA_casecontrol_outcome, file="1_results_IA.txt", quote = F, row.names = F, sep = "\t")

# one case does not have breastfeeding info, remove also this control
metadata_t1d <- metadata %>% 
  filter(!(is.na(T1D_casecontrol_ind))) %>%
  group_by(subject_id) %>% 
  filter(age_at_collection == min(age_at_collection)) %>%
  ungroup() %>% 
  group_by(T1D_casecontrol_ind) %>% 
  summarize(num_subjects=n()) %>%
  filter(num_subjects == 2) %>% 
  left_join(metadata)
  
# generate age wrt. T1D diagnosis
for (i in metadata_t1d$T1D_casecontrol_ind) {
  case_ages <- metadata_t1d$age_t1d[metadata_t1d$T1D_casecontrol_ind == i & metadata_t1d$T1D_casecontrol_outcome == 1]
  stopifnot(length(unique(case_ages)) == 1)
  seroconvertion_age <- unique(case_ages)
  metadata_t1d$age_wrt_t1d[ metadata_t1d$T1D_casecontrol_ind == i ] <- as.integer(metadata_t1d$age_at_collection[ metadata_t1d$T1D_casecontrol_ind == i ] - seroconvertion_age)
}
metadata_t1d %>% filter(age_wrt_t1d >=0, age_wrt_t1d < 30) %>% dplyr::select(age_wrt_t1d)
# many samples where age_wrt_firstpos between [0,30]

metadata_t1d <- metadata_t1d %>% filter(age_wrt_t1d <= 30)
species_t1d <- species[ metadata_t1d$sample_id , ]

# T1D case-control outcome
variables <- c("age_at_collection","delivery_simple","cc","breastfeeding","solid_food","T1D_casecontrol_outcome","num_wgs_reads")
#select_T1D_cohort <- which(!(is.na(metadata$T1D_casecontrol_ind)))

res <- Maaslin.wrapper(species_t1d, 
                       metadata_t1d, 
                       variables = variables, 
                       dMinSamp = 0.1, 
                       fOmitLogFile = T, 
                       strForcedPredictors = c("T1D_casecontrol_outcome"), 
                       strVerbosity = "ERROR", 
                       strRandomCovariates = c("subject_id"))

# write results:
write.table(res$T1D_casecontrol_outcome, file="1_results_T1D.txt", quote = F, row.names = F, sep = "\t")