require(microbiomics)
require(Maaslin)
require(tidyverse)
require(gridExtra)

data <- read.table("teddy_metacyc_pathways.tsv", header=T, sep="\t")
data <- data.frame(t(data))

load("teddy_metadata.RData")

metadata <- metadata[ !(is.na(metadata$num_wgs_reads)) , ]
metadata <- metadata[ metadata$num_wgs_reads > 5e6 , ]
metadata <- metadata[ !(is.na(metadata$breastfeeding)) , ]

metadata$sample_id <- paste("X", metadata$sample_id, sep = "")
metadata <- metadata[ metadata$sample_id %in% rownames(data) , ]
data <- data[ metadata$sample_id , ]
all(rownames(data) == metadata$sample_id)
data <- data / ( metadata$num_wgs_reads / 1e6 )

pseudocount <- 2^6

# Investigate the effect of pseudocount in median vs. variance
df <- data.frame(var=apply(log2(data+1),2,sd), avg=apply(log2(data+1),2,median))
p1 <- ggplot(df, aes(x=avg, y=var)) + geom_point() + theme_bw() + xlab("mean(log2(data+1))") + ylab("standard deviation(log2(data+1))") + coord_cartesian(xlim=c(0,10))
df <- data.frame(var=apply(log2(data+pseudocount),2,sd), avg=apply(log2(data+pseudocount),2,median))
p2 <- ggplot(df, aes(x=avg, y=var)) + geom_point() + theme_bw() + xlab("mean(log2(data+2^6))") + ylab("standard deviation(log2(data+2^6))") + coord_cartesian(xlim=c(6,8.5))
grid.arrange(p1,p2, ncol=2)

# IA casecontrol outcome
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
data_ia <- data[ metadata_ia$sample_id , ]

# filter out sparse pathways
data_ia <- data_ia[ , apply(data_ia,2,function(x) sum(x>0)) > (nrow(data_ia) * 0.3)  ]

# Covariates to control / test
variables <- c("age_at_collection",
               "delivery_simple",
               "cc",
               "breastfeeding",
               "solid_food",
               "IA_casecontrol_outcome",
               "num_wgs_reads")

res <- Maaslin.wrapper(log2(data_ia + pseudocount), 
                       metadata_ia, 
                       variables = variables, 
                       dMinSamp = 0.1, 
                       fOmitLogFile = T, 
                       strForcedPredictors = c("IA_casecontrol_outcome"), 
                       strVerbosity = "ERROR", 
                       strRandomCovariates = c("subject_id"), 
                       strTransform = "none")

# write/save results:
write.table(res$IA_casecontrol_outcome, file="2_results_IA.txt", quote = F, row.names = F, sep = "\t")
save(res , file="2_results_IA.RData")


# T1D casecontrol outcome
metadata_t1d <- metadata %>% 
  filter(!(is.na(T1D_casecontrol_ind))) %>%
  group_by(subject_id) %>% 
  filter(age_at_collection == min(age_at_collection)) %>%
  ungroup() %>% 
  group_by(T1D_casecontrol_ind) %>% 
  summarize(num_subjects=n()) %>%
  filter(num_subjects == 2) %>% 
  left_join(metadata)

for (i in metadata_t1d$T1D_casecontrol_ind) {
  case_ages <- metadata_t1d$age_t1d[metadata_t1d$T1D_casecontrol_ind == i & metadata_t1d$T1D_casecontrol_outcome == 1]
  stopifnot(length(unique(case_ages)) == 1)
  seroconvertion_age <- unique(case_ages)
  metadata_t1d$age_wrt_t1d[ metadata_t1d$T1D_casecontrol_ind == i ] <- as.integer(metadata_t1d$age_at_collection[ metadata_t1d$T1D_casecontrol_ind == i ] - seroconvertion_age)
}
metadata_t1d %>% filter(age_wrt_t1d >=0, age_wrt_t1d < 30) %>% dplyr::select(age_wrt_t1d)
# many samples where age_wrt_firstpos between [0,30]

metadata_t1d <- metadata_t1d %>% filter(age_wrt_t1d <= 30)
data_t1d <- data[ metadata_t1d$sample_id , ]

# filter out sparse pathways
data_t1d <- data_t1d[ , apply(data_t1d,2,function(x) sum(x>0)) > (nrow(data_t1d) * 0.3)  ]

# Covariates to control / test
variables <- c("age_at_collection",
               "delivery_simple",
               "cc",
               "breastfeeding",
               "solid_food",
               "T1D_casecontrol_outcome",
               "num_wgs_reads")

res <- Maaslin.wrapper(log2(data_t1d + pseudocount), 
                       metadata_t1d, 
                       variables = variables, 
                       dMinSamp = 0.1, 
                       fOmitLogFile = T, 
                       strForcedPredictors = c("T1D_casecontrol_outcome"), 
                       strVerbosity = "ERROR", 
                       strRandomCovariates = c("subject_id"), 
                       strTransform = "none")

# write/save results:
write.table(res$T1D_casecontrol_outcome, file="2_results_T1D.txt", quote = F, row.names = F, sep = "\t")
save(res , file="2_results_T1D.RData")
