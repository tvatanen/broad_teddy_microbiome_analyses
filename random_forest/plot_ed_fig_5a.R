library(ggplot2)
library(RColorBrewer)

# load pathway classification results
load("forests_T1D_pathways.RData")
df_errors <- data.frame(err=errors, test="pathways", stringsAsFactors = F)
median(errors)
# median error rate 0.449 for taxonomies

# load taxonomioc classification results
load("forests_T1D_taxa.RData")
df_errors <- rbind(df_errors,
                data.frame(err=errors, test="taxonomies", stringsAsFactors = F))
median(errors)
# median error rate 0.450 for pathways

# Extended Data Figure S5A
ggplot(df_errors, aes(y=err,x=test)) + 
  geom_boxplot(notch = T, outlier.colour = NA) + 
  ylab("Error rate per case-control pair") + theme_bw() + 
  geom_dotplot(binwidth = 0.025, stackdir = "center", binaxis = "y", alpha=0.5)
