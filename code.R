setwd("add directory with data")

library(tidyverse)
library(mixOmics)
library(ggpubr)
library(vegan)
library(caret)
library(rstatix)
library(UpSetR)
library(homologueDiscoverer)
library(ggdensity)
library(cowplot)
library(readxl)
library(lmerTest)

# Function for homologues analysis
plotAnnotatedStatic <- function(annotated, legend_setting = "bottom"){
  
  annotated <- mutate(annotated,
                      homologue_id = if_else(is.na(homologue_id), 0L, homologue_id),
                      homologue_id = as.factor(homologue_id))
  colourCount = length(unique(annotated$homologue_id))
  getPalette = colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
  ncolor = length(unique(annotated$homologue_id))
  annotated <- arrange(annotated, desc(homologue_id))
  g <- ggplot(annotated, aes(group = homologue_id)) +
    geom_point(aes(x = rt, y = mz, color = homologue_id, size = homologue_id,
                   shape = homologue_id, alpha = homologue_id)) +
    geom_line(data = filter(annotated, homologue_id != 0),
              aes(x = rt, y = mz, group = homologue_id), alpha = 0.4, linewidth = 0.1) +
    ggtitle("Annotated Peak Table") +
    scale_colour_manual(values = c("black", getPalette(colourCount-1))) +
    scale_alpha_manual(values =  c(0.2, rep(0.85, times = ncolor))) +
    scale_size_manual(values = c(0.1, rep(1.5, times = ncolor))) +
    scale_shape_manual(values = c(1, rep(19, times = ncolor))) +
    xlab("Retention Time (s)") +
    ylab("Mass to Charge Ratio") +
    theme(legend.position=legend_setting, text = element_text(family="mono"))
  return(g)
  
}


# Read data
feature_table <- read_csv("gnps_quant.csv")
metadata <- read_csv("metadata.csv")
sample_order <- read.csv("sequence.csv") %>% 
  dplyr::select(1,3,9)
annotations <- read.delim("fbmn_gnps2.tsv") %>%
  dplyr::filter(!str_detect(pattern = "REFRAME", LibraryName)) # remove drug library
annotations$X.Scan. <- as.character(annotations$X.Scan.)
canopus <- read_tsv("canopus.tsv")
canopus$featureId <- as.character(canopus$featureId)

info_feature <- feature_table %>% dplyr::select(1:3,7)
colnames(info_feature) <- c("Feature", "mz", "RT", "Corr_ID")
info_feature$Feature <- as.character(info_feature$Feature)

info_feature_complete <- info_feature %>% 
  left_join(annotations, by = c("Feature" = "X.Scan.")) %>%
  dplyr::select(1:4,18,24) %>% 
  left_join(canopus %>% distinct(featureId, .keep_all = TRUE), by = c("Feature" = "featureId"))

metadata_final <- metadata %>% left_join(sample_order, by = c("SampleID" = "File.Name"))


# Extract feature table 
data <- feature_table %>% 
  column_to_rownames("row ID") %>% dplyr::select(contains("Peak")) %>% 
  t() %>% as.data.frame() %>% rownames_to_column("SampleID") %>% 
  arrange(SampleID) %>% distinct(SampleID, .keep_all = TRUE)

data$SampleID <- gsub(".mzML Peak area", "", data$SampleID)

# Discard first blanks, pools and mix used for conditioning
data_final <- data %>% 
  dplyr::filter(!(SampleID %in% c("Blank_1", "Blank_2", "SixMix_0", 
                                  "PoolQC_1", "PoolQC_2", "PoolQC_3")))

# Calculate TIC
data_TIC <- data.frame(TIC = rowSums(data_final %>% column_to_rownames("SampleID"))) %>% 
  rownames_to_column("SampleID") %>% left_join(sample_order, by = c("SampleID" = "File.Name"))

data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|Six|Pool", SampleID))) %>%
  left_join(metadata_final) %>%
  ggscatter("Order", "TIC", color = "Donor", add = "reg.line") + ylim(0, 1e9) +
  stat_cor()
# TIC in DonorA is lower. No strong effect of run order


# Check sample type
sample_tic <- data_TIC %>% dplyr::filter(!(str_detect(pattern = "Blank|Pool|Mix", SampleID))) %>% summarise(mean(TIC))
pool_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", SampleID)) %>% summarise(mean(TIC))
six_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "SixMix", SampleID)) %>% summarise(mean(TIC))
blank_tic <- data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% summarise(mean(TIC))

# Ratios
ratio_tic_pb <- pool_tic/sample_tic
ratio_tic_sb <- sample_tic/blank_tic
ratio_tic_6p <- pool_tic/six_tic

# Check TIC overtime for QCpool, QCmix, and Blank
data_TIC %>% dplyr::filter(str_detect(pattern = "Pool", SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1.5e9) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "SixMix", SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 1.5e8) +
  stat_cor()

data_TIC %>% dplyr::filter(str_detect(pattern = "Blank", SampleID)) %>% 
  ggscatter("Order", "TIC", add = "reg.line") + ylim(0, 3e7) +
  stat_cor()


# Check possible acquired internal standards (IS)
fbmn_IS <- annotations %>% 
  dplyr::filter(str_detect(Compound_Name, regex("sulf", ignore_case = TRUE))) %>% 
  distinct(X.Scan., .keep_all = TRUE) %>% dplyr::filter(Organism != "BILELIB19")

# Extract Sulfamethazine (IS) from the table
table_IS <- data_final %>% column_to_rownames("SampleID") %>% t() %>% 
  as.data.frame() %>% rownames_to_column("ID") %>% dplyr::filter(ID %in% fbmn_IS$X.Scan.) %>% 
  column_to_rownames("ID") %>% t() %>% as.data.frame() %>% 
  rownames_to_column("filename") %>% left_join(metadata_final, by = c("filename" = "SampleID")) %>%
  dplyr::filter(!(str_detect(filename, "SixMix|Blank|Pool"))) %>%
  dplyr::select(filename, `974`, Order, Info)

colnames(table_IS)[2] <- "Sulfamethazine"

table_IS %>% ggscatter(x = "Order", y = "Sulfamethazine", add = c("reg.line")) +
  ylim(0, 2e6) +
  stat_cor()

table_IS %>% ggbarplot(x = "Order", y = "Sulfamethazine", xlab = "Run Order", 
                       ylab = "Peak Area Sulfamethazine", title = "Internal Standard Acquisition") +
  geom_hline(yintercept = mean(table_IS$Sulfamethazine, na.rm = TRUE), linetype = "dashed", color = "blue")

cv_is <- sd(table_IS$Sulfamethazine)/mean(table_IS$Sulfamethazine)
# CV is 46% because volumes of solvent added to each 
# sample was adjusted based on sample weight


# Check extracted features in Blanks, QCpool, and QCmix
data_blank <- data_final %>% dplyr::filter(str_detect(pattern = "Blank", SampleID))
data_pools <- data_final %>% dplyr::filter(str_detect(pattern = "Pool", SampleID))
data_sixmix <- data_final %>% dplyr::filter(str_detect(pattern = "Six", SampleID))

# Blanks
blanks_feature_info <- data.frame(Feature = colnames(data_blank)[-1],
                                  Mean_blank = data_blank %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_blank =  data_blank %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_blank = SD_blank/Mean_blank) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_blank, SD_blank, CV_blank) %>% 
  dplyr::filter(Mean_blank > 0) %>% arrange(desc(Mean_blank))

# Six mix
sixmix_feature_info <- data.frame(Feature = colnames(data_sixmix)[-1],
                                  Mean_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% colMeans(), 
                                  SD_sixmix = data_sixmix %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_sixmix = SD_sixmix/Mean_sixmix) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_sixmix, SD_sixmix, CV_sixmix) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% arrange(desc(Mean_sixmix))
# CVs of the six standards in QCmix are ranging between 6% to 9% --> good

# Pool
pools_feature_info <- data.frame(Feature = colnames(data_pools)[-1],
                                 Mean_pool = data_pools %>% column_to_rownames("SampleID") %>% colMeans(), 
                                 SD_pool =  data_pools %>% column_to_rownames("SampleID") %>% apply(2, sd)) %>%
  dplyr::mutate(CV_pool = SD_pool/Mean_pool) %>% left_join(annotations, by = c("Feature" = "X.Scan.")) %>% 
  left_join(info_feature) %>% dplyr::select(Feature, mz, RT, Compound_Name, SpectrumID, 
                                            Precursor_MZ, Mean_pool, SD_pool, CV_pool) %>% 
  dplyr::filter(Mean_pool > 0) %>% arrange(desc(Mean_pool))


# Blank subtraction - Remove features for which Pools/Blank < 10
feature_to_remove <- blanks_feature_info %>% left_join(pools_feature_info) %>%
  dplyr::filter(Mean_blank > 0) %>% 
  dplyr::mutate(Pool_Blank = Mean_pool/Mean_blank) %>% 
  dplyr::filter(Pool_Blank < 10 | is.na(Pool_Blank))

# Clean feature table
data_clean <- data_final %>% dplyr::select(-c(feature_to_remove$Feature)) %>%
  column_to_rownames("SampleID") %>%
  select_if(~sum(.) != 0) %>% rownames_to_column("SampleID")

# QCMix subtraction - Remove features for which Pool/QCmix < 10
feature_to_remove_mix <- sixmix_feature_info %>% left_join(pools_feature_info) %>% 
  dplyr::filter(Mean_sixmix > 0) %>% 
  dplyr::mutate(Pool_Mix = Mean_pool/Mean_sixmix) %>% 
  dplyr::filter(Pool_Mix < 10 | is.na(Pool_Mix)) %>% dplyr::filter(!(Feature %in% feature_to_remove$Feature))

# Clean feature table
data_clean2 <- data_clean %>% dplyr::select(-c(feature_to_remove_mix$Feature))


# PCA on raw data
PCA_raw <- mixOmics::pca(data_clean2 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Time"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal Metabolome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()


# Plastic contamination was observed in one donor - make sure to clean data from polymers
peg <- data_final %>% dplyr::select(SampleID, `1073`) %>% 
  left_join(metadata_final) 
colnames(peg)[2] <- "PEG"

peg %>% arrange(Order) %>% dplyr::mutate(LogPeg = log2(PEG+1)) %>%
  ggboxplot(x = "Donor", y = "PEG", add = "jitter",
            color = "Order", group = "Donor", legend = "right") + 
  scale_color_viridis_c()

# Check features
features_to_check <- info_feature %>% dplyr::select(1:3) %>% 
  dplyr::filter(Feature %in% colnames(data_clean)) %>% as_tibble()
colnames(features_to_check) <- c("peak_id", "mz", "rt")
features_to_check$mz <- as.double(features_to_check$mz)
features_to_check$rt <- as.double(features_to_check$rt)

feature_pool <- data_clean2 %>% 
  dplyr::filter(SampleID == "PoolQC_5") %>% # representative sample
  column_to_rownames("SampleID") %>%
  t() %>% as.data.frame() %>% rownames_to_column("peak_id")

features_to_check_pool <- features_to_check %>% left_join(feature_pool)

features_to_check_pool$peak_id <- as.integer(features_to_check_pool$peak_id)
colnames(features_to_check_pool)[4] <- "intensity"

# Detect PEGs
peak_table_PEG <- detectHomologues(features_to_check_pool, mz_steps = c(44.02628),
                                   rt_min = 0.1, rt_max = 100, ppm_tolerance = 10, 
                                   min_series_length = 5, search_mode = "targeted", 
                                   step_mode = "increment", verbose = TRUE)

plotAnnotatedStatic(peak_table_PEG)

sdb <- sdbCreate(peak_table_PEG, sample_origin = "data")

# Detect polymers in an untargeted fashion
peak_table_untargeted <- detectHomologues(features_to_check_pool, mz_min = 10, mz_max = 50, 
                                          rt_min = 0.5, rt_max = 100, ppm_tolerance = 5, 
                                          min_series_length = 5, search_mode = "untargeted", 
                                          step_mode = "increment", verbose = TRUE)

plotAnnotatedStatic(peak_table_untargeted)

# Combine obtained tables
overlapping_ids <- sdbCheckContained(sdb, peak_table_untargeted, ppm_tolerance = 0.5, rt_tolerance = 0.5)
sdb <- sdbPush(sdb, filter(peak_table_untargeted, !(homologue_id %in% overlapping_ids)))

sdb$sample_peak_id <- as.character(sdb$sample_peak_id)
features_contaminat <- info_feature_complete %>% dplyr::filter(Feature %in% sdb$sample_peak_id)


# Final cleaned table
data_clean3 <- data_clean2 %>% dplyr::select(-c(sdb$sample_peak_id))


# PCA on cleaned table
PCA_raw <- mixOmics::pca(data_clean3 %>% select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>%
                           column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_raw_scores <- data.frame(PCA_raw$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_raw_plot <- PCA_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Fecal Metabolome", i, sep = " "),
            xlab = paste("PC1 (", round(PCA_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# Keep only samples
data_sample <- data_clean3 %>% 
  dplyr::filter(!(str_detect(pattern = "Six|Blank|Pool", SampleID))) %>% 
  column_to_rownames("SampleID") %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% rownames_to_column("SampleID")

# RCLR transformation
data_sample_clr <- decostand(data_sample %>% column_to_rownames("SampleID"), method = "rclr")

# PCA
PCA_whole <- mixOmics::pca(data_sample_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_whole_scores <- data.frame(PCA_whole$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)
PCA_whole_scores$Time <- factor(PCA_whole_scores$Time, levels = c("im_2mn", "imd_4mn", "1dy", "1wk", "DLAB"))

i <- "Time"

PCA_plot <- PCA_whole_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - Extraction Method"),
            xlab = paste("PC1 (", round(PCA_whole$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_whole$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_whole_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed() + scale_color_viridis_d()

# PERMANOVA
dist_metabolites <- vegdist(data_sample_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, PCA_whole_scores$Donor)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ Donor * Time * Weight, PCA_whole_scores, na.action = na.omit)


#####################################
# Figure 2 - Comparison MeOH v EtOH #
#####################################

# Extract 20mg triplicates for 95% EtOH and 50% MeOH 
sample_20mg <- metadata_final %>% 
  dplyr::filter(Weight == "20mg") %>% 
  dplyr::filter(Time == "im_2mn" | Time == "DLAB")

data_20mg <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_20mg$SampleID) %>% 
  column_to_rownames("SampleID") %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  rownames_to_column("SampleID")

# RCLR transformation
data_20mg_clr <- decostand(data_20mg %>% column_to_rownames("SampleID"), method = "rclr")

# PCA of 20mg samples
PCA_20mg <- mixOmics::pca(data_20mg_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_20mg_scores <- data.frame(PCA_20mg$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Time"
j <- "Donor"

PCA_20mg_plot <- PCA_20mg_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6,
            title = paste("PCA - Fecal Metabolic Profiles (20 mg)"),
            xlab = paste("PC1 (", round(PCA_20mg$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_20mg$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_20mg_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  geom_point(data = PCA_20mg_scores %>% group_by((!!sym(j))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(j))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PERMANOVA
dist_metabolites <- vegdist(data_20mg_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, PCA_20mg_scores$Donor)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ Donor * Time, PCA_20mg_scores, na.action = na.omit)

#ggsave(plot = PCA_plot, filename = "Figure2A.svg", device = "svg", dpi = "retina", width = 4)


# Generate Upset plot to check metabolic features recovered by using the different extraction methods
data_20mg_presence <- data_20mg %>% column_to_rownames("SampleID")
data_20mg_presence[data_20mg_presence > 0] <- 1
data_20mg_presence_info <- data_20mg_presence %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Time)) %>% column_to_rownames("SampleID")

list_20mg <- list(`im_2mn` = data_20mg_presence_info %>% dplyr::filter(Time == "im_2mn") %>% 
                    dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% colnames(),
                  `DLAB` = data_20mg_presence_info %>% dplyr::filter(Time == "DLAB") %>% 
                    dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% colnames())

upset_20mg <- upset(fromList(list_20mg), order.by = "freq") # save pdf for Figure 2B


# Identified features only detected via 95% EtOH or 50% MeOH
features_etoh <- data_20mg_presence_info %>% dplyr::filter(Time != "DLAB") %>%
  dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

features_meoh <- data_20mg_presence_info %>% dplyr::filter(Time == "DLAB") %>%
  dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

exclusive_etoh <- features_etoh %>% anti_join(features_meoh) %>% left_join(info_feature_complete)
exclusive_meoh <- features_meoh %>% anti_join(features_etoh) %>% left_join(info_feature_complete)

intersect_extractions <- features_etoh %>% inner_join(features_meoh) %>% 
  left_join(info_feature_complete)

all_extractions <- features_etoh %>% full_join(features_meoh) %>% 
  left_join(info_feature_complete) %>% distinct(Feature, .keep_all = TRUE) %>%
  dplyr::filter(!(is.na(Compound_Name)))

# Superclass prediction via CANOPUS
etoh_pred <- exclusive_etoh %>% group_by(`NPC#superclass`) %>% 
  summarise(Count_etoh = n()) %>% arrange(desc(Count_etoh)) %>%
  dplyr::mutate(Ratio_etoh = Count_etoh/sum(Count_etoh))

colnames(etoh_pred)[1] <- "NPC Superclass"

meoh_pred <- exclusive_meoh %>% group_by(`NPC#superclass`) %>% 
  summarise(Count_meoh = n()) %>% arrange(desc(Count_meoh)) %>%
  dplyr::mutate(Ratio_meoh = Count_meoh/sum(Count_meoh))

colnames(meoh_pred)[1] <- "NPC Superclass"

pred_join <- etoh_pred %>% 
  full_join(meoh_pred) %>%
  dplyr::mutate(Delta_count = abs(Count_meoh - Count_etoh)) %>%
  dplyr::mutate(Delta_ratio = abs(Ratio_meoh - Ratio_etoh)) %>%
  arrange(desc(Delta_ratio))

#write_csv(x = pred_join, file = "SupplementaryTable1.csv")

# Pie chart for NPC Pathway
etoh_pred <- exclusive_etoh %>% group_by(`NPC#pathway`) %>% 
  summarise(Count_etoh = n()) %>% arrange(desc(Count_etoh)) %>%
  dplyr::mutate(Ratio_etoh = Count_etoh/sum(Count_etoh))

colnames(etoh_pred)[1] <- "NPC Pathway"

pie_etoh <- etoh_pred %>% ggpie(x = "Count_etoh", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

meoh_pred <- exclusive_meoh %>% group_by(`NPC#pathway`) %>% 
  summarise(Count_meoh = n()) %>% arrange(desc(Count_meoh)) %>%
  dplyr::mutate(Ratio_meoh = Count_meoh/sum(Count_meoh))

colnames(meoh_pred)[1] <- "NPC Pathway"

pie_meoh <- meoh_pred %>% ggpie(x = "Count_meoh", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

pie_combined <- ggarrange(pie_etoh, pie_meoh, common.legend = TRUE, legend = "right")

#ggsave(plot = pie_combined, filename = "Figure2B_pies.svg", device = "svg")


# Generate pairwise PLS-DA models between donors for each extraction method

# 50% MeOH - Subject A v Subject B
meoh_AB <- data_20mg %>% dplyr::filter(str_detect(pattern = "A|B", SampleID)) %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_meoh_AB <- mixOmics::pca(meoh_AB, ncomp = 2, center = TRUE, scale = TRUE)
PCA_meoh_AB_scores <- data.frame(PCA_meoh_AB$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_meoh_AB_plot <- PCA_meoh_AB_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 50% MeOH - Subject A v Subject B"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_meoh_AB$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_meoh_AB$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_meoh_AB_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_meoh_ab <- mixOmics::plsda(meoh_AB, PCA_meoh_AB_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_meoh_ab_scores <- data.frame(PLSDA_meoh_ab$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_meoh_ab_plot <- PLSDA_meoh_ab_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 50% MeOH - Subject A v Subject B", 
            xlab = paste("Component 1 (", round(PLSDA_meoh_ab$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_meoh_ab$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_meoh_ab_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_meoh_ab <- plotLoadings(PLSDA_meoh_ab, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_meoh_ab <- perf(PLSDA_meoh_ab, validation = "loo") 
plot(perf_plsda_meoh_ab, legend = FALSE)

VIPs_meoh_ab <- as.data.frame(mixOmics::vip(PLSDA_meoh_ab))
VIPs_meoh_ab_filter <- dplyr::filter(VIPs_meoh_ab, VIPs_meoh_ab$comp1 > 1)
VIPs_meoh_ab_filter$ID <- rownames(VIPs_meoh_ab_filter)
VIPs_meoh_ab_select <- VIPs_meoh_ab_filter %>% dplyr::select(ID, comp1)
VIPs_meoh_ab_Load <- VIPs_meoh_ab_select %>% 
  left_join(Loadings_meoh_ab, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

# 50% MeOH - Subject B v Subject C
meoh_BC <- data_20mg %>% dplyr::filter(str_detect(pattern = "B|C", SampleID)) %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_meoh_BC <- mixOmics::pca(meoh_BC, ncomp = 2, center = TRUE, scale = TRUE)
PCA_meoh_BC_scores <- data.frame(PCA_meoh_BC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_meoh_BC_plot <- PCA_meoh_BC_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 50% MeOH - Subject B v Subject C"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_meoh_BC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_meoh_BC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_meoh_BC_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_meoh_BC <- mixOmics::plsda(meoh_BC, PCA_meoh_BC_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_meoh_BC_scores <- data.frame(PLSDA_meoh_BC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_meoh_bc_plot <- PLSDA_meoh_BC_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 50% MeOH - Subject B v Subject C", 
            xlab = paste("Component 1 (", round(PLSDA_meoh_BC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_meoh_BC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_meoh_BC_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_meoh_bc <- plotLoadings(PLSDA_meoh_BC, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_meoh_bc <- perf(PLSDA_meoh_BC, validation = "loo") 
plot(perf_plsda_meoh_bc, legend = FALSE)

VIPs_meoh_bc <- as.data.frame(mixOmics::vip(PLSDA_meoh_BC))
VIPs_meoh_bc_filter <- dplyr::filter(VIPs_meoh_bc, VIPs_meoh_bc$comp1 > 1)
VIPs_meoh_bc_filter$ID <- rownames(VIPs_meoh_bc_filter)
VIPs_meoh_bc_select <- VIPs_meoh_bc_filter %>% dplyr::select(ID, comp1)
VIPs_meoh_bc_Load <- VIPs_meoh_bc_select %>% 
  left_join(Loadings_meoh_bc, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

# 50% MeOH - Subject A v Subject C
meoh_AC <- data_20mg %>% dplyr::filter(str_detect(pattern = "A|C", SampleID)) %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_meoh_AC <- mixOmics::pca(meoh_AC, ncomp = 2, center = TRUE, scale = TRUE)
PCA_meoh_AC_scores <- data.frame(PCA_meoh_AC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_meoh_AC_plot <- PCA_meoh_AC_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 50% MeOH - Subject A v Subject C"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_meoh_AC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_meoh_AC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_meoh_AC_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_meoh_ac <- mixOmics::plsda(meoh_AC, PCA_meoh_AC_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_meoh_ac_scores <- data.frame(PLSDA_meoh_ac$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_meoh_ac_plot <- PLSDA_meoh_ac_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 50% MeOH - Subject A v Subject C", 
            xlab = paste("Component 1 (", round(PLSDA_meoh_ac$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_meoh_ac$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_meoh_ac_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_meoh_ac <- plotLoadings(PLSDA_meoh_ac, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_meoh_ac <- perf(PLSDA_meoh_ac, validation = "loo") 
plot(perf_plsda_meoh_ac, legend = FALSE)

VIPs_meoh_ac <- as.data.frame(mixOmics::vip(PLSDA_meoh_ac))
VIPs_meoh_ac_filter <- dplyr::filter(VIPs_meoh_ac, VIPs_meoh_ac$comp1 > 1)
VIPs_meoh_ac_filter$ID <- rownames(VIPs_meoh_ac_filter)
VIPs_meoh_ac_select <- VIPs_meoh_ac_filter %>% dplyr::select(ID, comp1)
VIPs_meoh_ac_Load <- VIPs_meoh_ac_select %>% 
  left_join(Loadings_meoh_ac, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

combined_pca_meoh <- ggarrange(PCA_meoh_AB_plot, PCA_meoh_BC_plot, PCA_meoh_AC_plot, ncol = 3)
combined_plsda_meoh <- ggarrange(PLSDA_meoh_ab_plot, PLSDA_meoh_bc_plot, PLSDA_meoh_ac_plot, ncol = 3)


# 95% Etoh - Subject A v Subject B
etoh_AB <- data_20mg %>% left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  dplyr::filter(Donor != "DonorC") %>% dplyr::select(-Donor) %>%
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_etoh_AB <- mixOmics::pca(etoh_AB, ncomp = 2, center = TRUE, scale = TRUE)
PCA_etoh_AB_scores <- data.frame(PCA_etoh_AB$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_etoh_AB_plot <- PCA_etoh_AB_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 95% Etoh - Subject A v Subject B"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_etoh_AB$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_etoh_AB$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_etoh_AB_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_etoh_ab <- mixOmics::plsda(etoh_AB, PCA_etoh_AB_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_etoh_ab_scores <- data.frame(PLSDA_etoh_ab$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_etoh_ab_plot <- PLSDA_etoh_ab_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 95% Etoh - Subject A v Subject B", 
            xlab = paste("Component 1 (", round(PLSDA_etoh_ab$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_etoh_ab$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_etoh_ab_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_etoh_ab <- plotLoadings(PLSDA_etoh_ab, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_etoh_ab <- perf(PLSDA_etoh_ab, validation = "loo") 
plot(perf_plsda_etoh_ab, legend = FALSE)

VIPs_etoh_ab <- as.data.frame(mixOmics::vip(PLSDA_etoh_ab))
VIPs_etoh_ab_filter <- dplyr::filter(VIPs_etoh_ab, VIPs_etoh_ab$comp1 > 1)
VIPs_etoh_ab_filter$ID <- rownames(VIPs_etoh_ab_filter)
VIPs_etoh_ab_select <- VIPs_etoh_ab_filter %>% dplyr::select(ID, comp1)
VIPs_etoh_ab_Load <- VIPs_etoh_ab_select %>% 
  left_join(Loadings_etoh_ab, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

# 95% Etoh - Subject B v Subject C
etoh_BC <- data_20mg %>% left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  dplyr::filter(Donor != "DonorA") %>% dplyr::select(-Donor) %>%  
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_etoh_BC <- mixOmics::pca(etoh_BC, ncomp = 2, center = TRUE, scale = TRUE)
PCA_etoh_BC_scores <- data.frame(PCA_etoh_BC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_etoh_BC_plot <- PCA_etoh_BC_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 95% Etoh - Subject B v Subject C"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_etoh_BC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_etoh_BC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_etoh_BC_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_etoh_BC <- mixOmics::plsda(etoh_BC, PCA_etoh_BC_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_etoh_BC_scores <- data.frame(PLSDA_etoh_BC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_etoh_bc_plot <- PLSDA_etoh_BC_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 95% Etoh - Subject B v Subject C", 
            xlab = paste("Component 1 (", round(PLSDA_etoh_BC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_etoh_BC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_etoh_BC_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_etoh_bc <- plotLoadings(PLSDA_etoh_BC, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_etoh_bc <- perf(PLSDA_etoh_BC, validation = "loo") 
plot(perf_plsda_etoh_bc, legend = FALSE)

VIPs_etoh_bc <- as.data.frame(mixOmics::vip(PLSDA_etoh_BC))
VIPs_etoh_bc_filter <- dplyr::filter(VIPs_etoh_bc, VIPs_etoh_bc$comp1 > 1)
VIPs_etoh_bc_filter$ID <- rownames(VIPs_etoh_bc_filter)
VIPs_etoh_bc_select <- VIPs_etoh_bc_filter %>% dplyr::select(ID, comp1)
VIPs_etoh_bc_Load <- VIPs_etoh_bc_select %>% 
  left_join(Loadings_etoh_bc, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

# 95% Etoh - Subject A v Subject C
etoh_AC <- data_20mg %>% left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  dplyr::filter(Donor != "DonorB") %>% dplyr::select(-Donor) %>%
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  column_to_rownames("SampleID") %>% decostand(method = "rclr")

# PCA
PCA_etoh_AC <- mixOmics::pca(etoh_AC, ncomp = 2, center = TRUE, scale = TRUE)
PCA_etoh_AC_scores <- data.frame(PCA_etoh_AC$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"

PCA_etoh_AC_plot <- PCA_etoh_AC_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, alpha = 0.6,
            title = paste("PCA - 95% Etoh - Subject A v Subject C"), palette = "npg",
            xlab = paste("PC1 (", round(PCA_etoh_AC$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_etoh_AC$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_etoh_AC_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

# PLSDA
PLSDA_etoh_ac <- mixOmics::plsda(etoh_AC, PCA_etoh_AC_scores$Donor, ncomp = 2, scale = TRUE)
PLSDA_etoh_ac_scores <- data.frame(PLSDA_etoh_ac$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PLSDA_etoh_ac_plot <- PLSDA_etoh_ac_scores %>%
  ggscatter(x = "comp1", y = "comp2", color = "Donor", alpha = 0.6, title = "PLS-DA 95% Etoh - Subject A v Subject C", 
            xlab = paste("Component 1 (", round(PLSDA_etoh_ac$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("Component 2 (", round(PLSDA_etoh_ac$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            legend.title = "Group", ggtheme = theme_classic(), palette = c("#F8766D", "#619CFF")) +
  geom_point(data = PLSDA_etoh_ac_scores %>% group_by(Donor) %>% 
               summarise_at(vars(matches("comp")), mean), aes(comp1, comp2, color = Donor), size = 3, shape = 8) +
  theme(plot.title = element_text(size = 10),axis.title = element_text(size = 8),
        axis.text = element_text(size = 5)) + coord_fixed()

Loadings_etoh_ac <- plotLoadings(PLSDA_etoh_ac, plot = FALSE, contrib = "max") %>%
  rownames_to_column() %>% dplyr::select(rowname, GroupContrib, importance)

perf_plsda_etoh_ac <- perf(PLSDA_etoh_ac, validation = "loo") 
plot(perf_plsda_etoh_ac, legend = FALSE)

VIPs_etoh_ac <- as.data.frame(mixOmics::vip(PLSDA_etoh_ac))
VIPs_etoh_ac_filter <- dplyr::filter(VIPs_etoh_ac, VIPs_etoh_ac$comp1 > 1)
VIPs_etoh_ac_filter$ID <- rownames(VIPs_etoh_ac_filter)
VIPs_etoh_ac_select <- VIPs_etoh_ac_filter %>% dplyr::select(ID, comp1)
VIPs_etoh_ac_Load <- VIPs_etoh_ac_select %>% 
  left_join(Loadings_etoh_ac, by = c("ID" = "rowname")) %>% arrange(desc(comp1)) %>%
  left_join(info_feature_complete, by = c("ID" = "Feature")) %>% dplyr::select(1:8)

combined_pca_etoh <- ggarrange(PCA_etoh_AB_plot, PCA_etoh_BC_plot, PCA_etoh_AC_plot, ncol = 3)
combined_plsda_etoh <- ggarrange(PLSDA_etoh_ab_plot, PLSDA_etoh_bc_plot, PLSDA_etoh_ac_plot, ncol = 3)

combined_dim_plot <- ggarrange(combined_pca_meoh, combined_pca_etoh, 
                               combined_plsda_meoh, combined_plsda_etoh, ncol = 1)

#ggsave(plot = combined_dim_plot, filename = "SupplementaryFigure1.svg", device = "svg", dpi = "retina", height = 9, width = 8)


# Check concordance between the different models
colnames(VIPs_etoh_ab_Load)[1:3] <- c("ID", "comp1_etoh", "GroupContrib_etoh")

concordance_AB <- VIPs_meoh_ab_Load %>% 
  inner_join(VIPs_etoh_ab_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>%
  dplyr::mutate(Concordance = GroupContrib == GroupContrib_etoh,
                Delta_VIP = comp1 - comp1_etoh)

colnames(VIPs_etoh_bc_Load)[1:3] <- c("ID", "comp1_etoh", "GroupContrib_etoh")

concordance_BC <- VIPs_meoh_bc_Load %>% 
  inner_join(VIPs_etoh_bc_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>%
  dplyr::mutate(Concordance = GroupContrib == GroupContrib_etoh,
                Delta_VIP = comp1 - comp1_etoh)

colnames(VIPs_etoh_ac_Load)[1:3] <- c("ID", "comp1_etoh", "GroupContrib_etoh")

concordance_AC <- VIPs_meoh_ac_Load %>% 
  inner_join(VIPs_etoh_ac_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>%
  dplyr::mutate(Concordance = GroupContrib == GroupContrib_etoh,
                Delta_VIP = comp1 - comp1_etoh)

conc_ab <- concordance_AB %>% dplyr::filter(Concordance == "TRUE") %>% 
  ggscatter("comp1", "comp1_etoh", add = "reg.line", alpha = 0.4) + 
  stat_cor(label.x = 1, label.y = 1.46) +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject A v Subject B") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

conc_bc <-concordance_BC %>% dplyr::filter(Concordance == "TRUE") %>% 
  ggscatter("comp1", "comp1_etoh", add = "reg.line", alpha = 0.4) + 
  stat_cor(label.x = 1, label.y = 1.475) +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject A v Subject B") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

conc_ac <- concordance_AC %>% dplyr::filter(Concordance == "TRUE") %>% 
  ggscatter("comp1", "comp1_etoh", add = "reg.line", alpha = 0.4) + 
  stat_cor(label.x = 1, label.y = 1.37) +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject A v Subject B") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

conc_comb <- ggarrange(conc_ab, conc_bc, conc_ac, ncol = 3, 
                       common.legend = TRUE, legend = "right")

#ggsave(plot = conc_comb, filename = "SupplementaryFigure2.svg", device = "svg", dpi = "retina", width = 8, height = 4)


# Generate HDRs plots
prob_ab <- concordance_AB %>%
  dplyr::filter(Concordance == "TRUE") %>%
  ggplot(aes(comp1, comp1_etoh)) +
  geom_hdr(aes(fill = after_stat(probs)), 
           alpha = 1) +
  xlim(c(1, 1.4)) +
  ylim(c(1, 1.45)) +
  coord_fixed() +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject A v Subject B") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

prob_bc <- concordance_BC %>%
  dplyr::filter(Concordance == "TRUE") %>%
  ggplot(aes(comp1, comp1_etoh)) +
  geom_hdr(aes(fill = after_stat(probs)), 
           alpha = 1) +
  xlim(c(1, 1.4)) +
  ylim(c(1, 1.46)) +
  coord_fixed() +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject B v Subject C") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

prob_ac <- concordance_AC %>%
  dplyr::filter(Concordance == "TRUE") %>%
  ggplot(aes(comp1, comp1_etoh)) +
  geom_hdr(aes(fill = after_stat(probs)), 
           alpha = 1) +
  xlim(c(1, 1.31)) +
  ylim(c(1, 1.351)) +
  coord_fixed() +
  labs(x = "VIPs MeOH", y = "VIPs EtOH", title = "PLS-DA - Subject A v Subject C") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

prob_comb <- ggarrange(prob_ab, prob_bc, prob_ac, ncol = 3, 
                       common.legend = TRUE, legend = "right")

#ggsave(plot = prob_comb, filename = "Figure2C.svg", device = "svg", dpi = "retina", width = 8, height = 4)


# Check how many of the top100 features detected via MeOH are recovered via EtOH
VIPs_meoh_ab_Load_100 <- VIPs_meoh_ab_Load %>% 
  arrange(desc(comp1)) %>% head(100) %>% 
  left_join(VIPs_etoh_ab_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>% 
  dplyr::filter(is.na(comp1_etoh)) %>%
  left_join(intersect_extractions, by = c("ID" = "Feature"))

VIPs_meoh_bc_Load_100 <- VIPs_meoh_bc_Load %>% 
  arrange(desc(comp1)) %>% head(100) %>% 
  left_join(VIPs_etoh_bc_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>% 
  dplyr::filter(is.na(comp1_etoh)) %>%
  left_join(intersect_extractions, by = c("ID" = "Feature"))

VIPs_meoh_ac_Load_100 <- VIPs_meoh_ac_Load %>% 
  arrange(desc(comp1)) %>% head(100) %>% 
  left_join(VIPs_etoh_ac_Load %>% dplyr::select(ID, comp1_etoh, GroupContrib_etoh)) %>% 
  dplyr::filter(is.na(comp1_etoh)) %>%
  left_join(intersect_extractions, by = c("ID" = "Feature"))


# Check correlation of fold changes in univariate analysis
meoh_ab_ra <- data_20mg %>%
  dplyr::filter(str_detect(pattern = "A|B", SampleID)) %>% 
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_meoh_ab <- meoh_ab_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorA", "DonorB"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorB/DonorA) %>%
  mutate(Log2FC_meoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_meoh)) %>% dplyr::select(Metabolite, Log2FC_meoh)

meoh_bc_ra <- data_20mg %>%
  dplyr::filter(str_detect(pattern = "B|C", SampleID)) %>% 
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_meoh_bc <- meoh_bc_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorB", "DonorC"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorC/DonorB) %>%
  mutate(Log2FC_meoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_meoh)) %>% dplyr::select(Metabolite, Log2FC_meoh)

meoh_ac_ra <- data_20mg %>%
  dplyr::filter(str_detect(pattern = "A|C", SampleID)) %>% 
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_meoh_ac <- meoh_ac_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorA", "DonorC"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorC/DonorA) %>%
  mutate(Log2FC_meoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_meoh)) %>% dplyr::select(Metabolite, Log2FC_meoh)

etoh_ab_ra <- data_20mg %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(str_detect(pattern = "A|B", Donor)) %>% dplyr::select(-Donor) %>%
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_etoh_ab <- etoh_ab_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorA", "DonorB"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorB/DonorA) %>%
  mutate(Log2FC_etoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_etoh)) %>% dplyr::select(Metabolite, Log2FC_etoh)

etoh_bc_ra <- data_20mg %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(str_detect(pattern = "B|C", Donor)) %>% dplyr::select(-Donor) %>%
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_etoh_bc <- etoh_bc_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorB", "DonorC"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorC/DonorB) %>%
  mutate(Log2FC_etoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_etoh)) %>% dplyr::select(Metabolite, Log2FC_etoh)

etoh_ac_ra <- data_20mg %>%
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>%
  dplyr::filter(str_detect(pattern = "A|C", Donor)) %>% dplyr::select(-Donor) %>%
  column_to_rownames("SampleID") %>% select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Donor)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor) %>%
  summarise(across(everything(), mean))

fold_change_etoh_ac <- etoh_ac_ra %>%
  pivot_longer(cols = -Donor, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor, values_from = Value) %>% 
  mutate(across(c("DonorA", "DonorC"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorC/DonorA) %>%
  mutate(Log2FC_etoh = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_etoh)) %>% dplyr::select(Metabolite, Log2FC_etoh)

# Plots
fold_ab <- fold_change_meoh_ab %>% inner_join(fold_change_etoh_ab) %>% 
  dplyr::filter(Metabolite %in% concordance_AB$ID) %>%
  ggscatter(x = "Log2FC_meoh", y = "Log2FC_etoh",
            add = "reg.line", alpha = 0.5) + stat_cor() + coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log2FC MeOH", y = "Log2FC EtOH", title = "Fold Change - Subject A v Subject B") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

fold_bc <- fold_change_meoh_bc %>% inner_join(fold_change_etoh_bc) %>% 
  dplyr::filter(Metabolite %in% concordance_BC$ID) %>% 
  ggscatter(x = "Log2FC_meoh", y = "Log2FC_etoh",
            add = "reg.line", alpha = 0.5) + stat_cor() + coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log2FC MeOH", y = "Log2FC EtOH", title = "Fold Change - Subject B v Subject C") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

fold_ac <- fold_change_meoh_ac %>% inner_join(fold_change_etoh_ac) %>% 
  dplyr::filter(Metabolite %in% concordance_AC$ID) %>%
  ggscatter(x = "Log2FC_meoh", y = "Log2FC_etoh",
            add = "reg.line", alpha = 0.5) + stat_cor() + coord_fixed() +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "gray") +
  labs(x = "Log2FC MeOH", y = "Log2FC EtOH", title = "Fold Change - Subject A v Subject C") +
  theme(plot.title = element_text(size = 8),  axis.title.x = element_text(size = 8),  
        axis.title.y = element_text(size = 8), axis.text.x = element_text(size = 6),
        axis.text.y = element_text(size = 6))

fold_combined <- ggarrange(fold_ab, fold_bc, fold_ac, ncol = 3)

#ggsave(plot = fold_combined, filename = "Figure2D.svg", device = "svg", dpi = "retina", width = 8, height = 4)


######################
# Figure 3 - Storage #
######################

# Extract 20mg triplicates for storage analysis
sample_store <- metadata_final %>% 
  dplyr::filter(Weight == "20mg") %>% 
  dplyr::filter(Time == "im_2mn" | Time == "1dy" | Time == "1wk")

data_store <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_store$SampleID) %>% 
  column_to_rownames("SampleID") %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  rownames_to_column("SampleID")

# RCLR transformation
data_store_clr <- decostand(data_store %>% column_to_rownames("SampleID"), method = "rclr")

# PCA of 20mg
PCA_store <- mixOmics::pca(data_store_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_store_scores <- data.frame(PCA_store$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Time"
j <- "Donor"

PCA_store_plot <- PCA_store_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = j, shape = i, alpha = 0.6,
            title = paste("PCA - Human Fecal Metabolome"),
            xlab = paste("PC1 (", round(PCA_store$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_store$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_store_scores %>% group_by((!!sym(j))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(j))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed() + scale_color_viridis_d(option = "B")

# PERMANOVA
dist_metabolites <- vegdist(data_store_clr, method = "euclidean")
disper_donor <- betadisper(dist_metabolites, PCA_store_scores$Donor)
anova(disper_donor)
permanova <- adonis2(dist_metabolites ~ Donor * Time, PCA_store_scores, na.action = na.omit)

#ggsave(plot = PCA_store_plot, filename = "Figure3A.svg", device = "svg", dpi = "retina", width = 4)


# Generate Upset plot to check metabolic features recovered by storage time
data_storage_presence <- data_store %>% column_to_rownames("SampleID")
data_storage_presence[data_storage_presence > 0] <- 1
data_storage_presence_info <- data_storage_presence %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Time)) %>% column_to_rownames("SampleID")

list_storage <- list(`Immediate` = data_storage_presence_info %>% dplyr::filter(Time == "im_2mn") %>% 
                       dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% colnames(),
                     `1 Day` = data_storage_presence_info %>% dplyr::filter(Time == "1dy") %>% 
                       dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% colnames(),
                     `1 Week` = data_storage_presence_info %>% dplyr::filter(Time == "1wk") %>% 
                       dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% colnames())

upset_storage <- upset(fromList(list_storage), order.by = "freq") # save pdf for Figure 3B


# Identified features exclusively detected at different timepoints of storage
features_imm <- data_storage_presence_info %>% dplyr::filter(Time == "im_2mn") %>%
  dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

features_day <- data_storage_presence_info %>% dplyr::filter(Time == "1dy") %>%
  dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

features_week <- data_storage_presence_info %>% dplyr::filter(Time == "1wk") %>%
  dplyr::select(-Time) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

all_time <- features_imm %>% full_join(features_day) %>% full_join(features_week) %>%
  left_join(info_feature_complete) %>% distinct(Feature, .keep_all = TRUE) %>%
  dplyr::filter(!(is.na(Compound_Name)))

intersect_time <- features_imm %>% inner_join(features_day) %>% 
  inner_join(features_week) %>% 
  left_join(info_feature_complete)

exclusive_imm <- features_imm %>% anti_join(features_day) %>% 
  anti_join(features_week) %>% left_join(info_feature_complete)
exclusive_day <- features_day %>% anti_join(features_imm) %>% 
  anti_join(features_week) %>% left_join(info_feature_complete)
exclusive_week <- features_week %>% anti_join(features_imm) %>% 
  anti_join(features_day) %>% left_join(info_feature_complete)

intersect_imm_day <- features_imm %>% inner_join(features_day) %>%
  anti_join(features_week) %>% left_join(info_feature_complete)
intersect_imm_week <- features_imm %>% inner_join(features_week) %>%
  anti_join(features_day) %>% left_join(info_feature_complete)
intersect_day_week <- features_day %>% inner_join(features_week) %>%
  anti_join(features_imm) %>% left_join(info_feature_complete)

# Superclass prediction via CANOPUS
exclusive_imm_pred <- exclusive_imm %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_exclusive")

colnames(exclusive_imm_pred)[1] <- "NPC Superclass"

exclusive_day_pred <- exclusive_day %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Day_exclusive")

colnames(exclusive_day_pred)[1] <- "NPC Superclass"

exclusive_week_pred <- exclusive_week %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Week_exclusive")

colnames(exclusive_week_pred)[1] <- "NPC Superclass"

imm_day_pred <- intersect_imm_day %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_day")

colnames(imm_day_pred)[1] <- "NPC Superclass"

imm_week_pred <- intersect_imm_week %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_week")

colnames(imm_week_pred)[1] <- "NPC Superclass"

day_week_pred <- intersect_day_week %>% group_by(`NPC#superclass`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Day_week")

colnames(day_week_pred)[1] <- "NPC Superclass"

pred_combine <- rbind(exclusive_imm_pred, exclusive_day_pred, exclusive_week_pred,
                      imm_day_pred, imm_week_pred, day_week_pred)

#write_csv(x = pred_combine, file = "SupplementaryTable2.csv")

# Pie chart for NPC Pathway
exclusive_imm_pred <- exclusive_imm %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_exclusive")

colnames(exclusive_imm_pred)[1] <- "NPC Pathway"

pie_imm <- exclusive_imm_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

exclusive_day_pred <- exclusive_day %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Day_exclusive")

colnames(exclusive_day_pred)[1] <- "NPC Pathway"

pie_day <- exclusive_day_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

exclusive_week_pred <- exclusive_week %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Week_exclusive")

colnames(exclusive_week_pred)[1] <- "NPC Pathway"

pie_week <- exclusive_week_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

imm_day_pred <- intersect_imm_day %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_day")

colnames(imm_day_pred)[1] <- "NPC Pathway"

pie_imm_day <- imm_day_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

imm_week_pred <- intersect_imm_week %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Immediate_week")

colnames(imm_week_pred)[1] <- "NPC Pathway"

pie_imm_week <- imm_week_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

day_week_pred <- intersect_day_week %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count)) %>%
  dplyr::mutate(Time = "Day_week")

colnames(day_week_pred)[1] <- "NPC Pathway"

pie_day_week <- day_week_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

pie_combined_storage <- ggarrange(pie_imm_day, pie_week, pie_day_week, pie_imm_week,
                                  pie_imm, pie_day, common.legend = TRUE, legend = "right", nrow = 1)

#ggsave(plot = pie_combined_storage, filename = "Figure3B_pies.svg", device = "svg")


# Generate 36 pairwise PLS-DA models
metadata_storage <- metadata_final %>% dplyr::filter(Weight == "20mg") %>%
  dplyr::filter(!(str_detect(pattern = "DLAB|4mn", Time))) %>% 
  dplyr::mutate(Donor_time = paste(Donor, Time, sep = "_")) %>% 
  dplyr::select(SampleID, Donor_time)

data_store_info <- data_store %>% left_join(metadata_storage) %>% 
  column_to_rownames("SampleID")

# Check pairwise models 

# Function
generate_models <- function(data, type_column) {
  
  # Get all unique pairs of sample types
  sample_types <- unique(data[[type_column]])
  sample_pairs <- combn(sample_types, 2, simplify = FALSE)
  
  pca_plots <- list()
  plsda_plots <- list()
  vip_combined <- data.frame()
  error_plsda_combined <- data.frame()
  
  for (pair in sample_pairs) {
    
    # Filter data for the current pair of sample types
    data_filtered <- data %>%
      dplyr::filter((!!sym(type_column) == pair[1]) | (!!sym(type_column) == pair[2])) %>%
      dplyr::select(-Donor_time)
    
    # PCA
    pca_model <- mixOmics::pca(data_filtered %>%
                                 select_at(vars(-one_of(nearZeroVar(., names = TRUE)))), 
                               ncomp = 2, center = TRUE, scale = TRUE)
    
    pca_scores <- data.frame(pca_model$variates$X) %>%
      rownames_to_column("SampleID") %>%
      left_join(metadata_storage)
    
    pca_plot <- pca_scores %>%
      ggscatter(x = "PC1", y = "PC2", color = type_column, alpha = 0.6,
                title = paste("PCA -", pair[1], "vs", pair[2]),
                xlab = paste("PC1 (", round(pca_model$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
                ylab = paste("PC2 (", round(pca_model$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
                ggtheme = theme_classic()) +
      geom_point(data = pca_scores %>%
                   group_by(!!sym(type_column)) %>%
                   summarise(across(starts_with("PC"), ~ mean(as.numeric(.), na.rm = TRUE))),
                 aes(PC1, PC2, color = !!sym(type_column)), size = 4, shape = 8) +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 4)) +
      coord_fixed()
    
    pca_plots[[paste(pair[1], pair[2], sep = "_")]] <- pca_plot
    
    # PLSDA
    metadata_filter <- data_frame(SampleID = rownames(data_filtered)) %>% left_join(metadata_storage)
    
    plsda_model <- mixOmics::plsda(data_filtered %>%
                                     select_at(vars(-one_of(nearZeroVar(., names = TRUE)))),
                                   metadata_filter[[type_column]], ncomp = 2, scale = TRUE)
    
    plsda_scores <- data.frame(plsda_model$variates$X) %>%
      rownames_to_column("SampleID") %>%
      left_join(metadata_storage)
    
    plsda_plot <- plsda_scores %>%
      ggscatter(x = "comp1", y = "comp2", color = type_column, alpha = 0.6,
                title = paste("PLSDA -", pair[1], "vs", pair[2]),
                xlab = paste("Component 1 (", round(plsda_model$prop_expl_var$X[1] * 100, digits = 1), "%)", sep = ""),
                ylab = paste("Component 2 (", round(plsda_model$prop_expl_var$X[2] * 100, digits = 1), "%)", sep = ""),
                legend.title = "Group", ggtheme = theme_classic()) +
      geom_point(data = plsda_scores %>%
                   select(comp1, comp2, !!sym(type_column)) %>%
                   group_by(!!sym(type_column)) %>%
                   summarise(across(matches("comp"), mean)),
                 aes(comp1, comp2, color = !!sym(type_column)), size = 3, shape = 8) +
      theme(plot.title = element_text(size = 10),
            axis.title = element_text(size = 8),
            axis.text = element_text(size = 5)) +
      coord_fixed()
    
    plsda_loadings <- plotLoadings(plsda_model, plot = FALSE, contrib = "max") %>%
      rownames_to_column() %>% dplyr::select(rowname, GroupContrib)
    
    plsda_plots[[paste(pair[1], pair[2], sep = "_")]] <- plsda_plot
    
    # Generate VIP data frames
    vip_data <- as.data.frame(mixOmics::vip(plsda_model))
    vip_filtered <- dplyr::filter(vip_data, vip_data$comp1 > 1)
    vip_filtered$ID <- rownames(vip_filtered)
    vip_select <- vip_filtered %>%
      dplyr::select(ID, comp1) %>% dplyr::arrange(desc(comp1)) %>%
      left_join(plsda_loadings, by = c("ID" = "rowname")) %>% 
      left_join(info_feature_complete, by = c("ID" = "Feature"))
    
    vip_select <- vip_select %>%
      dplyr::mutate(model_pair = paste(pair[1], pair[2], sep = "_"))
    
    vip_combined <- rbind(vip_combined, vip_select)
    
    # Evaluate model
    perf_plsda <- perf(plsda_model, validation = "loo", progressBar = FALSE, auc = FALSE)
    error_plsda <- data_frame(CER = min(perf_plsda$error.rate$BER[2,])) %>%
      dplyr::mutate(model_pair = paste(pair[1], pair[2], sep = ";"))
    
    error_plsda_combined <- rbind(error_plsda_combined, error_plsda)
  }
  
  # Combine all PCA plots into one
  pca_combined_plot <- plot_grid(plotlist = pca_plots, ncol = 10)
  
  # Combine all PLSDA plots into one
  plsda_combined_plot <- plot_grid(plotlist = plsda_plots, ncol = 10)
  
  assign("VIP_combined", vip_combined, envir = .GlobalEnv)
  assign("Error_PLSDA_combined", error_plsda_combined, envir = .GlobalEnv)
  assign("PCA_combined_plot", pca_combined_plot, envir = .GlobalEnv)
  assign("PLSDA_combined_plot", plsda_combined_plot, envir = .GlobalEnv)
  
}

# Example usage
generate_models(data_store_info, "Donor_time")

# Extacrt VIP of reference and compare
VIP_AB <- VIP_combined %>% 
  dplyr::filter(str_detect(pattern = "DonorA", model_pair)) %>% 
  dplyr::filter(str_detect(pattern = "DonorB", model_pair))

VIP_BC <- VIP_combined %>% 
  dplyr::filter(str_detect(pattern = "DonorB", model_pair)) %>% 
  dplyr::filter(str_detect(pattern = "DonorC", model_pair))

VIP_AC <- VIP_combined %>% 
  dplyr::filter(str_detect(pattern = "DonorA", model_pair)) %>% 
  dplyr::filter(str_detect(pattern = "DonorC", model_pair))

# Select top 100 features for reference
VIP_AB_ref <- VIP_AB %>% 
  dplyr::filter(model_pair == "DonorB_im_2mn_DonorA_im_2mn") %>%
  head(100)
VIP_BC_ref <- VIP_BC %>% 
  dplyr::filter(model_pair == "DonorC_im_2mn_DonorB_im_2mn") %>%
  head(100)
VIP_AC_ref <- VIP_AC %>% 
  dplyr::filter(model_pair == "DonorC_im_2mn_DonorA_im_2mn") %>%
  head(100)


# Generate function
process_model_pairs <- function(VIP, VIP_ref) {
  
  unique_pairs <- unique(VIP$model_pair)
  
  result_df <- data.frame(ModelPair = character(),
                          RowCount = integer(),
                          stringsAsFactors = FALSE)
  
  for (pair in unique_pairs) {
    
    VIP_filter <- VIP %>% 
      dplyr::filter(model_pair == pair) %>%
      dplyr::filter(ID %in% VIP_ref$ID) %>%
      left_join(VIP_ref %>% dplyr::select(ID, GroupContrib) %>% 
                  rename(GroupContrib_ref = GroupContrib)) %>%
      dplyr::mutate(GroupContrib = gsub("_.*", "", GroupContrib)) %>%
      dplyr::mutate(GroupContrib_ref = gsub("_.*", "", GroupContrib_ref)) %>%
      dplyr::mutate(Concordance = GroupContrib == GroupContrib_ref) %>%
      dplyr::filter(Concordance == TRUE)
    
    result_df <- rbind(result_df, data.frame(ModelPair = pair, RowCount = nrow(VIP_filter)))
  }
  
  return(result_df)
}

# Example usage
results_AB <- process_model_pairs(VIP_AB, VIP_AB_ref) %>% 
  arrange(desc(RowCount)) %>% dplyr::mutate(Ratio = RowCount/100)
results_BC <- process_model_pairs(VIP_BC, VIP_BC_ref) %>% 
  arrange(desc(RowCount)) %>% dplyr::mutate(Ratio = RowCount/100)
results_AC <- process_model_pairs(VIP_AC, VIP_AC_ref) %>% 
  arrange(desc(RowCount)) %>% dplyr::mutate(Ratio = RowCount/100)

lollipop_AB <- results_AB %>% ggdotchart(x = "ModelPair", y = "RowCount", dot.size = 3, color = "#050507", sorting = "none",
                                         order = c('DonorB_im_2mn_DonorA_im_2mn', 'DonorB_1dy_DonorA_im_2mn', 'DonorA_1dy_DonorB_im_2mn', 
                                                   'DonorA_1dy_DonorB_1dy', 'DonorA_1wk_DonorB_im_2mn', 'DonorB_1wk_DonorA_im_2mn', 
                                                   'DonorA_1wk_DonorB_1dy', 'DonorB_1wk_DonorA_1dy', 'DonorA_1wk_DonorB_1wk'),
                                         add = "segments", add.params = list(color = "lightgray", size = 0.1), ylab = "Recovery (%)",
                                         title = "PLSDA Storage - Subject A v Subject B") + ylim(80, 100) +
  theme_bw() + theme(plot.title = element_text(size = 10),
                     axis.title = element_text(size = 8),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                     axis.text = element_text(size = 6),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

lollipop_BC <- results_BC %>% ggdotchart(x = "ModelPair", y = "RowCount", dot.size = 3, color = "#ba3955", sorting = "none",
                                         order = c('DonorC_im_2mn_DonorB_im_2mn', 'DonorC_1dy_DonorB_im_2mn', 'DonorB_1dy_DonorC_im_2mn',
                                                   'DonorB_1dy_DonorC_1dy', 'DonorC_1wk_DonorB_im_2mn', 'DonorB_1wk_DonorC_im_2mn', 
                                                   'DonorC_1wk_DonorB_1dy', 'DonorB_1wk_DonorC_1dy', 'DonorB_1wk_DonorC_1wk'),
                                         add = "segments", add.params = list(color = "lightgray", size = 0.1), ylab = "Recovery (%)",
                                         title = "PLSDA Storage - Subject B v Subject C") + ylim(80, 100) +
  theme_bw() + theme(plot.title = element_text(size = 10),
                     axis.title = element_text(size = 8),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                     axis.text = element_text(size = 6),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

lollipop_AC <- results_AC %>% ggdotchart(x = "ModelPair", y = "RowCount", dot.size = 3, color = "#fbb040", sorting = "none",
                                         order = c('DonorC_im_2mn_DonorA_im_2mn', 'DonorC_1dy_DonorA_im_2mn', 'DonorA_1dy_DonorC_im_2mn', 
                                                   'DonorA_1dy_DonorC_1dy', 'DonorA_1wk_DonorC_im_2mn', 'DonorC_1wk_DonorA_im_2mn', 
                                                   'DonorA_1wk_DonorC_1dy', 'DonorC_1wk_DonorA_1dy', 'DonorA_1wk_DonorC_1wk'),
                                         add = "segments", add.params = list(color = "lightgray", size = 0.1), ylab = "Recovery (%)",
                                         title = "PLSDA Storage - Subject A v Subject C") + ylim(80, 100) +
  theme_bw() + theme(plot.title = element_text(size = 10),
                     axis.title = element_text(size = 8),
                     axis.title.x = element_blank(),
                     axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
                     axis.text = element_text(size = 6),
                     panel.grid.major = element_blank(), 
                     panel.grid.minor = element_blank())

lollipop_combined <- ggarrange(lollipop_AB, lollipop_BC, lollipop_AC, nrow = 1)

#ggsave(plot = lollipop_combined, filename = "Figure3C.svg", device = "svg", dpi = "retina", width = 7.5, height = 3)


# Calculate fold changes of features per subject after 1 week storage
subA_week_ra <- data_store_info %>%
  dplyr::filter(str_detect(pattern = "DonorA_im_2mn|DonorA_1wk", Donor_time)) %>% 
  dplyr::select(-Donor_time) %>% dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_storage %>% dplyr::select(SampleID, Donor_time)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor_time) %>%
  summarise(across(everything(), mean))

fc_subA_week <- subA_week_ra %>%
  pivot_longer(cols = -Donor_time, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor_time, values_from = Value) %>% 
  #mutate(across(c("DonorA_1wk", "DonorA_im_2mn"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorA_1wk/DonorA_im_2mn) %>%
  mutate(Log2FC_A = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_A)) %>% dplyr::select(Metabolite, Log2FC_A)

subB_week_ra <- data_store_info %>%
  dplyr::filter(str_detect(pattern = "DonorB_im_2mn|DonorB_1wk", Donor_time)) %>% 
  dplyr::select(-Donor_time) %>% dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_storage %>% dplyr::select(SampleID, Donor_time)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor_time) %>%
  summarise(across(everything(), mean))

fc_subB_week <- subB_week_ra %>%
  pivot_longer(cols = -Donor_time, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor_time, values_from = Value) %>% 
  #mutate(across(c("DonorB_1wk", "DonorB_im_2mn"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorB_1wk/DonorB_im_2mn) %>%
  mutate(Log2FC_B = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_B)) %>% dplyr::select(Metabolite, Log2FC_B)

subC_week_ra <- data_store_info %>%
  dplyr::filter(str_detect(pattern = "DonorC_im_2mn|DonorC_1wk", Donor_time)) %>% 
  dplyr::select(-Donor_time) %>% dplyr::select(where(~ sum(.) != 0)) %>%
  dplyr::mutate(row_sum = rowSums(.)) %>%
  dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
  dplyr::select(-row_sum) %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_storage %>% dplyr::select(SampleID, Donor_time)) %>% 
  dplyr::select(-SampleID) %>% group_by(Donor_time) %>%
  summarise(across(everything(), mean))

fc_subC_week <- subC_week_ra %>%
  pivot_longer(cols = -Donor_time, names_to = "Metabolite", values_to = "Value") %>%
  pivot_wider(names_from = Donor_time, values_from = Value) %>% 
  #mutate(across(c("DonorC_1wk", "DonorC_im_2mn"), ~ifelse(. == 0, 1e-9, .))) %>%
  mutate(Fold_Change = DonorC_1wk/DonorC_im_2mn) %>%
  mutate(Log2FC_C = log2(Fold_Change)) %>%
  arrange(desc(Log2FC_C)) %>% dplyr::select(Metabolite, Log2FC_C)

# Check combined FC and annotations
fc_sub_combined <- fc_subA_week %>% inner_join(fc_subB_week) %>%
  inner_join(fc_subC_week) %>% 
  dplyr::mutate(FC_time = (rowSums(cbind(Log2FC_A, Log2FC_B, Log2FC_C)))/3) %>%
  dplyr::mutate(Concordance = case_when(Log2FC_A > 0 & Log2FC_B > 0 & Log2FC_C > 0 ~ "Yes",
                                        Log2FC_A < 0 & Log2FC_B < 0 & Log2FC_C < 0 ~ "Yes",
                                        TRUE ~ "No")) %>%
  dplyr::mutate(Group = "All")

fc_sub_AB <- fc_subA_week %>% full_join(fc_subB_week) %>% 
  dplyr::filter(Metabolite %in% VIP_AB_ref$ID) %>% 
  dplyr::mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  dplyr::mutate(FC_time = case_when(is.na(Log2FC_A) ~ Log2FC_B,
                                    is.na(Log2FC_B) ~ Log2FC_A,
                                    TRUE ~ rowMeans(cbind(abs(Log2FC_A), abs(Log2FC_B)), na.rm = TRUE))) %>%
  dplyr::mutate(Group = "AB")

fc_sub_BC <- fc_subB_week %>% full_join(fc_subC_week) %>% 
  dplyr::filter(Metabolite %in% VIP_BC_ref$ID) %>% 
  dplyr::mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  dplyr::mutate(FC_time = case_when(is.na(Log2FC_B) ~ Log2FC_C,
                                    is.na(Log2FC_C) ~ Log2FC_B,
                                    TRUE ~ rowMeans(cbind(abs(Log2FC_B), abs(Log2FC_C)), na.rm = TRUE))) %>%
  dplyr::mutate(Group = "BC")

fc_sub_AC <- fc_subA_week %>% full_join(fc_subC_week) %>% 
  dplyr::filter(Metabolite %in% VIP_AC_ref$ID) %>% 
  dplyr::mutate(across(where(is.numeric), ~ na_if(., Inf))) %>%
  dplyr::mutate(FC_time = case_when(is.na(Log2FC_A) ~ Log2FC_C,
                                    is.na(Log2FC_C) ~ Log2FC_A,
                                    TRUE ~ rowMeans(cbind(abs(Log2FC_A), abs(Log2FC_C)), na.rm = TRUE))) %>%
  dplyr::mutate(Group = "AC")


fc_comparison <- rbind(fc_sub_AB %>% dplyr::select(Group, FC_time), 
                       fc_sub_BC %>% dplyr::select(Group, FC_time),
                       fc_sub_AC %>% dplyr::select(Group, FC_time)) %>%
  dplyr::filter(!(is.infinite(FC_time))) %>%
  dplyr::mutate(FC_time = abs(FC_time))

# Generate boxplot
fc_comparison_plot <- fc_comparison %>% ggboxplot(x = "Group", y = "FC_time", add = "jitter", ylab = "abc(Log2FC)",
                                                  title = "Time Log2FC of Top 100 Features per pairwise PLS-DA model") +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text = element_text(size = 6))

#ggsave(plot = fc_comparison_plot, filename = "SupplementaryFigure3B.svg", device = "svg", dpi = "retina", width = 6, height = 3)

fc_comparison %>% group_by(Group) %>% summarise(mean_group = median(FC_time))


# Check all features beased on classification
fc_subA_week_info <- fc_subA_week %>% 
  dplyr::filter(!is.infinite(Log2FC_A)) %>%
  left_join(info_feature_complete, by = c("Metabolite" = "Feature")) %>%
  rename_with(~ gsub("#", "_", .)) %>%
  dplyr::mutate(Log2FC_A = abs(Log2FC_A)) %>%
  ggboxplot("NPC_pathway", y = "Log2FC_A", add = "jitter", ylab = "abc(Log2FC)",
            title = "Time Log2FC of Subject A Features", add.params = list(alpha = 0.3),
            order = c("Alkaloids", "Amino acids and Peptides", "Carbohydrates",
                      "Fatty acids", "Polyketides", "Shikimates and Phenylpropanoids", 
                      "Terpenoids")) + ylim(0, 11) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text = element_text(size = 6))

fc_subB_week_info <- fc_subB_week %>% 
  dplyr::filter(!is.infinite(Log2FC_B)) %>%
  left_join(info_feature_complete, by = c("Metabolite" = "Feature")) %>%
  rename_with(~ gsub("#", "_", .)) %>%
  dplyr::mutate(Log2FC_B = abs(Log2FC_B)) %>%
  ggboxplot("NPC_pathway", y = "Log2FC_B", add = "jitter", ylab = "abc(Log2FC)",
            title = "Time Log2FC of Subject B Features", add.params = list(alpha = 0.3),
            order = c("Alkaloids", "Amino acids and Peptides", "Carbohydrates",
                      "Fatty acids", "Polyketides", "Shikimates and Phenylpropanoids", 
                      "Terpenoids")) + ylim(0, 11) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text = element_text(size = 6))

fc_subC_week_info <- fc_subC_week %>% 
  dplyr::filter(!is.infinite(Log2FC_C)) %>%
  left_join(info_feature_complete, by = c("Metabolite" = "Feature")) %>%
  rename_with(~ gsub("#", "_", .)) %>%
  dplyr::mutate(Log2FC_C = abs(Log2FC_C)) %>%
  ggboxplot("NPC_pathway", y = "Log2FC_C", add = "jitter", ylab = "abc(Log2FC)",
            title = "Time Log2FC of Subject C Features", add.params = list(alpha = 0.3),
            order = c("Alkaloids", "Amino acids and Peptides", "Carbohydrates",
                      "Fatty acids", "Polyketides", "Shikimates and Phenylpropanoids", 
                      "Terpenoids")) + ylim(0, 11) +
  theme(plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.title.x = element_blank(),
        axis.text.x = element_text(size = 6, angle = 45, hjust = 1),
        axis.text = element_text(size = 6))

fc_sub_week_combined <- ggarrange(fc_subA_week_info, fc_subB_week_info, fc_subC_week_info, nrow = 3)

#ggsave(plot = fc_sub_week_combined, filename = "SupplementaryFigure3A.svg", device = "svg", dpi = "retina", width = 8, height = 8)


#######################
# Figure 4 - Sampling #
#######################

# Extract 10, 20, 30mg triplicates for 95% EtOH and 50% MeOH 
sample_sampling <- metadata_final %>% 
  dplyr::filter(Time == "1wk" | Time == "DLAB")

data_sampling_etoh <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_sampling$SampleID) %>% 
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>%
  column_to_rownames("SampleID") %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  rownames_to_column("SampleID")

data_sampling_meoh <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_sampling$SampleID) %>% 
  dplyr::filter(str_detect(pattern = "A|B|C", SampleID)) %>%
  column_to_rownames("SampleID") %>% 
  select_at(vars(-one_of(nearZeroVar(., names = TRUE)))) %>% 
  rownames_to_column("SampleID")

# RCLR transformation
data_sampling_etoh_clr <- decostand(data_sampling_etoh %>% column_to_rownames("SampleID"), method = "rclr")
data_sampling_meoh_clr <- decostand(data_sampling_meoh %>% column_to_rownames("SampleID"), method = "rclr")

# PCA sampling EtOH and MeOH
PCA_sampling_etoh <- mixOmics::pca(data_sampling_etoh_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_etoh_scores <- data.frame(PCA_sampling_etoh$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PCA_sampling_meoh <- mixOmics::pca(data_sampling_meoh_clr, ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_meoh_scores <- data.frame(PCA_sampling_meoh$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"
j <- "Weight"

# Gonna use Manu color palette ("G-Thomson/Manu")

PCA_sampling_etoh_plot <- PCA_sampling_etoh_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#DD3C51", "#313657", "#1F6683"),
            title = paste("PCA - Human Fecal Metabolome EtOH"),
            xlab = paste("PC1 (", round(PCA_sampling_etoh$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_etoh$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_etoh_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

PCA_sampling_meoh_plot <- PCA_sampling_meoh_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#214d65", "#287DAB", "#E5BF86"),
            title = paste("PCA - Human Fecal Metabolome MeOH"),
            xlab = paste("PC1 (", round(PCA_sampling_meoh$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_meoh$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_meoh_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

pca_combined <- ggarrange(PCA_sampling_etoh_plot, PCA_sampling_meoh_plot, nrow = 1)

#ggsave(plot = pca_combined, filename = "Figure4A.svg", device = "svg", dpi = "retina", width = 8)

# PERMANOVA
dist_metabolites_etoh <- vegdist(data_sampling_etoh_clr, method = "euclidean")
disper_donor_etoh <- betadisper(dist_metabolites_etoh, PCA_sampling_etoh_scores$Donor)
anova(disper_donor_etoh)
permanova_etoh <- adonis2(dist_metabolites_etoh ~ Donor * Weight, 
                          PCA_sampling_etoh_scores, na.action = na.omit)

dist_metabolites_meoh <- vegdist(data_sampling_meoh_clr, method = "euclidean")
disper_donor_meoh <- betadisper(dist_metabolites_meoh, PCA_sampling_meoh_scores$Donor)
anova(disper_donor_meoh)
permanova_meoh <- adonis2(dist_metabolites_meoh ~ Donor * Weight, 
                          PCA_sampling_meoh_scores, na.action = na.omit)

#ggsave(plot = PCA_plot, filename = "Figure2A.svg", device = "svg", dpi = "retina", width = 4)


# Check not transformed data
PCA_sampling_etoh_raw <- mixOmics::pca(data_sampling_etoh %>% 
                                     column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_etoh_raw_scores <- data.frame(PCA_sampling_etoh_raw$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PCA_sampling_meoh_raw <- mixOmics::pca(data_sampling_meoh %>% 
                                     column_to_rownames("SampleID"), ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_meoh_raw_scores <- data.frame(PCA_sampling_meoh_raw$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"
j <- "Weight"

PCA_sampling_etoh_raw_plot <- PCA_sampling_etoh_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#DD3C51", "#313657", "#1F6683"),
            title = paste("PCA - Human Fecal Metabolome EtOH"),
            xlab = paste("PC1 (", round(PCA_sampling_etoh_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_etoh_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_etoh_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

PCA_sampling_meoh_raw_plot <- PCA_sampling_meoh_raw_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#214d65", "#287DAB", "#E5BF86"),
            title = paste("PCA - Human Fecal Metabolome MeOH"),
            xlab = paste("PC1 (", round(PCA_sampling_meoh_raw$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_meoh_raw$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_meoh_raw_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

pca_raw_combined <- ggarrange(PCA_sampling_etoh_raw_plot, PCA_sampling_meoh_raw_plot, nrow = 1)

#ggsave(plot = pca_raw_combined, filename = "SupplementaryFFigure4A.svg", device = "svg", dpi = "retina", width = 8)

# PERMANOVA
dist_metabolites_etoh_raw <- vegdist(data_sampling_etoh %>% column_to_rownames("SampleID"), method = "euclidean")
disper_donor_etoh_raw <- betadisper(dist_metabolites_etoh_raw, PCA_sampling_etoh_raw_scores$Donor)
anova(disper_donor_etoh_raw)
permanova_etoh_raw <- adonis2(dist_metabolites_etoh_raw ~ Donor * Weight, 
                          PCA_sampling_etoh_raw_scores, na.action = na.omit)

dist_metabolites_meoh_raw <- vegdist(data_sampling_meoh %>% column_to_rownames("SampleID"), method = "euclidean")
disper_donor_meoh_raw <- betadisper(dist_metabolites_meoh_raw, PCA_sampling_meoh_raw_scores$Donor)
anova(disper_donor_meoh_raw)
permanova_meoh_raw <- adonis2(dist_metabolites_meoh_raw ~ Donor * Weight, 
                              PCA_sampling_meoh_raw_scores, na.action = na.omit)


# Check relative abundance data
PCA_sampling_etoh_ra <- mixOmics::pca(data_sampling_etoh %>% column_to_rownames("SampleID") %>%
                                        dplyr::mutate(row_sum = rowSums(.)) %>%
                                        dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
                                        dplyr::select(-row_sum),
                                   ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_etoh_ra_scores <- data.frame(PCA_sampling_etoh_ra$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

PCA_sampling_meoh_ra <- mixOmics::pca(data_sampling_meoh %>% column_to_rownames("SampleID") %>%
                                        dplyr::mutate(row_sum = rowSums(.)) %>%
                                        dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
                                        dplyr::select(-row_sum),
                                   ncomp = 2, center = TRUE, scale = TRUE)
PCA_sampling_meoh_ra_scores <- data.frame(PCA_sampling_meoh_ra$variates$X) %>% 
  rownames_to_column("SampleID") %>% left_join(metadata_final)

i <- "Donor"
j <- "Weight"

PCA_sampling_etoh_ra_plot <- PCA_sampling_etoh_ra_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#DD3C51", "#313657", "#1F6683"),
            title = paste("PCA - Human Fecal Metabolome EtOH"),
            xlab = paste("PC1 (", round(PCA_sampling_etoh_ra$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_etoh_ra$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_etoh_ra_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

PCA_sampling_meoh_ra_plot <- PCA_sampling_meoh_ra_scores %>%
  ggscatter(x = "PC1", y = "PC2", color = i, shape = j, alpha = 0.6, palette = c("#214d65", "#287DAB", "#E5BF86"),
            title = paste("PCA - Human Fecal Metabolome MeOH"),
            xlab = paste("PC1 (", round(PCA_sampling_meoh_ra$prop_expl_var$X[1]*100, digits = 1),"%)", sep = ""), 
            ylab = paste("PC2 (", round(PCA_sampling_meoh_ra$prop_expl_var$X[2]*100, digits = 1),"%)", sep = ""),
            ggtheme = theme_classic()) +
  geom_point(data = PCA_sampling_meoh_ra_scores %>% group_by((!!sym(i))) %>% 
               summarise_at(vars(matches("PC")), mean), aes(PC1, PC2, color = (!!sym(i))), size = 4, shape = 8) +
  theme(plot.title = element_text(size = 10), axis.title = element_text(size = 8),
        axis.text = element_text(size = 4)) + coord_fixed()

pca_combined_ra <- ggarrange(PCA_sampling_etoh_ra_plot, PCA_sampling_meoh_ra_plot, nrow = 1)

#ggsave(plot = pca_combined, filename = "SupplementaryFigure4B.svg", device = "svg", dpi = "retina", width = 8)

# PERMANOVA
dist_metabolites_etoh_ra <- vegdist(data_sampling_etoh %>% column_to_rownames("SampleID") %>%
                                   dplyr::mutate(row_sum = rowSums(.)) %>%
                                   dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
                                   dplyr::select(-row_sum), method = "euclidean")
disper_donor_etoh_ra <- betadisper(dist_metabolites_etoh_ra, PCA_sampling_etoh_ra_scores$Donor)
anova(disper_donor_etoh_ra)
permanova_etoh_ra <- adonis2(dist_metabolites_etoh_ra ~ Donor * Weight, 
                          PCA_sampling_etoh_ra_scores, na.action = na.omit)

dist_metabolites_meoh_ra <- vegdist(data_sampling_meoh %>% column_to_rownames("SampleID") %>%
                                   dplyr::mutate(row_sum = rowSums(.)) %>%
                                   dplyr::mutate(across(everything(), ~ . / row_sum)) %>%
                                   dplyr::select(-row_sum), method = "euclidean")
disper_donor_meoh_ra <- betadisper(dist_metabolites_meoh_ra, PCA_sampling_meoh_ra_scores$Donor)
anova(disper_donor_meoh_ra)
permanova_meoh_ra <- adonis2(dist_metabolites_meoh_ra ~ Donor * Weight, 
                          PCA_sampling_meoh_ra_scores, na.action = na.omit)


# Generate Upset plot to check metabolic features recovered by weight
data_sampling_presence <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_sampling$SampleID) %>% 
  column_to_rownames("SampleID")
data_sampling_presence[data_sampling_presence > 0] <- 1
data_sampling_presence_info <- data_sampling_presence %>% rownames_to_column("SampleID") %>% 
  left_join(metadata_final %>% dplyr::select(SampleID, Weight)) %>% column_to_rownames("SampleID")

list_weight <- list(`10mg` = data_sampling_presence_info %>% dplyr::filter(Weight == "10mg") %>% 
                       dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% colnames(),
                     `20mg` = data_sampling_presence_info %>% dplyr::filter(Weight == "20mg") %>% 
                       dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% colnames(),
                     `30mg` = data_sampling_presence_info %>% dplyr::filter(Weight == "30mg") %>% 
                       dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% colnames())

upset_weight <- upset(fromList(list_weight), order.by = "freq") # save pdf for Figure 4B


# Identified features detected in all samples
features_10mg <- data_sampling_presence_info %>% dplyr::filter(Weight == "10mg") %>%
  dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

features_20mg <- data_sampling_presence_info %>% dplyr::filter(Weight == "20mg") %>%
  dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

features_30mg <- data_sampling_presence_info %>% dplyr::filter(Weight == "30mg") %>%
  dplyr::select(-Weight) %>% select(where(~sum(.) >= 1)) %>% t() %>% as.data.frame() %>% 
  rownames_to_column("Feature") %>% dplyr::select(1)

intersect_sampling <- features_10mg %>% inner_join(features_20mg) %>% 
  inner_join(features_30mg) %>% left_join(info_feature_complete)

all_extractions <- features_10mg %>% full_join(features_20mg) %>% 
  full_join(features_30mg) %>% left_join(info_feature_complete) %>% 
  distinct(Feature, .keep_all = TRUE) %>%
  dplyr::filter(!(is.na(Compound_Name)))

# Pie chart for NPC Pathway
intersect_pred <- intersect_sampling %>% group_by(`NPC#pathway`) %>% 
  summarise(Count = n()) %>% arrange(desc(Count)) %>%
  dplyr::mutate(Ratio = Count/sum(Count))

colnames(intersect_pred)[1] <- "NPC Pathway"

pie_intersect <- intersect_pred %>% 
  ggpie(x = "Count", fill = "NPC Pathway", legend = "right") + scale_fill_viridis_d()

#ggsave(plot = pie_intersect, filename = "Figure4B_pie.svg", device = "svg", dpi = "retina")


# Check correlation to gDNA
dna <- read_csv("dna_yield.csv") %>%
  dplyr::mutate(Info = gsub("\\.", "_", Sample)) %>% 
  left_join(metadata_final) %>% dplyr::filter(Time == "1wk") %>%
  dplyr::mutate(Weight = as.numeric(gsub("mg", "", Weight)))

dna_plot <- dna %>% ggscatter(x = "Weight", y = "Concentration", add = "reg.line", 
                              color = "Donor", legend = "right", xlab = "Weight (mg)",
                              ylab = "gDNA (ng/uL)", palette = c("#DD3C51", "#313657", "#1F6683"),
                              title = "Reation between sample weight and extracted gDNA") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = dna_plot, filename = "SupplementaryFigure5A.svg", device = "svg", dpi = "retina", width = 3, height = 3)


# Linear mixed-effects model to take into account repeated measures on subjects
model_concentration <- lmerTest::lmer(Concentration ~ Weight + (1 | Donor), data = dna)

# Summarize the model to see the results
summary(model_concentration)

# Check correlation to extracted peak areas
data_weight <- data_sample %>% 
  dplyr::filter(SampleID %in% sample_sampling$SampleID) %>% 
  dplyr::filter(!(str_detect(pattern = "A|B|C", SampleID))) %>%
  column_to_rownames("SampleID")

data_weight_info <- data_weight %>% rowSums() %>% as.data.frame() %>%
  rownames_to_column("SampleID") %>% left_join(dna)

colnames(data_weight_info)[2] <- "TIC"

data_weight_info <- data_weight_info %>% dplyr::mutate(LogTIC = log10(TIC))

tic_plot <- data_weight_info %>%
  ggscatter(x = "Weight", y = "LogTIC", add = "reg.line", 
            color = "Donor", legend = "right", xlab = "Weight (mg)",
            ylab = "log10(Total Extracted Area)", palette = c("#DD3C51", "#313657", "#1F6683"),
            title = "Reation between sample weight and extracted peak areas") + 
  theme(legend.position = "none",
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = tic_plot, filename = "SupplementaryFigure5B.svg", device = "svg", dpi = "retina", width = 3, height = 3)


# Linear mixed-effects model 
model_tic <- lmerTest::lmer(LogTIC ~ Weight + (1 | Donor), data = data_weight_info)

# Summarize the model to see the results
summary(model_tic)

# Check correlation between gDNA and extracted peak areas
tic_gDNA <- data_weight_info %>%
  ggscatter(x = "Concentration", y = "LogTIC", add = "reg.line", 
            color = "Donor", legend = "right", xlab = "gDNA (ng/uL)",
            ylab = "Log10(Total Extracted Area)",
            title = "Relation between extracted gDNA and peak areas",
            palette = c("#DD3C51", "#313657", "#1F6683")) + 
  theme(legend.position = "none",
        plot.title = element_text(size = 10),
        axis.title = element_text(size = 8),
        axis.text = element_text(size = 6))

#ggsave(plot = tic_gDNA, filename = "Figure4C.svg", device = "svg", dpi = "retina", width = 3, height = 3)

# Linear mixed-effects model 
model_tic_conc <- lmerTest::lmer(LogTIC ~ Concentration + (1 | Donor), data = data_weight_info)

# Summarize the model to see the results
summary(model_tic_conc)

