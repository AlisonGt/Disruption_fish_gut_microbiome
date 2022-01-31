
#packageVersion("phyloseq")

library(seqinr)
library(vegan)
library(VennDiagram)
library(phyloseq)
library(forcats)
library(ggplot2)
library(ggpubr)
library(patchwork)
#library(picante)
library(RColorBrewer)
library(mixOmics)
library(RVAideMemoire)
library(tidyverse)

setwd("~/Documents/THESE/3_Projets/1_Medaka_28j/3_Analyses/3_analyses_microbio/")

#*******************
# METABARCODING ----
#*******************
  
#* Functions ####

#relative abundance function
abund.col = function(col){ 
  return((col/sum(col))*100)
}
abund.col.1 = function(col){ 
  return((col/sum(col)))
}
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}
#transform a dataframe to a matrix
tb2matrix <- function(tb){
  M <- as.matrix(tb[, -1])
  row.names(M) <- as_vector(tb[, 1])
  return(M)
}

#* Colors ####

#colors
col.comp <- c(gut = '#DCC190', water = '#8DC8E0', culture = '#C5E08D', biofilm = '#A9A6A3', faeces = '#BB8207')
col.gut.d28 <- c(`d28_0` = '#77B5FE', `d28_Z8` = '#F3D617', `d28_1` = '#BBF0B1',
                 `d28_10` = '#03C119', `d28_100` = '#018010')
col.gut.d28.d33.d39 <- c(`d28_0` = '#77B5FE', `d28_Z8` = '#F3D617', `d28_1` = '#BBF0B1',
                         `d28_10` = '#03C119', `d28_100` = '#018010', `d33` = '#77B5FE', `d39` = '#018010')
col.gut.d28.0.100.d33.d39 <- c(`d28_0` = "#77B5FE", `d28_100` = "#018010",
                               `d33_0_0` = "#77B5FE", `d33_Z8_0` = "#77B5FE", `d33_1_0` = "#77B5FE",
                               `d33_10_0` = "#77B5FE", `d33_100_0` = "#77B5FE",
                               `d39_0_0_100` = "#018010", `d39_Z8_0_100` = "#018010", `d39_1_0_100` = "#018010",
                               `d39_10_0_100` = "#018010", `d39_100_0_100` = "#018010")
col.gut.d33.d39 <- c(`d33` = '#77B5FE', `d39` = '#018010')
col.gut.0.100 = c(d28_0 = "#77B5FE", d28_100 = "#018010",
                  `d33 (0-100)` = "#79b6fe",`d39 (0-100)` = "#028011")
col.tissue <- c(gut = '#DCC190', liver = 'maroon', muscle = 'pink')


#* Tables ####

#metadata
metadata <- read_csv(file = "2_R/data/2_metabarcoding/metadata.csv") %>% 
  filter(publish != "d28_Z8c_g3" & publish != "d28_10b_g4" & publish != "food") #remove samples under rarefaction threshold
##new columns
metadata$pub.treat <- as_factor(paste(metadata$time, metadata$treatment, sep = "_")) #figure
metadata$pub.treat.hist <- as_factor(paste(metadata$time, metadata$historic, sep = "_")) #figure2
#factor
metadata$group <- factor(metadata$group, levels = c("culture","water","biofilm","gut","faeces"))

#ASV table: rarefied at 6978 reads (on QIIME2) (1_qiime2/8_core.metrics.nochl2/rarefied-table/rarefied-table.txt)
#and remove contaminated samples and rename column names with published names
df <- read_tsv(file = "2_R/data/2_metabarcoding/rarefied-table.txt", skip = 1) %>%
  dplyr::rename(OTUID = "#OTU ID") %>% #rename
  dplyr::select(c("OTUID", metadata$id)) %>% #reorder
  rename_at(vars(`0Ag1`:`39Ew`), ~ str_replace_all(., setNames(metadata$publish, metadata$id))) 

#taxonomy
taxo <- read_tsv(file = "2_R/data/2_metabarcoding/taxonomy.tsv") %>%
  dplyr::rename(FeatureID = "Feature ID")

#sequences : create dataframe from a fasta file (install.packages("seqinr"))
sequences <- read.fasta(file = "./2_R/data/2_metabarcoding/sequences.fasta", seqtype = "AA", as.string = T, set.attributes = F) 
sequences <- data.frame(sequences, check.names = F)
sequences <- t(sequences) #transpose
sequences <- data.frame(OTUID = rownames(sequences), sequences = sequences[,1]) #dataframe with 2 columns (OTUID as header)


#* Merging tables ####

#table + taxonomy + sequences + new asv names
rar_df <- merge(df, taxo, by.x = "OTUID", by.y = "FeatureID")
rar_df <- merge(rar_df, sequences, by.x = "OTUID", by.y = "OTUID")
asv_newnames <- paste0("asv", seq(1, nrow(rar_df)))
rar_df <- rar_df %>%
  as_tibble(.) %>% 
  mutate(ASV = asv_newnames) %>% #create new column
  dplyr::select(Confidence, OTUID, ASV, Taxon, sequences, everything()) %>% #reorganise sample columns
  filter(rowSums(.[,-1] > 0) > 0) #keep lines with at least 1 value
#list OTU and ASV
asv.otu <- tibble("OTUID" = rar_df$OTUID, "ASV" = rar_df$ASV)

#* Export table for publication ####
  #write_csv(rar_df,"../../4_Paper/soumissions/1_Nature_com/data/metabarcoding_data.csv")
  
  
#* Phyloseq tables ####
  
#features (matrix)
physeq.otu <- read_tsv(file = "./2_R/data/2_metabarcoding/rarefied-table.txt", skip = 1) %>% 
  dplyr::rename(OTUID = "#OTU ID") %>% 
  dplyr::select(c("OTUID", metadata$id)) %>% #reorder
  rename_at(vars(`0Ag1`:`39Ew`), ~ str_replace_all(., setNames(metadata$publish, metadata$id))) #change names for publication
physeq.otu <- tb2matrix(physeq.otu)

#taxonomy (matrix)
taxo <- taxo[-1,-3] %>% 
  dplyr::rename(OTUID = "FeatureID") %>% 
  separate(., Taxon, into = c("Kingdom","Phylum","Class","Order","Family","Genus","Species"), sep = "; ") %>% 
  mutate_if(is.character, str_replace_all, pattern = "[a-z]__", replacement = "") #D_[0-9]__
physeq.taxo <- tb2matrix(taxo)

#tree
physeq.tree <- read_tree("./2_R/data/2_metabarcoding/tree.nwk") #rooted tree
physeq.tree <- ape::multi2di(physeq.tree)

#metadata (data.frame)
physeq.meta <- tb2matrix(metadata) %>% as.data.frame(.)

#phyloseq preparation
physeq.otu <- otu_table(physeq.otu, taxa_are_rows = TRUE) #physeq format
physeq.taxo <- tax_table(physeq.taxo) #physeq format
physeq.meta <- sample_data(physeq.meta) #physeq format
physeq.tree <- phyloseq::phy_tree(physeq.tree)

#phyloseq object
physeq <- merge_phyloseq(physeq.otu, physeq.tree, physeq.taxo, physeq.meta) #merge to delete asv from taxonomy according to asv table "features"

#phyloseq object into relative abundances then a large dataframe for ggplot
mdf <- transform_sample_counts(physeq, function(x) x/sum(x)) %>% 
  psmelt(.) #long step


#********************
# 1 - COMPARTMENTS ####
#**********************

#* Tables ####

#ASV table in abundance 
comp.abund <- rar_df %>% 
  dplyr::select(ASV, metadata$publish) %>% 
  mutate_if(is.double, abund.col) #check colSums(comp.abund[,-1])

#* Index tables ####

#asv table : filter but don't remove ASV with zeros everywhere because it changes index (samples deleted)
rar_df_alpha <- rar_df[,-c(1:5)] %>% t(.)
#richness
richness <- vegan::specnumber(rar_df_alpha) #hist(richness)
#diversity including abundance (shannon)
shannon <- vegan::diversity(rar_df_alpha, index = "shannon") #hist(shannon)
#equitability
evenness <- shannon / log(richness) #hist(evenness)

#* Export table for publication ####

#create index tibble
Sample <- names(richness)
indices <- tibble(Sample, `Species richness`=richness, Shannon=shannon, Evenness=evenness)
#import the table with number of reads
read.sum <- read_csv(file = "./2_R/data/2_metabarcoding/read_sum.csv") %>% 
  #dplyr::select(-`d28_Z8c_g3`,-`d28_10b_g4`,-c(`food`:`Total_features`)) %>%
  dplyr::select(-c(`Total_frequency`:`Total_features`)) %>%
  rename_at(vars(`d0_0_g1`:`d39_100_0_100_w`), 
            ~ str_replace_all(., setNames(metadata$publish, metadata$id))) %>% 
  dplyr::select(-1) %>% t(.)
names <- rownames(read.sum)
#import previous metadata
tmp <- read_csv(file = "./2_R/data/1_cleaning/metadata.csv") %>% 
  dplyr::rename(Sample = "publish") %>%
  dplyr::select(Sample, group, time)
#join both tables
read.sum <- read.sum %>% 
  as_tibble(.) %>% 
  mutate(Sample = as.character(names)) %>% 
  dplyr::relocate(Sample, .before = 1) %>% 
  dplyr::rename(`Raw reads` = "V1", `Quality filtered reads` = "V2") %>%
  left_join(., indices, by = "Sample") %>% 
  inner_join(., tmp, by = "Sample") %>% 
  dplyr::select(Sample,Compartment=group,`Sampling time (days⁻¹)`=time,everything())
#>>>> Sup-table3 ####
# write_csv(read.sum, "../../4_Paper/data/Sup-table3.csv", col_names = T)

#* Index values ####

#table
tmp <- as.data.frame(cbind(richness, shannon, evenness))
tmp <- tmp %>%
  tibble(sample = rownames(tmp)) %>% 
  left_join(., metadata, by = c("sample"="publish")) 

#richness
tmp %>%
  dplyr::select(sample, richness, group) %>%
  group_by(group) %>%
  summarise(mean = mean(richness), sd = sd(richness))

#evenness
tmp %>%
  dplyr::select(sample, evenness, group) %>%
  group_by(group) %>%
  summarise(mean = mean(evenness), sd = sd(evenness))

#* Index plots ####

#richness
richness_plot <- data.frame(richness = richness, group = metadata$group)
p1 <- ggplot(data = richness_plot, aes(x = group, y = richness, fill = group)) +
  theme_pubr() +
  geom_violin(alpha = 0.5, trim = F) +
  geom_jitter(position = position_jitter(0.05), alpha = 0.4, color = "dimgrey", size = 1.5) +
  labs(y = "Species richness") +
  scale_fill_manual(values = col.comp) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text = element_text(size = 13), axis.title.y = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA)) +
  stat_summary(fun = mean, geom = "point", show.legend = F, size = 2.5, color = "maroon")

#shannon
shannon_plot <- data.frame(shannon = shannon, group = metadata$group)
p2 <- ggplot(data = shannon_plot, aes(x = group, y = shannon, fill = group)) +
  theme_pubr() +
  geom_violin(alpha = 0.5, trim = F) +
  geom_jitter(position = position_jitter(0.05), alpha = 0.4, color = "dimgrey", size = 1.5) +
  labs(y = "Shannon") +
  scale_fill_manual(values = col.comp) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text = element_text(size = 13), axis.title.y = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA)) +
  stat_summary(fun = mean, geom = "point", show.legend = F, size = 2.5, color = "maroon")

#evenness
evenness_plot <- data.frame(evenness = evenness, group = metadata$group)
p3 <- ggplot(data = evenness_plot, aes(x = group, y = evenness, fill = group)) +
  theme_pubr() +
  geom_violin(alpha = 0.5, trim = F) +
  geom_jitter(position = position_jitter(0.05), alpha = 0.4, color = "dimgrey", size = 1.5) +
  labs(y = "Evenness") +
  scale_fill_manual(values = col.comp) +
  theme(legend.position = "none", axis.title.x = element_blank(),
        axis.text = element_text(size = 13), axis.title.y = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA)) +
  stat_summary(fun = mean, geom = "point", show.legend = F, size = 2.5, color = "maroon")


#* PCoA ####

#relevel
sample_data(physeq)$group <- fct_relevel(sample_data(physeq)$group, 
                                              "culture","water","biofilm","gut","faeces")
#distance matrix
pcoa.w <- ordinate(physeq = physeq, method = "PCoA", distance = "wunifrac") # weighted unif

#PCoA weighted unifrac
p4 <- plot_ordination(physeq = physeq, ordination = pcoa.w, color = "group") + 
  geom_point(size = 4) + 
  theme_pubr() +
  theme(legend.position = "none", axis.text = element_text(size = 13), 
        axis.title = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA)) +
  scale_color_manual(values = alpha(col.comp, 0.9)) +
  labs(colour = "Group")


#>>>> Sup-figure1 ####

#design
layout <- '
AABB
CCDD'

#plot
p <- (p1 | p2) / (p3 | p4) +
    patchwork::plot_annotation(tag_levels = 'a') & 
    theme(plot.tag = element_text(size = 16, hjust = 0, face = "bold"),
          plot.tag.position = "topleft")
#save
#ggsave("Sup-figure1.pdf", p, "pdf", "../../4_Paper/figures", width = 10, height = 8)


#* Signif. ASVs ####

#LEfSe ASVs = dominant asv retained from the 2 LEfSe analyses gut d28 (5 treatments or 2 treatments) on Galaxy (LDA thresholds > 3.5)
LEfSe.asv.dom <- c("asv1349","asv1564","asv1644","asv1662","asv2042","asv2363")
#Significant ASVs = dominant asv retained from LEfSe analyses gut d28 (5 treatments or 2 treatments) on Galaxy (LDA thresholds > 3.5) + one dominant (Cetboacterium)
Signif.asv <- c("asv1349","asv1564","asv1620","asv1644","asv1662","asv2042","asv2363")
#Significant ASVs + ASVs retained network d28
Signif.network.asv <- c("asv475","asv1349","asv1543","asv1564","asv1620","asv1644",
                        "asv1662","asv1995","asv2350","asv2363") #"asv2042",
#Significant ASVs + ASVs retained network d28 0-100
Signif.network.asv.2 <- c("asv1349","asv1536","asv1941","asv1662","asv1543","asv2042")
#Significant ASVs + ASVs retained network d33-39
Signif.network.asv.3 <- c("asv1406","asv2161","asv1030","asv2363")
#Dominant ASVs previously described
Dom.asv <- c("asv1620","asv1349","asv1564","asv1644","asv1662","asv2363")
#Other discriminant ASVs previously described
Dis.asv <- c("asv475","asv1543","asv1995","asv2350")
Dis.asv2 <- c("asv475","asv1030","asv1406","asv1543","asv1995","asv2161","asv2350")
#All ASVs
all.asv <- c("asv475","asv1030","asv1349","asv1406","asv1543","asv1564","asv1620","asv1644",
             "asv1662","asv1995","asv2161","asv2350","asv2363")

#* Dominant ASV ----

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% Dom.asv) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  mutate(na = 0) %>% 
  dplyr::rename(Compartment = "group") %>% 
  group_by(Compartment)
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1620","asv1349","asv1564","asv1644","asv1662","asv2363"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1620**<br>(g. *Cetobacterium*)","**ASV1349**<br>(g. *Flavobacterium*)",
          "**ASV1564**<br>(g. *Aeromonas*)","**ASV1644**<br>(d. Bacteria)",
          "**ASV1662**<br>(g. *ZOR0006*)","**ASV2363**<br>(g. *Reyranella*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels: rar_df %>% dplyr::filter(ASV %in% Dom.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p2 <- ggplot(data = tmp, aes(x = Compartment, y = abund, fill = Compartment)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_wrap(~ ASV, scales = "free_y", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.comp) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 9),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.x = element_markdown(size = 9)) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)
p2
#>>>> Sup-figure 6 ----
#ggsave("Sup-figure6", p2, "pdf","../../4_Paper/figures", width = 7.5, height = 5.5)

##design
#layout <- '
#AAB'
#
##plot
#p <- p1+p2+
#  plot_layout(design = layout) & 
#  plot_annotation(tag_levels = 'a') & 
#  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
#        plot.tag.position = "topleft")
##ggsave("Figure6a-b", p, "pdf","../../4_Paper/figures", width = 7, height = 8.5)



#* Discriminant ASV ----

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% Dis.asv) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  mutate(na = 0) %>% 
  dplyr::rename(Compartment = "group") %>% 
  group_by(Compartment)
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv475","asv1543","asv1995","asv2350"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV475**<br>(g. *Shewanella*)","**ASV1543**<br>(g. *Epulopiscium*)",
          "**ASV1995**<br>(g. *Shewanella*)","**ASV2350**<br>(f. Barnesiellaceae)")
names(labs) <- levels(tmp$ASV)
#check taxa levels: rar_df %>% dplyr::filter(ASV %in% Dis.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p4 <- ggplot(data = tmp, aes(x = Compartment, y = abund, fill = Compartment)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ na, scales = "free_y", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.comp) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)  #+
p4
#>>>> Sup-figure 7a-b ----

#design
layout <- '
AAB'

#plot
p <- p3+p4+
  plot_layout(design = layout) & 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#  ggsave("Sup-figure7a-b", p, "pdf","../../4_Paper/figures", width = 7, height = 6.5)




#* Discriminant ASV 2 ----

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% Dis.asv2) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  mutate(na = 0) %>% 
  dplyr::rename(Compartment = "group") %>% 
  group_by(Compartment)
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv475","asv1030","asv1406","asv1543","asv1995",
                                      "asv2161","asv2350"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV475**<br>(g. *Shewanella*)","**ASV1030**<br>(g. *Reyranella*)",
          "**ASV1406**<br>(g. *Vibrio*)","**ASV1543**<br>(g. *Epulopiscium*)",
          "**ASV1995**<br>(g. *Shewanella*)","**ASV2161**<br>(g. *Vibrio*)",
          "**ASV2350**<br>(f. Barnesiellaceae)")
names(labs) <- levels(tmp$ASV)
#check taxa levels: rar_df %>% dplyr::filter(ASV %in% Dis.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p4 <- ggplot(data = tmp, aes(x = Compartment, y = abund, fill = Compartment)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ na, scales = "free_y", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.comp) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 8.5),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.2)  #+
p4
#>>>> Sup-figure 7a-b ----

#design
layout <- '
AAB'

#plot
p <- p3+p4+
  plot_layout(design = layout) & 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 14, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#ggsave("Sup-figure7a-b_new", p, "pdf","../../4_Paper/figures", width = 7, height = 9.8)


  
#* Other ASV plots ----

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% Signif.network.asv.2) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  mutate(na = 0) %>% 
  dplyr::rename(Compartment = "group") %>% 
  group_by(Compartment)
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1349","asv1536","asv1941","asv1662","asv1543","asv2042"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1349**<br>(g. *Flavobacterium*)","**ASV1536**<br>(g. *Xanthobacter*)",
          "**ASV19415**<br>(o. *Rhizobiales*)","**ASV1662**<br>(g. *ZOR0006*)",
          "**ASV1543**<br>(g. *Epulopiscium*)","**ASV2042**<br>(g. *Romboutsia*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels: rar_df %>% dplyr::filter(ASV %in% Dis.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p4 <- ggplot(data = tmp, aes(x = Compartment, y = abund, fill = Compartment)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ na, scales = "free_y", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.comp) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)  #+
p4
#>>>> Sup-figure X ----

#design
layout <- '
AAB'

#plot
p <- p3+p4+
  plot_layout(design = layout) & 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#  ggsave("Sup-figureX", p, "pdf","../../4_Paper/figures", width = 7, height = 8.5)


#new table
Signif.network.asv.3

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% Signif.network.asv.3) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  mutate(na = 0) %>% 
  dplyr::rename(Compartment = "group") %>% 
  group_by(Compartment)
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1406","asv2161","asv1030","asv2363"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1406**<br>(g. *Vibrio*)","**ASV2161**<br>(g. *Vibrio*)",
          "**ASV1030**<br>(g. *Reyranella*)","**ASV2363**<br>(g. *Reyranella*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels: rar_df %>% dplyr::filter(ASV %in% Dis.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p4 <- ggplot(data = tmp, aes(x = Compartment, y = abund, fill = Compartment)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ na, scales = "free_y", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.comp) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 10),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)  #+
p4
#>>>> Sup-figure XX ----

#design
layout <- '
AAB'

#plot
p <- p3+p4+
  plot_layout(design = layout) & 
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#  ggsave("Sup-figureXX", p, "pdf","../../4_Paper/figures", width = 7, height = 7)





#* Proportions ----

#table
tmp <- comp.abund %>% #check colSums(comp.abund[,-1])
  filter(ASV %in% all.asv) %>% 
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  dplyr::rename(Compartment = "group") %>%
  group_by(ASV, Compartment)

#food compartment
food <- read_delim(file = "2_R/data/1_cleaning/asv-table.txt", delim = "\t", col_names = T,
                   skip = 1) %>%
  rename(OTUID = "#OTU ID") %>% 
  dplyr::select(OTUID, food) %>%
  filter(OTUID %in% asv.otu$OTUID) %>% #keep ASV after cleaning and filtering
  mutate_if(is.double, abund.col) %>%  #check colSums(food[,-1])
  left_join(., asv.otu, by = "OTUID") %>% 
  dplyr::select(ASV, food) %>%
  filter(rowSums(.[, -1] > 0) > 0)
#* FOOD (absent in all but ASV1620)
food[food$ASV %in% all.asv,]

#culture compartment
tmp %>% filter(Compartment == "culture") %>% 
  summarise(m = mean(abund), sd = sd(abund))
#water compartment
tmp %>% filter(Compartment == "water") %>% 
  summarise(m = mean(abund), sd = sd(abund))
#biofilm compartment
tmp %>% filter(Compartment == "biofilm") %>% 
  summarise(m = mean(abund), sd = sd(abund))
#gut compartment
tmp %>% filter(Compartment == "gut") %>% 
  summarise(m = mean(abund), sd = sd(abund))
#faeces compartment
tmp %>% filter(Compartment == "faeces") %>% 
  summarise(m = mean(abund), sd = sd(abund))

#some asv
tmp %>% filter(ASV == "asv2363") %>% 
  summarise(m = mean(abund), sd = sd(abund))
tmp %>% filter(ASV == "asv2363") %>% 
  dplyr::filter(Compartment == "water")
tmp %>% filter(ASV == "asv1662") %>% 
  summarise(m = mean(abund), sd = sd(abund))


#*****************
# 2 - GUT d28 ####
#*****************

#* Tables ####
  
#ASV table
gut.d28 <- rar_df %>% 
  dplyr::select(ASV, metadata$publish[metadata$group == "gut" & metadata$time == "d28"]) 
#ASV table in abundance
gut.d28.abund <- merge(taxo, rar_df, by = "OTUID")
gut.d28.abund <- gut.d28.abund %>%
  as_tibble(.) %>% 
  dplyr::select(ASV, OTUID, Kingdom, Phylum, Class, Order, Family, Genus, Species,
                sequences, metadata$publish[metadata$group == "gut" & metadata$time == "d28"]) %>% 
  mutate_if(is.double, abund.col) #check colSums(gut.d28.abund[,-(1:10)])
#abundant ASV table
gut.d28.abundant <- gut.d28.abund %>%
  filter(rowSums(.[, -(1:10)] >= 1) > 0) #keep lines with at least one TRUE (>= 1%)
nrow(gut.d28.abundant) #70
#dominant ASV table
gut.d28.dominant <- gut.d28.abund %>%
  filter(rowSums(.[, -(1:10)] >= 10) > 0) #keep lines with at least one TRUE (>= 10%)
nrow(gut.d28.dominant) #8

#metadata
metadata.gut.d28 <- metadata %>% 
  dplyr::filter(time == "d28" & group == "gut")
metadata.gut.d28$pub.treat <- as_factor(paste(metadata.gut.d28$time, metadata.gut.d28$treatment, sep = "_"))
##relabel
#metadata.gut.d28$pub.treat <- factor(metadata.gut.d28$pub.treat,
#                                     levels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100"))
#metadata.gut.d28 <- metadata.gut.d28 %>% arrange(pub.treat)

#phyloseq table
physeq.gut.d28 <- subset_samples(physeq, group == "gut" & time == "d28")
#relevel
sample_data(physeq.gut.d28)$pub.treat <- fct_relevel(sample_data(physeq.gut.d28)$pub.treat, 
                                                     "d28_0","d28_Z8","d28_1","d28_10","d28_100")
#in abundance
physeq.gut.d28.abund = transform_sample_counts(physeq.gut.d28, function(x){x / sum(x)}) %>% 
  psmelt(.) #long step


#* Index tables ####

#asv table : filter but don't remove ASV with zeros everywhere because it changes index (samples deleted)
rar_df_alpha <- gut.d28[,-1] %>% t(.)
#richness
richness <- specnumber(rar_df_alpha) #min(richness) = 42 #max(richness) = 219
#diversity including abundance (shannon)
shannon <- diversity(rar_df_alpha, index = "shannon") #hist(shannon)
#equitability
evenness <- shannon / log(richness) #hist(evenness)


#* Index values ####

#table
tmp <- as.data.frame(cbind(richness, shannon, evenness))
tmp <- tmp %>%
  tibble(sample = rownames(tmp)) %>% 
  left_join(., metadata, by = c("sample"="publish"))
#richness
tmp %>%
  dplyr::select(sample, richness, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(richness), sd = sd(richness))

#shannon
tmp %>%
  dplyr::select(sample, shannon, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(shannon), sd = sd(shannon))

#evenness
tmp %>%
  dplyr::select(sample, evenness, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(evenness), sd = sd(evenness))

 
#* Plot Bars ####

#Phylum

#subset data with metadata
tab <- physeq.gut.d28.abund %>% 
  mutate(Phylum = fct_explicit_na(Phylum, na_level = "Unassigned")) %>% #missing phylum in one variable (mdf1$Phylum[mdf1$Phylum == "NA"])
  group_by(Phylum, Sample) %>% #group by phyla then sample
  summarise(Abundance = sum(Abundance*100), Treatment = unique(pub.treat)) %>% #calculate sum of reads on each phylum and put in a new column + keep condition column
  ungroup()
#phylum <1%
tab$Phylum <- as.character(tab$Phylum) #convert to character
tab$Phylum[tab$Abundance < 1] <- "<1%" #change Phylum names if rel. abund is below 0.01%
#rearrange phylum levels according to mean of their relative abundances
tab$Phylum <- fct_rev(fct_reorder(tab$Phylum, tab$Abundance, mean)) %>%
  fct_drop() #automatically relevel according to current level
#color plot with associated phylum names
fig.col <- c(brewer.pal(11,"Set3"),"wheat3","aliceblue", brewer.pal(10, "Paired"), brewer.pal(3, "Set2"))
names(fig.col) <- c("Fusobacteriota","Firmicutes","Proteobacteria","Bacteroidota","Unassigned","Planctomycetota",
                    "Actinobacteriota","Cyanobacteria","Verrucomicrobiota","Myxococcota",
                    "Chloroflexi","Dependentiae","<1%","Chlamydiae","Dadabacteria",
                    "Acidobacteria","Armatimonadetes","Spirochaetes","Deinococcus-Thermus",
                    "Rokubacteria","Patescibacteria","Omnitrophicaeota","Nitrospirae","Gemmatimonadetes","Fibrobacteres",
                    "Elusimicrobia")
#plot
p1 <- ggplot(tab, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = NA, width = 0.8) +
  scale_y_continuous(name = "Relative abundance (%)", expand = c(0,0)) +
  scale_fill_manual(values = fig.col, breaks = names(fig.col)) +
  facet_wrap(~Treatment, nrow = 1, scales = "free_x") +
  theme_pubr() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), axis.title.y = element_text(size = 15),
        axis.text.y = element_text(size = 13),
        legend.justification = "left",
        legend.position = "right", legend.text = element_text(size = 13),
        legend.title = element_text(size = 14), strip.text = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank())
p1


#Firmicutes, families

#subset data with metadata
tab <- physeq.gut.d28.abund %>% 
  filter(Phylum == "Firmicutes") %>% 
  mutate(Family = fct_explicit_na(Family, na_level = "Unassigned")) %>%
  mutate(Abundance = Abundance*100) #Family
#Others Firmicutes
tab$Family <- as.character(tab$Family) 
tab$Family[tab$Abundance < 1] <- "Others Firmicutes" 
#rearrange family levels according to mean of their relative abundances
tab$Family <- fct_rev(fct_reorder(tab$Family, tab$Abundance, mean)) %>% fct_drop()
#color plot
fig.col <- c("#f6f076","lightgoldenrod3","gold","lightgoldenrod4")
names(fig.col) <- levels(tab$Family)
#plot
p2 <- ggplot(tab, aes(x = Sample, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "stack", color = NA, width = 0.8) +
  scale_y_continuous(name = "Relative abundance (%)", limits = c(0,100), expand = c(0,0)) +
  scale_fill_manual(values = fig.col, breaks = names(fig.col)) +
  facet_wrap(~pub.treat, nrow = 1, scales = "free_x") +
  theme_pubr() +
  theme(axis.text.x = element_blank(), axis.title.x = element_blank(),
        axis.ticks.x = element_blank(), axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 15), legend.text = element_text(size = 13),
        legend.position = "right", legend.justification = "left",
        legend.title = element_text(size = 14), strip.text = element_text(size = 15),
        panel.border = element_rect(colour = "black", fill = NA),
        strip.background = element_blank())
p2


#* PCoA ####

#distance matrix
pcoa.unw <- ordinate(physeq = physeq.gut.d28, method = "PCoA", distance = "unifrac", weighted = F) # unweighted unifrac
pcoa.w <- ordinate(physeq = physeq.gut.d28, method = "PCoA", distance = "wunifrac") # weighted unifrac

#PCoA weighted unifrac
p3 <- plot_ordination(physeq = physeq.gut.d28, ordination = pcoa.w, color = "pub.treat") + 
  theme_pubr() +
  geom_point(size = 4.5) + 
  scale_color_manual(values = col.gut.d28) +
  labs(title = "PCoA (wU) on the gut microbiome at d28", colour = "Treatment") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        legend.position = "right", legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = rel(1.4), face = "bold.italic", 
                                  colour = "gray20", hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA))
#p3 <- p3 + coord_equal(ratio = 0.8)

#PCoA weighted unifrac
p4 <- plot_ordination(physeq = physeq.gut.d28, ordination = pcoa.unw, color = "pub.treat") + 
  theme_pubr() +
  geom_point(size = 4.5) + 
  scale_color_manual(values = col.gut.d28) +
  labs(title = "PCoA (unwU) on the gut microbiome at d28", colour = "Treatment") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        legend.position = "right", legend.title = element_text(size = 14),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = rel(1.4), face = "bold.italic", 
                                  colour = "gray20", hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA))
#p4 <- p4 + coord_equal(ratio = 0.8)


#>>>> Figure 2a-b-c-d ####

#design
layout <- '
AB
AB
CC
DD'

#plot
p <- p3+p4+p1+p2+
  plot_layout(design = layout, guides = 'collect') +
  plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 20, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
p
#save
#ggsave("Figure2a-b-c-d.pdf", p, "pdf", "../../4_Paper/figures", width = 14, height = 12)



#* Permanova ####

#weighted UniFrac distance
set.seed(1473)
dist <- phyloseq::distance(physeq.gut.d28, method = "wunifrac")

#permanova (comparisons of intra and inter-group variatiosn)
vegan::adonis(dist ~ pub.treat, data = metadata.gut.d28, permutations = 999) #p-value = 0.001

#pairwise comparisons (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, metadata.gut.d28$pub.treat, nperm = 999,
                                    p.method = "bonferroni")

#* Permdisp ####

#dispersion (comparisons of intra-group variations)
disp <- vegan::betadisper(dist, group = metadata.gut.d28$pub.treat)
permutest(disp,  pairwise = T)


#* Proportions Phylums in 0 ####

#subset only d28_0 and calculate mean/sd of each main phylums across d28_0 samples
physeq.gut.d28.abund %>%
  dplyr::filter(pub.treat == "d28_0") %>%
  filter(Phylum %in% c("Fusobacteriota","Firmicutes","Proteobacteria","Bacteroidota")) %>%
  group_by(Phylum, Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  group_by(Phylum) %>% 
  summarise(mean = mean(sum), sd = sd(sum))

#subset only d28_0 and calculate mean across samples of the sum of the main phylums per sample
physeq.gut.d28.abund %>%
  dplyr::filter(pub.treat == "d28_0") %>% 
  dplyr::filter(Phylum %in% c("Fusobacteriota","Firmicutes","Proteobacteria","Bacteroidota")) %>% 
  group_by(Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  summarise(mean = mean(sum), sd = sd(sum))
  

#* Proportions Firmicutes ####

#subset Firmicutes and calculate mean/sd of Firmicutes in each treatment
physeq.gut.d28.abund %>%
  filter(Phylum %in% "Firmicutes") %>%
  group_by(pub.treat,Sample) %>% 
  summarise(sum = sum(Abundance)*100) %>%
  group_by(pub.treat) %>% 
  summarise(mean = mean(sum), sd = sd(sum))


#* LEfSe ####

#~> LEfSe 5 groups ####

# Tables

#metadata
metadata.gut.d28

#table with no level (only last level displayed) + new ASV names
LEfSe.table <- rar_df %>% 
  dplyr::select(Taxon, ASV, metadata.gut.d28$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "") %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everything before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  filter(rowSums(.[,-c(1:2)] > 0) > 0) #check colSums(LEfSe.table[,-c(1:2)])

#(Export table) ####

#add the ASV name after the Taxon id
LEfSe.Galaxy <- LEfSe.table %>% 
  unite(Taxon, c("Taxon","ASV"))
#export metadata first
LEfSe <- as.data.frame(t(metadata.gut.d28[,c(9,1)]))
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t28-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = T)
#asv table then added to the previous file with only metadata
LEfSe <- as.data.frame(LEfSe.Galaxy)
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t28-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = F, append = T)

# Results

#asv retained from LEfSe analysis on Galaxy (LDA threshold > 3.5)
LEfSe.asv <- c("asv1644","asv1564","asv1941","asv2006","asv2145","asv789","asv857",
               "asv2042","asv952","asv1662","asv1995","asv1349","asv1602","asv1116",
               "asv2350","asv1419")

#which LEfSe asvs are dominant ?
merge(taxo, rar_df, by = "OTUID") %>% 
  as_tibble(.) %>%
  dplyr::filter(ASV %in% LEfSe.asv) %>% 
  dplyr::filter(ASV %in% gut.d28.dominant$ASV) %>% .$ASV

# Dot Plot

#table with cleaned taxa names
LEfSe <- rar_df %>% 
  dplyr::select(Taxon, ASV, metadata.gut.d28$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "") 
#asv retained from LEfSe analysis on Galaxy (LDA threshold > 3.5)
LEfSe.asv
#table
LEfSe.table <- LEfSe %>% 
  dplyr::filter(ASV %in% LEfSe.asv) %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  mutate_if(is.character, str_replace_all, pattern = "__", replacement = ". ") %>% #remove __
  filter(rowSums(.[,-c(1:2)] > 0) > 0)
#factor
LEfSe.table$ASV <- factor(LEfSe.table$ASV, levels = rev(c("asv1644","asv1564","asv1941","asv2006","asv2145","asv789","asv857",
                                                          "asv2042","asv952","asv1662","asv1995","asv1349","asv1602","asv1116",
                                                          "asv2350","asv1419")))
LEfSe.table <- LEfSe.table %>% 
  arrange(ASV) %>% 
  dplyr::mutate(taxon = forcats::as_factor(paste0(Taxon," (",ASV,")"))) %>%
  #dplyr::mutate(taxon = interaction(Taxon, ASV, sep = " ("), Taxon = NULL, ASV = NULL) %>% #join asv and taxa names and keep factor (unite doesn't keep it)
  dplyr::select(taxon, everything()) %>% 
  dplyr::select(-Taxon,-ASV) %>% 
  pivot_longer(cols = !taxon, names_to = "sample") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>%
  select(taxon, sample, value, pub.treat) %>% 
  group_by(pub.treat, taxon) %>%
  summarise(min = min(value, na.rm = T),
            mean = mean(value, na.rm = T),
            max = max(value, na.rm = T))

#change names with uncultured and choose the last level
LEfSe %>% dplyr::filter(ASV %in% LEfSe.asv) %>% dplyr::select(ASV, Taxon)
#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("g. *Macellibacteroides* (ASV1419)","f. Barnesiellaceae (ASV2350)","g. *Hyphomicrobium* (ASV1116)",        
          "g. *Hyphomicrobium* (ASV1602)","**g. *Flavobacterium* (ASV1349)**","g. *Shewanella* (ASV1995)",        
          "**g. *ZOR0006* (ASV1662)**","s. *cf. Leptolyngbya* (ASV952)","**g. *Romboutsia* (ASV2042)**",        
          "o. Rhizobiales (ASV857)","g. *Aurantisolimonas* (ASV789)","s. *Vibrio cholerae* (ASV2145)",        
          "g. *Shewanella* (ASV2006)","o. Rhizobiales (ASV1941)","**g. *Aeromonas* (ASV1564)**",        
          "**d. Bacteria (ASV1644)**")
names(labs) <- levels(LEfSe.table$taxon)

##legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
#labs <- c("s. *uncultured bacterium* (ASV1419)","s. *uncultured bacterium* (ASV2350)","g. *Hyphomicrobium* (ASV1116)",        
#          "g. *Hyphomicrobium* (ASV1602)","**g. *Flavobacterium* (ASV1349)**","g. *Shewanella* (ASV1995)",        
#          "**s. *uncultured bacterium* (ASV1662)**","s. *cf. Leptolyngbya* (ASV952)","**g. *Romboutsia* (ASV2042)**",        
#          "o. Rhizobiales (ASV857)","s. *uncultured Bacteroidetes* (ASV789)","s. *Vibrio cholerae* (ASV2145)",        
#          "g. *Shewanella* (ASV2006)","o. Rhizobiales (ASV1941)","**g. *Aeromonas* (ASV1564)**",        
#          "**d. Bacteria (ASV1644)**")
#names(labs) <- levels(LEfSe.table$taxon)

#plot
library(ggtext)
ggplot(data = LEfSe.table, aes(x = mean, y = taxon, color = pub.treat)) +
  theme_pubr() +
  geom_pointrange(aes(xmin=min, xmax=max), size = 0.7,
                  position = position_dodge(width = 0.6)) +
  scale_color_manual(values = col.gut.d28) +
  scale_y_discrete(labels = labs) +
  theme(legend.position = "right",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = element_markdown(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 14)) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray60") +
  labs(x = "Relative abundance (%)", col = "Treatment") -> p
#>>>> Figure 3 ----
#ggsave("Figure3.pdf", p, "pdf","../../4_Paper/figures", height = 9, width = 9)


#~> LEfSe 0-100 ####

# Tables

#metadata
metadata.gut.d28.0.100 <- metadata.gut.d28 %>% 
  dplyr::filter(pub.treat == "d28_0" | pub.treat == "d28_100")

#table with no level (only last level displayed) + new ASV names
LEfSe.table <- rar_df %>% 
  dplyr::select(Taxon, ASV, metadata.gut.d28.0.100$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "") %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  filter(rowSums(.[,-c(1:2)] > 0) > 0) #check colSums(LEfSe_asv[,-c(1:2)])

#(Export table) ####

#add the ASV name after the Taxon id
LEfSe.Galaxy <- LEfSe.table %>% 
  unite(Taxon, c("Taxon","ASV"))
#export metadata first
LEfSe <- as.data.frame(t(metadata.gut.d28.0.100[,c(9,1)]))
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t28-0-100-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = T)
#asv table then added to the previous file with only metadata
LEfSe <- as.data.frame(LEfSe.Galaxy)
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t28-0-100-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = F, append = T)

# Results

#asv retained from LEfSe analysis on Galaxy (LDA threshold > 3.5)
LEfSe.asv <- c("asv2363","asv1564","asv1349","asv1493",
               "asv1941","asv2400","asv1402","asv253","asv857",
               "asv789","asv2190","asv1995","asv2506","asv1536",
               "asv1598","asv1543","asv2042","asv1936","asv1817",
               "asv2577","asv2379","asv1554","asv1046","asv1419",
               "asv1320","asv1662")

#which LEfSe asvs are dominant ?
merge(taxo, rar_df, by = "OTUID") %>%
  as_tibble(.) %>%
  filter(ASV %in% LEfSe.asv) %>% 
  filter(ASV %in% gut.d28.dominant$ASV) %>% .$ASV

#* (Dot plot) ----

#table with cleaned taxa names
LEfSe <- rar_df %>% 
  dplyr::select(Taxon, ASV, metadata.gut.d28.0.100$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "") 
#asv retained from LEfSe analysis on Galaxy (LDA threshold > 3.5)
LEfSe.asv
#table
LEfSe.table <- LEfSe %>% 
  dplyr::filter(ASV %in% LEfSe.asv) %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  mutate_if(is.character, str_replace_all, pattern = "__", replacement = ". ") %>% #remove everthing before the first ; containing NA
  filter(rowSums(.[,-c(1:2)] > 0) > 0)
#factor
LEfSe.table$ASV <- factor(LEfSe.table$ASV, levels = rev(c("asv2363","asv1564","asv1349","asv1493",
                                                          "asv1941","asv2400","asv1402","asv253","asv857",
                                                          "asv789","asv2190","asv1995","asv2506","asv1536",
                                                          "asv1598","asv1543","asv2042","asv1936","asv1817",
                                                          "asv2577","asv2379","asv1554","asv1046","asv1419",
                                                          "asv1320","asv1662")))
LEfSe.table <- LEfSe.table %>% 
  arrange(ASV) %>% 
  dplyr::mutate(taxon = forcats::as_factor(paste0(Taxon," (",ASV,")"))) %>%
  dplyr::select(taxon, everything()) %>%
  dplyr::select(-Taxon,-ASV) %>% 
  pivot_longer(cols = !taxon, names_to = "sample") %>% 
  left_join(., metadata, by = c("sample"="publish")) %>%
  select(taxon, sample, value, pub.treat) %>% 
  group_by(pub.treat, taxon) %>%
  summarise(min = min(value, na.rm = T),
            mean = mean(value, na.rm = T),
            max = max(value, na.rm = T)) 

#change names with uncultured and choose the last level
LEfSe %>% dplyr::filter(ASV %in% LEfSe.asv) %>% dplyr::select(Taxon) -> tmp
LEfSe %>% dplyr::filter(ASV == "asv1543") %>% dplyr::select(ASV, Taxon)
#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**g. *ZOR0006* (ASV1662)**","f. Caulobacteraceae (ASV1320)", "g. *Macellibacteroides* (ASV1419)",  
          "g. *Rhodococcus* (ASV1046)","g. *Brevundimonas* (ASV1554)","g. *Noviherbaspirillum* (ASV2379)", 
          "o. Enterobacterales (ASV2577)","s. *uncultured planctomycete* (ASV1817)","f. Caulobacteraceae (ASV1936)", 
          "g. *Romboutsia* (ASV2042)", "g. *Epulopiscium* (ASV1543)" , "g. *EV818SWSAP88* (ASV1598)",  
          "g. *Xanthobacter* (ASV1536)", "g. *Reyranella* (ASV2506)" , "g. *Shewanella* (ASV1995)",  
          "f. Polyangiaceae (ASV2190)", "s. *uncultured Bacteroidetes* (ASV789)", "o. Rhizobiales (ASV857)",  
          "g. *Gemmobacter* (ASV253)", "p. Bdellovibrionota (ASV1402)", "g. *Schlesneria* (ASV2400)",  
          "o. Rhizobiales (ASV1941)", "g. *Flavobacterium* (ASV1493)", "**g. *Flavobacterium* (ASV1349)**",  
          "**g. *Aeromonas* (ASV1564)**", "**g. *Reyranella* (ASV2363)**")
names(labs) <- levels(LEfSe.table$taxon)

##legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
#labs <- c("**s. *uncultured bacterium* (ASV1662)**","f. Caulobacteraceae (ASV1320)", "s. *uncultured bacterium* (ASV1419)",  
#          "g. *Rhodococcus* (ASV1046)","g. *Brevundimonas* (ASV1554)","g. *Noviherbaspirillum* (ASV2379)", 
#          "o. Enterobacterales (ASV2577)","s. *uncultured planctomycete* (ASV1817)","f. Caulobacteraceae (ASV1936)", 
#          "g. *Romboutsia* (ASV2042)", "s. *uncultured bacterium* (ASV1543)" , "g. *EV818SWSAP88* (ASV1598)",  
#          "g. *Xanthobacter* (ASV1536)", "g. *Reyranella* (ASV2506)" , "g. *Shewanella* (ASV1995)",  
#          "g. *uncultured* (ASV2190)", "s. *uncultured Bacteroidetes* (ASV789)", "o. Rhizobiales (ASV857)",  
#          "g. *Gemmobacter* (ASV253)", "p. Bdellovibrionota (ASV1402)", "g. *Schlesneria* (ASV2400)",  
#          "o. Rhizobiales (ASV1941)", "g. *Flavobacterium* (ASV1493)", "**g. *Flavobacterium* (ASV1349)**",  
#          "**g. *Aeromonas* (ASV1564)**", "**g. *Reyranella* (ASV2363)**")
#names(labs) <- levels(LEfSe.table$taxon)

#plot
p1 <- ggplot(data = LEfSe.table, aes(x = mean, y = taxon, color = pub.treat)) +
  theme_pubr() +
  geom_pointrange(aes(xmin=min, xmax=max), size = 0.8,
                  position = position_dodge(width = 0.6)) +
  scale_color_manual(values = col.gut.d28) +
  scale_y_discrete(labels = labs) +
  theme(legend.position = "right",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(face = "bold.italic", hjust = 0.5,
                                  color = "gray20"),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = ggtext::element_markdown(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 14)) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray60") +
  labs(title = "Significant ASVs", 
       x = "Relative abundance (%)", col = "Treatment")
#>>>> Sup-figure 3a ----


#* Shotgun Bins ####

# Tables

#bin table from metagenomic analysis
df.shotg <- read_tsv("../6_analyses_metagenome/2_analyses/2_R/data/RESULTS.matrix.bin-sample-RelABDnew.txt") %>% 
  dplyr::rename(Bin = "X1") %>% 
  mutate_if(is.character, str_replace_all, pattern = "Bin_0", replacement = "bin") %>%
  mutate_if(is.character, str_replace_all, pattern = "Bin_", replacement = "bin")

#metadata with only samples from metagenomic
metadata.shotgun <- read_csv(file = "../6_analyses_metagenome/2_analyses/2_R/data/metadata.csv") %>% 
  dplyr::filter(id %in% colnames(df.shotg)[-1]) #keep samples for metagenomic
#factor
metadata.shotgun$pub.treat <- as_factor(paste(metadata.shotgun$time, metadata.shotgun$treatment, sep = "_"))

#change with publication names
df.shotg <- df.shotg %>% 
  rename_at(vars(`28D1g1`:`28E3g2`), ~ str_replace_all(., setNames(metadata.shotgun$publish, metadata.shotgun$id)))
#factor
df.shotg$Bin <- factor(df.shotg$Bin, levels = df.shotg$Bin)

#bin taxo
taxo.shotg <- read_csv("../6_analyses_metagenome/2_analyses/2_R/data/32bins_taxa.csv") %>% 
  dplyr::rename(Bin = "id") %>%
  mutate_if(is.character, str_replace_all, pattern = "Bin_0", replacement = "bin") %>%
  mutate_if(is.character, str_replace_all, pattern = "Bin_", replacement = "bin") %>% 
  mutate(phylum = na_if(phylum, "not classified"),
         class = na_if(class, "not classified"),
         order = na_if(order, "not classified"),
         family = na_if(family, "not classified"),
         genus = na_if(genus, "not classified"),
         species = na_if(species, "not classified")) %>% 
  dplyr::mutate(Taxon = paste0("p__",phylum,";c__",class,";o__",order,";f__",family,";g__",genus,";s__",species)) %>% 
  mutate_if(is.character, str_replace_all, pattern = ";.__NA.*$", replacement = "") %>% #remove everthing before the first ; containing NA
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "__", replacement = ". ") %>% 
  dplyr::select(Bin, Taxon)

#Merge bin table + taxo and add the Bin name after the Taxon id
df.shotg.tax <- df.shotg %>% 
  left_join(., taxo.shotg, by = "Bin")


# Statistical tests

library(rstatix)
library(purrr)
#library(magrittr)

#prepare table
tmp <- df.shotg %>%
  gather(key = key, value = value, 2:ncol(.)) %>% #transposing
  dplyr::rename(sample = "key") %>% 
  left_join(., metadata.shotgun, by = c("sample"="publish")) %>% 
  dplyr::select(sample, Bin, value, pub.treat)

#non parametric test between two groups (wilcoxon test)
test <- tmp %>% 
  split(., .$Bin) %>% 
  purrr::map(~ rstatix::wilcox_test(.x, value ~ pub.treat)) %>% 
  purrr::map_dfr(.x = ., .f = magrittr::extract, .id = "Bin") #extract data
#write_tsv(test, "output/wilcoxon_bins_new.tsv")

#filter Bin with p value below 0/05
bin.0.05 <- test %>% 
  dplyr::filter(p < 0.05)


# Dot plot

#all bin
bin <- bin.0.05$Bin
#bin up in 0
bin.up.0 <- c("bin11","bin22","bin24")
#bin up in 100
bin.up.100 <- setdiff(bin, bin.up.0)

#filter table
tmp1 <- df.shotg.tax %>%
  dplyr::filter(Bin %in% bin.0.05$Bin)
#factor
tmp1$Bin <- factor(tmp1$Bin, levels = rev(c(bin.up.100,bin.up.0)))
#table of significant abundant bins (5%)
tmp1 <- tmp1 %>% 
  arrange(Bin) %>% 
  dplyr::mutate(taxon = forcats::as_factor(paste0(Taxon," (",Bin,")"))) %>%
  #dplyr::mutate(taxon = interaction(Taxon, Bin, sep = "_"), Taxon = NULL, Bin = NULL) %>% #join bin and taxa names and keep factor (unite doesn't keep it)
  dplyr::select(taxon, everything()) %>% 
  dplyr::select(-Taxon, -Bin) %>% 
  pivot_longer(., cols = 2:ncol(.), "sample") %>% #transposing
  left_join(., metadata.shotgun, by = c("sample"="publish")) %>% 
  dplyr::select(sample, taxon, value, pub.treat) %>% 
  group_by(pub.treat, taxon) %>% 
  summarise(min = min(value, na.rm = T),
            mean = mean(value, na.rm = T),
            max = max(value, na.rm = T)) #calculate means

#legend text
labs <- c("p. Firmicutes (bin24)","o. Bacteroidales (bin22)",
          "f. Porphyromonadaceae (bin11)","o. Rhodobacterales (bin26)",
          "o. Rhizobiales (bin25)","o. Rhizobiales (bin23)",
          "g. *Aeromonas* (bin19)","s. *Flavobacterium suncheonense* (bin16)",
          "s. *Gemmobacter aquaticus* (bin12)","o. Rhodospirillales (bin10)",
          "g. *Flavobacterium* (bin9)")
names(labs) <- levels(tmp1$taxon)

#bar plot
p2 <- ggplot(data = tmp1, aes(x = mean, y = taxon, color = pub.treat)) +
  theme_pubr() +
  geom_pointrange(aes(xmin=min, xmax=max), size = 0.8,
                  position = position_dodge(width = 0.5)) +
  scale_color_manual(values = col.gut.d28) +
  scale_y_discrete(labels = labs) +
  theme(legend.position = "right",
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14),
        axis.line = element_blank(),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(face = "bold.italic", hjust = 0.5,
                                  color = "gray20"),
        axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        axis.text.y = ggtext::element_markdown(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 14)) +
  labs(title = "Significant Bins", x = "Relative abundance (%)",
       col = "Treatment")
#>>>> Sup-figure 3b ----
#ggsave("Sup-figure3b.pdf", p2, "pdf", "../../../../4_Paper/figures/", height = 5.5, width = 7.5)

#>>>> Sup-figure 3a-b ----

p <- (p1 | (p2 / guide_area())) +
  patchwork::plot_layout(guides = "collect") +
  patchwork::plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#save
#ggsave("Sup-figure3a-b.pdf", p, "pdf", "../../4_Paper/figures", width = 16, height = 10)



#* Signif. ASVs ####

# Tables

#LEfSe dominant ASVs retained from the 2 LEfSe analyses (5 treatments or 2 treatments) on Galaxy (LDA thresholds > 3.5)
LEfSe.asv.dom <- c("asv1349","asv1564","asv1644","asv1662","asv2042","asv2363")
LEfSe.otu.dom <- asv.otu %>% 
  dplyr::filter(ASV %in% LEfSe.asv.dom) %>% .$OTUID
#Significant dominant ASVs retained from 2 LEfSe analyses (5 treatments or 2 treatments) on Galaxy (LDA thresholds > 3.5) + one dominant (Cetboacterium)
Signif.asv <- c("asv1349","asv1564","asv1620","asv1644","asv1662","asv2042","asv2363")
#ASVs retained network d28
Network.asv <- c("asv475","asv1349","asv1543","asv1564","asv1995","asv2350")
#Significant ASVs + ASVs retained network d28
Signif.network.asv <- c("asv475","asv1349","asv1543","asv1564","asv1620","asv1644",
                        "asv1662","asv1995","asv2042","asv2350","asv2363")
#Dominant ASVs previously described
Dom.asv <- c("asv1349","asv1564","asv1620","asv1644","asv1662","asv2363")

# Proportions among treatments

#subset ASVs and calculate min and max across samples of the sum of the ASVs per sample
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% LEfSe.otu.dom) %>% 
  group_by(Sample) %>%
  summarise(sum = sum(Abundance)*100) %>%
  summarise(min = min(sum), max = max(sum))

#subset each ASV and calculate their mean/sd in each treatment
#Flavobacterium & Aeromonas
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% c("asv1349","asv1564")]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#Reyranella
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv2363"]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#ZOR_0006
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv1662"]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#Bacteria
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv1644"]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#Romboutsia
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv2042"]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100), max = max(Abundance*100))
#Network ASVs: asv475
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% Network.asv[1]]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100), max = max(Abundance*100))
#Network ASVs: asv1995
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% Network.asv[5]]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100), max = max(Abundance*100))
#Network ASVs: asv1543
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% Network.asv[3]]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100), max = max(Abundance*100))
#Network ASVs: asv2350
physeq.gut.d28.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% Network.asv[6]]) %>% 
  group_by(OTU, pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100), max = max(Abundance*100))


# Plot multiple ASVs

#ASV
Network.asv[c(1,3,5,6)]
rar_df %>% 
  dplyr::filter(ASV %in% Network.asv[c(1,3,5,6)]) %>% dplyr::select(ASV, Taxon)

#network ASV across treatments
tmp <- gut.d28.abundant %>% 
  dplyr::select(-(2:10)) %>% 
  dplyr::filter(ASV %in% Network.asv) %>% 
  dplyr::slice(-c(2,4)) %>% 
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28, by = c("sample"="publish"))
#factor
tmp$pub.treat <- factor(tmp$pub.treat,
                        levels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100"))
tmp$ASV <- factor(tmp$ASV,
                  levels = Network.asv[c(1,5,3,6)])
#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV475**<br>(g. *Shewanella*)","**ASV1995**<br>(g. *Shewanella*)",
          "**ASV1543**<br>(g. *Epulopiscium*)","**ASV2350**<br>(f. Barnesiellaceae)")
names(labs) <- levels(tmp$ASV)

#plot
library(ggtext)
p5 <- ggplot(data = tmp, aes(x = pub.treat, y = abund, fill = pub.treat)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  geom_jitter(shape = 16, position = position_jitter(0.2), size = 1.5, color = "grey", alpha = 0.6) +
  facet_wrap(~ ASV, ncol = 2, scales = "free", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28) +
  labs(y = "Relative abundance (%)", fill = "Treatment") +
  theme(legend.position = "top",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.title.x = element_blank(), axis.ticks.x = element_blank(),
        strip.text.x = element_markdown(size = 11),
        axis.text.x = element_blank()) +
  guides(fill = guide_legend(override.aes = list(linetype = 0, shape = NA, alpha = 1)))
#>>>> Figure 4c ----


#*************************
# 3 - GUT d28/d33/d39 ####
#*************************

#* Tables ####

#metadata
metadata.gut.d28.d33.d39 <- metadata %>%
  dplyr::filter((time == "d28" | time == "d33" | time == "d39") & group == "gut") %>% 
  mutate(pub.hist = .$publish) %>% #new column for barplot
  mutate(pub.treat = .$publish) #new column for PCoA
#change sample names
metadata.gut.d28.d33.d39$pub.hist <- c(substr(metadata.gut.d28.d33.d39$pub.hist[1:56], 1, 
                                              nchar(metadata.gut.d28.d33.d39$pub.hist[1:56])-4), #delete the last four characters
                                       substr(metadata.gut.d28.d33.d39$pub.hist[57:91], 1, 
                                              nchar(metadata.gut.d28.d33.d39$pub.hist[57:91])-3)) #delete the last three characters
metadata.gut.d28.d33.d39$pub.treat <- c(substr(metadata.gut.d28.d33.d39$pub.treat[1:56], 1, 
                                               nchar(metadata.gut.d28.d33.d39$pub.treat[1:56])-4), #delete the last four characters
                                        substr(metadata.gut.d28.d33.d39$pub.treat[57:71], 1, 3), #keep first three characters
                                        substr(metadata.gut.d28.d33.d39$pub.treat[72:91], 1, 3)) #keep first three characters
#factor
metadata.gut.d28.d33.d39$pub.hist <- as_factor(metadata.gut.d28.d33.d39$pub.hist) #levels(metadata$pub.hist)
metadata.gut.d28.d33.d39$pub.treat <- as_factor(metadata.gut.d28.d33.d39$pub.treat)

#ASV table
gut.d28.d33.d39 <- rar_df %>% 
  dplyr::select(ASV, metadata.gut.d28.d33.d39$publish)

#phyloseq table
physeq.gut.d28.d33.d39 <- subset_samples(physeq, group == "gut" & 
                                           (time == "d28" | time == "d33" | time == "d39"))
#change names then factor
sample_data(physeq.gut.d28.d33.d39)$pub.treat  <- c(as.character(sample_data(physeq.gut.d28.d33.d39)$pub.treat[1:56]),
                                                    substr(sample_data(physeq.gut.d28.d33.d39)$pub.treat[57:71], 1, 3), #keep first three characters
                                                    substr(sample_data(physeq.gut.d28.d33.d39)$pub.treat[72:91], 1, 3)) #keep first three characters
sample_data(physeq.gut.d28.d33.d39)$pub.hist  <- c(as.character(sample_data(physeq.gut.d28.d33.d39)$pub.treat[1:56]),
                                                   substr(rownames(sample_data(physeq.gut.d28.d33.d39))[57:91], 1,
                                                          nchar(rownames(sample_data(physeq.gut.d28.d33.d39))[57:91])-3)) #delete the last three characters
#factor
sample_data(physeq.gut.d28.d33.d39)$pub.treat <- as_factor(sample_data(physeq.gut.d28.d33.d39)$pub.treat)
sample_data(physeq.gut.d28.d33.d39)$pub.hist <- as_factor(sample_data(physeq.gut.d28.d33.d39)$pub.hist)
#in abundance
physeq.gut.d28.d33.d39.abund = transform_sample_counts(physeq.gut.d28.d33.d39, function(x){x / sum(x)}) %>% 
  psmelt(.) #long step


#* Index tables ####

#asv table : filter but don't remove ASV with zeros everywhere because it changes index (samples deleted)
rar_df_alpha <- gut.d28.d33.d39[,-1] %>% t(.)
#richness
richness <- specnumber(rar_df_alpha) #min(richness) = 42 #max(richness) = 219
#diversity including abundance (shannon)
shannon <- diversity(rar_df_alpha, index = "shannon") #hist(shannon)
#equitability
evenness <- shannon / log(richness) #hist(evenness)


#* Index values ####

#table
tmp <- as.data.frame(cbind(richness, shannon, evenness))
tmp <- tmp %>%
  tibble(sample = rownames(tmp)) %>% 
  left_join(., metadata.gut.d28.d33.d39, by = c("sample"="publish"))

#richness
tmp %>%
  dplyr::select(sample, richness, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(richness), sd = sd(richness))

#shannon
tmp %>%
  dplyr::select(sample, shannon, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(shannon), sd = sd(shannon))

#evenness
tmp %>%
  dplyr::select(sample, evenness, pub.treat) %>%
  group_by(pub.treat) %>%
  summarise(mean = mean(evenness), sd = sd(evenness))


#* Index plots ####

#table
tmp <- tmp %>% 
  pivot_longer(., cols = c(1:3), names_to = "index")
#plot
ggplot(data = tmp, aes(x = pub.treat, y = value, fill = pub.treat)) +
  theme_pubr() +
  geom_boxplot() +
  facet_wrap(~ index, scales = "free_y") +
  scale_fill_manual(values = col.gut.d28.d33.d39)


#* Phylum d33-d39 ####

#subset data with metadata
tab1 <- physeq.gut.d28.d33.d39.abund %>% 
  filter(group == "gut" & time == "d33") %>% 
  mutate(Phylum = fct_explicit_na(Phylum, na_level = "Unassigned")) %>% #missing phylum in one variable (mdf1$Phylum[mdf1$Phylum == "NA"])
  group_by(Phylum, Sample) %>% #group by cond then sample then phylum
  summarise(Abundance = sum(Abundance*100), Treatment = unique(pub.hist), time = unique(time)) %>%
  ungroup()
tab2 <- physeq.gut.d28.d33.d39.abund %>% 
  filter(group == "gut" & time == "d39") %>% 
  mutate(Phylum = fct_explicit_na(Phylum, na_level = "Unassigned")) %>% #missing phylum in one variable (mdf1$Phylum[mdf1$Phylum == "NA"])
  group_by(Phylum, Sample) %>% #group by cond then sample then phylum
  summarise(Abundance = sum(Abundance*100), Treatment = unique(pub.hist), time = unique(time)) %>%
  ungroup()
#phylum <1%
tab1$Phylum <- as.character(tab1$Phylum) 
tab1$Phylum[tab1$Abundance < 1] <- "<1%" 
tab2$Phylum <- as.character(tab2$Phylum) 
tab2$Phylum[tab2$Abundance < 1] <- "<1%" 
#rearrange factor levels
tab1$Phylum <- fct_rev(fct_reorder(tab1$Phylum, tab1$Abundance, mean)) %>%
  fct_drop() #automatically relevel according to current level
tab2$Phylum <- fct_rev(fct_reorder(tab2$Phylum, tab2$Abundance, mean)) %>%
  fct_drop() #automatically relevel according to current level
#color plot with associated phylum names
fig.col <- c(brewer.pal(11,"Set3"),"wheat3","aliceblue", brewer.pal(10, "Paired"), brewer.pal(3, "Set2"))
names(fig.col) <- c("Fusobacteriota","Firmicutes","Proteobacteria","Bacteroidota","Unassigned","Planctomycetota",
                     "Actinobacteriota","Cyanobacteria","Verrucomicrobiota","Myxococcota",
                     "Chloroflexi","Dependentiae","<1%","Chlamydiae","Dadabacteria",
                     "Acidobacteria","Armatimonadetes","Spirochaetes","Deinococcus-Thermus",
                     "Rokubacteria","Patescibacteria","Omnitrophicaeota","Nitrospirae","Gemmatimonadetes","Fibrobacteres",
                     "Elusimicrobia")

#plot1
p1 <- ggplot(tab1, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.8) +
  scale_y_continuous(name = "Relative abundance (%)", expand = c(0,0)) +
  scale_fill_manual(values = fig.col, breaks = names(fig.col)) +
  facet_wrap(~Treatment, nrow = 1, scales = "free_x") +
  theme_classic2() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "none",
        strip.background = element_blank())

#plot2
p2 <- ggplot(tab2, aes(x = Sample, y = Abundance, fill = Phylum)) +
  geom_bar(stat = "identity", position = "stack", color = "black", width = 0.8) +
  scale_y_continuous(name = "Relative abundance (%)", expand = c(0,0)) +
  scale_fill_manual(values = fig.col, breaks = names(fig.col)) +
  facet_wrap(~Treatment, nrow = 1, scales = "free_x") +
  theme_classic2() +
  theme(axis.title.x = element_blank(), axis.text.x = element_blank(),
        axis.ticks.x = element_blank(), legend.position = "right",
        strip.background = element_blank())

#>>>> Figure 5a-b ####


#* PCoA ----

#table
physeq.gut.d28.d33.d39

#distance matrix
pcoa.w <- ordinate(physeq = physeq.gut.d28.d33.d39, method = "PCoA", distance = "wunifrac") # weighted unifrac

#PCoA weighted unifrac
p1 <- plot_ordination(physeq = physeq.gut.d28.d33.d39, ordination = pcoa.w, color = "pub.treat",
                    shape = "pub.treat") + 
  theme_pubr() +
  geom_point(mapping = aes(stroke = 1), size = 4) + #increase circle point width with stroke
  scale_shape_manual(name = "Treament", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = c(19,19,19,19,19,1,1)) +
  scale_color_manual(name = "Treament", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = col.gut.d28.d33.d39) +
  labs(title = "Gut microbiota (d28, d33, d39)", colour = "Treatment") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        legend.position = "right", legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        plot.title = element_text(size = rel(1.3), face = "bold.italic", 
                                  colour = "gray20", hjust = 0.5),
        panel.border = element_rect(colour = "black", fill = NA)) +
  guides(color = guide_legend(override.aes = list(stroke = 1))) #increase circle point width in legend with stroke
#remove the overlap layer with point (allow to use circle point)
p1$layers <- p1$layers[-1]

#>>>> Sup-figure 4a ####
#ggsave("Sup-figure5.pdf", p, "pdf", "../../4_Paper/figures", width = 7, height = 5)


#* Permanova ####

#weighted UniFrac distance
set.seed(1473)
dist <- phyloseq::distance(physeq.gut.d28.d33.d39, method = "wunifrac")

#permanova (comparisons of intra and inter-group variatiosn)
vegan::adonis(dist ~ pub.treat, data = metadata.gut.d28.d33.d39, permutations = 999) #p-value = 0.001

#pairwise comparisons (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, metadata.gut.d28.d33.d39$pub.treat, nperm = 999,
                                    p.method = "bonferroni")

#* Permdisp ####

#dispersion (comparisons of intra-group variations)
disp <- vegan::betadisper(dist, group = metadata.gut.d28.d33.d39$pub.treat)
permutest(disp, pairwise = T)

#* LEfSe d33-d39 ----

# Tables

#metadata
metadata.d33.d39 <- metadata %>% 
  dplyr::filter((time == "d33" | time == "d39") & group == "gut")
#factor
metadata.d33.d39$pub.treat <- factor(metadata.d33.d39$pub.treat)
levels(metadata.d33.d39$pub.treat) <- c("d33","d39")

#no level (last level displayed) + new ASV names
LEfSe.table <- rar_df %>% 
  select(Taxon, ASV, metadata.d33.d39$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "") %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  filter(rowSums(.[,-c(1:2)] > 0) > 0) #check colSums(LEfSe.table[,-c(1:2)])

#color
col.gut.d33.d39

#(Export table) ####

#add the ASV name after the Taxon id
LEfSe.Galaxy <- LEfSe.table %>%
  unite(Taxon, c("Taxon","ASV"))
#export metadata first
LEfSe <- as.data.frame(t(metadata.d33.d39[,c(6,1)])) #test LEfSe between treatment d33(water) d39(100) so choose column treatment
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t33-t39-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = T)
#asv table then added to the previous file with only metadata
LEfSe <- as.data.frame(LEfSe.Galaxy)
#write.table(LEfSe, "2_R/output/LEfSe-galaxy-gut-t33-t39-species-asv-class-sample.txt", sep = "\t", quote = F, 
#            col.names = F, row.names = F, append = T)

# Results

#asv retained from LEfSe analysis on Galaxy (LDA threshold > ??)
LEfSe.asv <- c("asv1644","asv1536","asv2363","asv2138","asv320","asv2362","asv1538",
               "asv2056","asv1029","asv777","asv780","asv1030","asv2506","asv2458",
               "asv406","asv2400","asv2042","asv2577","asv1543","asv2350","asv1662")

#* (Dot plot) ----

#table with cleaned taxa names
LEfSe <- rar_df %>% 
  dplyr::select(Taxon, ASV, metadata.d33.d39$publish) %>% 
  mutate_if(is.double, abund.col) %>% 
  mutate_if(is.character, str_replace_all, pattern = " ", replacement = "")
#asv retained from LEfSe analysis on Galaxy (LDA threshold > 3.5)
LEfSe.asv
#table
LEfSe.table <- LEfSe %>% 
  dplyr::filter(ASV %in% LEfSe.asv) %>%
  mutate_if(is.character, str_replace_all, pattern = ".*;", replacement = "") %>% #remove everthing before the last ;
  mutate_if(is.character, str_replace_all, pattern = "\\.", replacement = "") %>% #remove any "." in Taxa names (specify with \\)
  mutate_if(is.character, str_replace_all, pattern = "__", replacement = ". ") #remove __
#factor
LEfSe.table$ASV <- factor(LEfSe.table$ASV,
                          levels = rev(c("asv1644","asv1536","asv2363","asv2138","asv320","asv2362","asv1538",
                                         "asv2056","asv1029","asv777","asv780","asv1030","asv2506","asv2458",
                                         "asv406","asv2400","asv2042","asv2577","asv1543","asv2350","asv1662")))

LEfSe.table <- LEfSe.table %>% 
  arrange(ASV) %>% 
  dplyr::mutate(taxon = forcats::as_factor(paste0(Taxon," (",ASV,")"))) %>%
  dplyr::select(taxon, everything()) %>%
  dplyr::select(-Taxon,-ASV) %>% 
  pivot_longer(cols = !taxon, names_to = "sample") %>% 
  left_join(., metadata.d33.d39, by = c("sample"="publish")) %>%
  select(taxon, sample, value, pub.treat) %>% 
  group_by(pub.treat, taxon) %>%
  summarise(min = min(value, na.rm = T),
            mean = mean(value, na.rm = T),
            max = max(value, na.rm = T)) 

#change names with uncultured and choose the last level
LEfSe %>% dplyr::filter(ASV %in% LEfSe.asv) %>% dplyr::select(Taxon, ASV)
LEfSe %>% dplyr::filter(ASV == "asv2350") %>% dplyr::select(ASV, Taxon)
#legend text
labs <- c("g. *ZOR0006* (ASV1662)","f. Barnesiellaceae (ASV2350)","g. *Epulopiscium* (ASV1543)",   
          "o. Enterobacterales (ASV2577)","g. *Romboutsia* (ASV2042)","g. *Schlesneria* (ASV2400)",
          "f. Rhizobiaceae (ASV406)","g. *Bosea* (ASV2458)","g. *Reyranella* (ASV2506)",  
          "g. *Reyranella* (ASV1030)","g. *Babeliaceae* (ASV780)","s. *uncultured planctomycete* (ASV777)",
          "g. *Devosia* (ASV1029)","f. Rhodobacteraceae (ASV2056)","g. *Legionella* (ASV1538)",             
          "s. *Alsobacter sp* (ASV2362)","g. *Methylopila* (ASV320)","g. *Microcystis PCC-7914* (ASV2138)",   
          "g. *Reyranella* (ASV2363)","g. *Xanthobacter* (ASV1536)","d. Bacteria (ASV1644)")
names(labs) <- levels(LEfSe.table$taxon)

##legend text
#labs <- c("s. *uncultured bacterium* (asv1662)","s. *uncultured bacterium* (asv2350)","s. *uncultured bacterium* (asv1543)",   
#          "o. Enterobacterales (asv2577)","g. *Romboutsia* (asv2042)","g. *Schlesneria* (asv2400)",
#          "f. Rhizobiaceae (asv406)","g. *Bosea* (asv2458)","g. *Reyranella* (asv2506)",  
#          "g. *Reyranella* (asv1030)","s. *uncultured organism* (asv780)","s. *uncultured planctomycete* (asv777)",
#          "g. *Devosia* (asv1029)","f. Rhodobacteraceae (asv2056)","g. *Legionella* (asv1538)",             
#          "s. *Alsobacter sp* (asv2362)","g. *Methylopila* (asv320)","g. *Microcystis PCC-7914* (asv2138)",   
#          "g. *Reyranella* (asv2363)","g. *Xanthobacter* (asv1536)","d. Bacteria (asv1644)")
#names(labs) <- levels(LEfSe.table$taxon)

#plot
library(ggtext)
p <- ggplot(data = LEfSe.table, aes(x = mean, y = taxon, color = pub.treat)) +
  theme_pubr() +
  geom_pointrange(aes(xmin=min, xmax=max), size = 0.8,
                  position = position_dodge(width = 0.6)) +
  scale_color_manual(values = col.gut.d33.d39) +
  scale_y_discrete(labels = labs) +
  theme(legend.position = "right", 
        legend.text = element_text(size = 13),
        legend.title = element_text(size = 14), 
        panel.border = element_rect(colour = "black", fill = NA),
        axis.line = element_blank(),
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank(),
        axis.text.y = element_markdown(size = 13),
        axis.text.x = element_text(size = 13),
        axis.title.x = element_text(size = 14)) +
  geom_vline(xintercept = 10, linetype = "dashed", color = "gray60") +
  labs(x = "Relative abundance (%)", col = "Treatment")
#>>>> Sup-figure 6 ----
#ggsave("Sup-figure6.pdf", p, "pdf","../../4_Paper/figures", height = 7, width = 8)


#* 4 - GUT d28_0-100/d33/d39 ----

#* Tables ----

#metadata
metadata.gut.d28.0.100.d33.d39 <- metadata.gut.d28.d33.d39 %>%
  dplyr::filter((pub.hist == "d28_0" | pub.hist == "d28_100") & 
                  time == "d28" | time == "d33" | time == "d39")
#factor
metadata.gut.d28.0.100.d33.d39$pub.hist <- factor(metadata.gut.d28.0.100.d33.d39$pub.hist)

#ASV table
gut.d28.0.100.d33.d39 <- rar_df %>% 
  dplyr::select(ASV, metadata.gut.d28.0.100.d33.d39$publish) %>% 
  mutate_if(is.double, abund.col) #check colSums(gut.d28.0.100.d33.d39[,-1])

#color
col.gut.d28.0.100.d33.d39


#* Signif. ASVs ####

#Dominant ASVs previously described
Dom.asv <- c("asv1349","asv1564","asv1620","asv1644","asv1662","asv2363")
#Other discriminant ASVs previously described
Dis.asv <- c("asv475","asv1543","asv1995","asv2350")
Dis.asv2 <- c("asv475","asv1030","asv1406","asv1543","asv1995","asv2161","asv2350")
#others
Signif.network.asv.2

# Proportions among treatments

#subset each ASV and calculate their mean/sd in each treatment
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV %in% Dom.asv]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(ASV, pub.hist) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100)) -> tmp
#asv2363
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv2363"]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#vibrio
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv1406"]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#vibrio
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv2161"]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#reyranella
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv1030"]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))
#reyranella
physeq.gut.d28.d33.d39.abund %>%
  dplyr::filter(OTU %in% asv.otu$OTUID[asv.otu$ASV == "asv2363"]) %>% 
  left_join(., asv.otu, by = c("OTU" = "OTUID")) %>% 
  group_by(pub.treat) %>%
  summarise(mean = mean(Abundance*100), sd = sd(Abundance*100))


#* Dominant ASV ----

#table
tmp <- gut.d28.0.100.d33.d39 %>% 
  filter(ASV %in% Dom.asv) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28.0.100.d33.d39, by = c("sample"="publish"))
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1349","asv1564","asv1620","asv1644","asv1662","asv2363"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1349**<br>(g. *Flavobacterium*)","**ASV1564**<br>(g. *Aeromonas*)",
          "**ASV1620**<br>(g. *Cetobacterium*)","**ASV1644**<br>(d. Bacteria)",
          "**ASV1662**<br>(g. *ZOR0006*)","**ASV2363**<br>(g. *Reyranella*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels
rar_df %>% dplyr::filter(ASV %in% Dom.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p1 <- ggplot(data = tmp, aes(x = pub.hist, y = abund, fill = pub.hist)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ time, scales = "free", space = "free_x", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28.0.100.d33.d39) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)

#>>>> Figure 6a ----
#ggsave("Figure6a",p1,"pdf","../../4_Paper/figures/Figure6a.pdf", width = 6, height = 8)


#* Discriminant ASV ----

#table
tmp <- gut.d28.0.100.d33.d39 %>% 
  filter(ASV %in% Dis.asv) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28.0.100.d33.d39, by = c("sample"="publish"))
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv475","asv1543","asv1995","asv2350"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV475**<br>(g. *Shewanella*)","**ASV1543**<br>(g. *Epulopiscium*)",
          "**ASV1995**<br>(g. *Shewanella*)","**ASV2350**<br>(f. Barnesiellaceae)")
names(labs) <- levels(tmp$ASV)
#check taxa levels
rar_df %>% dplyr::filter(ASV %in% Dis.asv) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p3 <- ggplot(data = tmp, aes(x = pub.hist, y = abund, fill = pub.hist)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ time, scales = "free", space = "free_x", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28.0.100.d33.d39) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)
#>>>> Sup-figure 7a ----
#ggsave("Figure7a",p1,"pdf","../../4_Paper/figures/Figure6a.pdf", width = 6, height = 8)


#* Discriminant ASV ----

#table
tmp <- gut.d28.0.100.d33.d39 %>% 
  filter(ASV %in% Dis.asv2) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28.0.100.d33.d39, by = c("sample"="publish"))
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv475","asv1030","asv1406","asv1543","asv1995",
                                      "asv2161","asv2350"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV475**<br>(g. *Shewanella*)","**ASV1030**<br>(g. *Reyranella*)",
          "**ASV1406**<br>(g. *Vibrio*)","**ASV1543**<br>(g. *Epulopiscium*)",
          "**ASV1995**<br>(g. *Shewanella*)","**ASV2161**<br>(g. *Vibrio*)",
          "**ASV2350**<br>(f. Barnesiellaceae)")
names(labs) <- levels(tmp$ASV)
#check taxa levels
rar_df %>% dplyr::filter(ASV %in% Dis.asv2) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p3 <- ggplot(data = tmp, aes(x = pub.hist, y = abund, fill = pub.hist)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ time, scales = "free", space = "free_x", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28.0.100.d33.d39) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 8.5),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 8.5),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.2)
p3
#>>>> Sup-figure 7a ----
#ggsave("Figure7a",p1,"pdf","../../4_Paper/figures/Figure6a.pdf", width = 6, height = 8)



#* Other ASV plots ----

#table
tmp <- gut.d28.0.100.d33.d39 %>% 
  filter(ASV %in% Signif.network.asv.2) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28.0.100.d33.d39, by = c("sample"="publish"))
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1349","asv1536","asv1941","asv1662","asv1543","asv2042"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1349**<br>(g. *Flavobacterium*)","**ASV1536**<br>(g. *Xanthobacter*)",
          "**ASV19415**<br>(o. *Rhizobiales*)","**ASV1662**<br>(g. *ZOR0006*)",
          "**ASV1543**<br>(g. *Epulopiscium*)","**ASV2042**<br>(g. *Romboutsia*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels
rar_df %>% dplyr::filter(ASV %in% Signif.network.asv.2) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p3 <- ggplot(data = tmp, aes(x = pub.hist, y = abund, fill = pub.hist)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ time, scales = "free", space = "free_x", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28.0.100.d33.d39) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)
p3
#>>>> Sup-figure X ----


#other table
Signif.network.asv.3

#table
tmp <- gut.d28.0.100.d33.d39 %>% 
  filter(ASV %in% Signif.network.asv.3) %>%
  pivot_longer(cols = !ASV, names_to = "sample", values_to = "abund") %>% 
  left_join(., metadata.gut.d28.0.100.d33.d39, by = c("sample"="publish"))
#factor
tmp$ASV <- factor(tmp$ASV, levels = c("asv1406","asv2161","asv1030","asv2363"))

#legend text (Rmarkdown symbols used to put names in *italic* or **bold**)
labs <- c("**ASV1406**<br>(g. *Vibrio*)","**ASV2161**<br>(g. *Vibrio*)",
          "**ASV1030**<br>(g. *Reyranella*)","**ASV2363**<br>(g. *Reyranella*)")
names(labs) <- levels(tmp$ASV)
#check taxa levels
rar_df %>% dplyr::filter(ASV %in% Signif.network.asv.3) %>% select(ASV, Taxon)

#plot with grid
library(ggtext)
p3 <- ggplot(data = tmp, aes(x = pub.hist, y = abund, fill = pub.hist)) +
  theme_pubr() +
  geom_boxplot(alpha = 0.6, outlier.shape = NA) +
  facet_grid(ASV ~ time, scales = "free", space = "free_x", labeller = labeller(ASV = labs)) +
  scale_fill_manual(values = col.gut.d28.0.100.d33.d39) +
  labs(y = "Relative abundance (%)") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8),
        axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.text.y = element_text(size = 9),
        axis.title.y = element_text(size = 11),
        strip.text.y = element_markdown(size = 9),
        strip.text.x = element_blank()) +
  geom_jitter(shape = 16, position = position_jitter(0.2), 
              alpha = 0.2, color = "dimgrey", size = 1.1)
p3
#>>>> Sup-figure XX      ----





#* d28_0-100 Vs d33-d39 ----

#common discriminant ASVs between the two LEfSe analyses (d28_0/d28_100 and d33/d39)

asv.d33.d39 <- c("asv1644","asv1536","asv2363","asv2138","asv320","asv2362","asv1538",
                 "asv2056","asv1029","asv777","asv780","asv1030","asv2506","asv2458",
                 "asv406","asv2400","asv2042","asv2577","asv1543","asv2350","asv1662")
asv.d28.0.100 <- c("asv2363","asv1564","asv1349","asv1493",
                   "asv1941","asv2400","asv1402","asv253","asv857",
                   "asv789","asv2190","asv1995","asv2506","asv1536",
                   "asv1598","asv1543","asv2042","asv1936","asv1817",
                   "asv2577","asv2379","asv1554","asv1046","asv1419",
                   "asv1320","asv1662")
length(asv.d33.d39)
length(asv.d28.0.100)
intersect(asv.d33.d39,asv.d28.0.100)




#*******************************
# 4 - GUT d28/d33/d39 0-100 ####
#*******************************

#* Tables ----

#physoleq file
physeq.gut.0.100 <- subset_samples(physeq.gut.d28.d33.d39,
                   historic %in% c("0","0_0","0_0_100","100","100_0","100_0_100"))
#factor
sample_data(physeq.gut.0.100)$pub.treat <- as_factor(sample_data(physeq.gut.0.100)$pub.treat)
levels(sample_data(physeq.gut.0.100)$pub.treat) <- c("d28_0", "d28_100", "d33 (0-100)", "d39 (0-100)")

#color
col.gut.0.100

#* PCoA ----

#distance matrix
gut.pcoa.w <- ordinate(physeq = physeq.gut.0.100, method = "PCoA", distance = "wunifrac") # weighted unifrac

#PCoA weighted unifrac (1 variable)
p3 <- plot_ordination(physeq = physeq.gut.0.100, ordination = gut.pcoa.w,
                      color = "pub.treat", shape = "pub.treat") + 
  theme_pubr() +
  geom_point(mapping = aes(stroke = 1), size = 3) + #increase circle point width with stroke
  scale_shape_manual(values = c(19,19,1,1)) + #c(3,4,17,16)
  scale_color_manual(values = col.gut.0.100) +
  labs(title = "GUT MICROBIOME") +
  theme(legend.position = "none",
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = rel(1), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  stat_ellipse(type = "t", level = 0.95) +
  guides(color = guide_legend(override.aes = list(linetype = 0, stroke = 1))) #remove lines beacause of ellipse
#remove the overlap layer with point (allow to use circle point)
p3$layers <- p3$layers[-1]

#>>>> Figure 5a-b-c ####

#* Permanova ####

#weighted UniFrac distance
set.seed(1473)
dist <- phyloseq::distance(physeq.gut.0.100, method = "wunifrac")

#permanova (comparisons of intra and inter-group variatiosn)
vegan::adonis(dist ~ pub.treat, data = metadata.gut.0.100, permutations = 999) #p-value = 0.001

#pairwise comparisons (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, metadata.gut.0.100$pub.treat, nperm = 9999,
                                    p.method = "bonferroni")

#* Permdisp ####

#dispersion (comparisons of intra-group variations)
disp <- vegan::betadisper(dist, group = metadata.gut.0.100$pub.treat)
permutest(disp, pairwise = T)


#***************** 
# ------------------
# ------------------
# METABOLOMIC ----
#*****************

#* Tables ----

#metadata for independant analysis of guts with only samples for the publication
metadata.gut.d28
metadata.gut.d28.d33.d39

#metadata for common analysis between tissues
metadata.tissue <- read_csv(file = "../4_analyses_metabo/R/data/metadata_tissues_common.csv") %>% 
  mutate(pub.treat = paste(.$time, .$treatment, sep = "_")) %>% #new column for PCoA
  mutate(pub.treat = replace(pub.treat, pub.treat == "d33_0", "d33")) %>% #replace values of column names
  mutate(pub.treat = replace(pub.treat, pub.treat == "d39_100", "d39"))
metadata.tissue.d28 <- metadata.tissue %>% 
  dplyr::filter(time == "d28")
#factor
metadata.tissue$pub.treat <- as_factor(metadata.tissue$pub.treat)
metadata.tissue.d28$pub.treat <- as_factor(metadata.tissue.d28$pub.treat)

#color
col.comp
col.gut.d28
col.gut.d28.d33.d39


# 1 - GUT d28 ####
#*****************

#* Tables ----
setwd("~/Documents/THESE/3_Projets/1_Medaka_28j/3_Analyses/3_analyses_microbio/")

#metadata with only samples for the publication
metadata.gut.d28

#data matrix (carefull with doublons in pepmass column)
df.metabo <- read_delim("../4_analyses_metabo/R/data/gut_t28_t33_t39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.gut.d28$id)) %>% #reorder
  rename_at(vars(`28D1g1`:`28E3g3`), ~ str_replace_all(., setNames(metadata.gut.d28$publish, metadata.gut.d28$id)))
df.metabo %>% 
  filter(rowSums(.[,-1] > 0) > 0) #keep lines with at least 1 value
#scaling
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(df.metabo[,-1], "pareto")
rownames(X) <- df.metabo$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28


#* PCA ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.gut.d28$pub.treat,
          title = "", ind.names = FALSE, legend = TRUE,  ellipse = F, 
          legend.title = "Treatment",
          pch = c(18,18,18,18,18), col.per.group = col.gut.d28) #PCA samples
#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var
#cleaned plot
p1 <- ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.gut.d28$pub.treat), shape = 18, size = 5) +
  theme_bw() +
  scale_color_manual(values = col.gut.d28) +
  labs(title = "Gut metabolome (d28)", x = "PC1: 35%", y = "PC2: 20%",
       col = "Treatment") +
  theme(panel.grid = element_blank(), legend.position = "none",
        plot.title = element_text(size = rel(1.4), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5))
#>>>> Figure 4a-c ----
layout <- '
A
B'
p <- ((p1 / p5)) +
  patchwork::plot_layout(design = layout, guides = 'keep', tag_level = 'new') +
  patchwork::plot_annotation(tag_levels = list(c('a', 'c'))) & 
  theme(plot.tag = element_text(size = 20, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#save
#ggsave("Figure4a-b-c.pdf", p, "pdf", "../../4_Paper/figures", width = 6, height = 8.5)


#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata, by = c("sample"="publish")) %>% 
  dplyr::select(-("id":"aquarium")) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,1:1674]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

##manhattan distance
#set.seed(2)
#dist <- vegdist(df.metabo.activ, method = "manhattan")

#euclidean distance
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#post hoc   
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")
  
#* Permdisp ----
test <- vegan::betadisper(dist, group = tmp$pub.treat)
permutest(test, pairwise = T)


#************************* 
# 2 - GUT d28/d33/d39 ####
#*************************

#* Tables ----

#metadata with only samples for the publication
metadata.gut.d28.d33.d39$pub.treat

#data matrix (carefull with doublons in pepmass column)
df.metabo.d28.d33.d39 <- read_delim("../4_analyses_metabo/R/data/gut_t28_t33_t39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.gut.d28.d33.d39$id)) %>% #reorder
  rename_at(vars(`28D1g1`:`39Eg4`), ~ str_replace_all(., setNames(metadata.gut.d28.d33.d39$publish, 
                                                                  metadata.gut.d28.d33.d39$id)))
#>>>> Sup_Data_4 ----
#write_csv(df.metabo.d28.d33.d39, "../../4_Paper/soumissions/1_Nature_com/data/Supplementary_Data_4.csv")

#scaling
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(df.metabo.d28.d33.d39[,-1], "pareto")
rownames(X) <- df.metabo.d28.d33.d39$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28.d33.d39

#* Diversity ----
div <- df.metabo.d28.d33.d39[,-1] %>% t(.)
specnumber(div, groups = metadata.gut.d28.d33.d39$pub.treat)


#* PCA (unsupervised) ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.gut.d28.d33.d39$pub.treat,
          title = "", ind.names = FALSE, legend = TRUE,  ellipse = T, 
          legend.title = "Treatment")

#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var
#cleaned plot
p2 <- ggplot(tmp, mapping = aes(x = PC1, y = PC2, col = metadata.gut.d28.d33.d39$pub.treat,
                                shape = metadata.gut.d28.d33.d39$pub.treat)) +
  theme_pubr() +
  geom_point(mapping = aes(stroke = 1), size = 4) +
  scale_shape_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = c(18,18,18,18,18,5,5)) +
  scale_color_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = col.gut.d28.d33.d39) +
  labs(title = "Gut metabolome (d28, d33, d39)", x = "PC1: 33%", y = "PC2: 17%",
       col = "Treatment") +
  theme(axis.text = element_text(size = 13), axis.title = element_text(size = 15),
        legend.position = "none", legend.title = element_text(size = 15),
        legend.text = element_text(size = 13),
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = rel(1.3), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(override.aes = list(stroke = 1))) #remove lines because of ellipse

#>>>> Sup-figure 4b ####
#plot
#p <- (p1 / p2) +
#  patchwork::plot_layout(guides = 'collect') +
#  patchwork::plot_annotation(tag_levels = 'a') & 
#  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
#        plot.tag.position = "topleft")
#save
#ggsave("Sup-figure5a-b.pdf", p1, "pdf", "../../4_Paper/figures", width = 7, height = 8.5)


#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.gut.d28.d33.d39, by = c("sample"="publish")) %>% 
  dplyr::select(-("id":"aquarium")) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,1:1674]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

##manhattan distance
#set.seed(2)
#dist <- vegdist(df.metabo.activ, method = "manhattan")

#euclidean distance
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#post hoc
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")

#* Permdisp ----
test <- vegan::betadisper(dist, group = tmp$pub.treat)
permutest(test, pairwise = T)

  

#*******************************
# 3 - GUT d28/d33/d39 0-100 ####
#*******************************

#* Tables ----

#metadata with only samples for the publication
metadata.gut.0.100 <- metadata %>% 
  dplyr::filter(group == "gut" & (time == "d28" | time == "d33" | time == "d39")) %>% 
  dplyr::filter((historic == "0" | historic == "0_0" | historic == "0_0_100" |
                   historic == "100" | historic == "100_0" | historic == "100_0_100"))
#factor
metadata.gut.0.100$pub.treat <- factor(metadata.gut.0.100$pub.treat)
levels(metadata.gut.0.100$pub.treat) <- c("d28_0","d28_100","d33 (0-100)","d39 (0-100)")

#data matrix only 0 versus 100 (carefull with doublons in pepmass column)
df.metabo.0.100 <- read_delim("../4_analyses_metabo/R/data/gut_t28_t33_t39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.gut.0.100$id)) %>% #reorder
  rename_at(vars(`28D1g1`:`39Eg4`), ~ str_replace_all(., setNames(metadata.gut.0.100$publish, 
                                                                  metadata.gut.0.100$id)))
#scaling after filtering 0 versus 100
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(df.metabo.0.100[,-1], "pareto")
rownames(X) <- df.metabo.0.100$PEPMASS
X <- t(X) #transposing

#color
col.gut.0.100


#* PCA (unsupervised) ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.gut.0.100$pub.treat,
          title = "", ind.names = FALSE, legend = TRUE,  ellipse = T, 
          legend.title = "Treatment", col.per.group = col.gut.0.100,
          pch = c(3,4,17,16))

#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var
#cleaned plot
p4 <- ggplot(tmp, mapping = aes(x = PC1, y = PC2, col = metadata.gut.0.100$pub.treat,
                          shape = metadata.gut.0.100$pub.treat)) +
  theme_pubr() +
  geom_point(mapping = aes(stroke = 1), size = 3) + #increase circle point width with stroke
  scale_shape_manual(name = "Treatment", labels = c("d28_0","d28_100","d33 (0-100)","d39 (0-100)"),
                     values = c(19,19,1,1)) +
  scale_color_manual(name = "Treatment", labels = c("d28_0","d28_100","d33 (0-100)","d39 (0-100)"), 
                     values = col.gut.0.100) +
  labs(title = "GUT METABOLOME", x = "PC1: 34%", y = "PC2: 25%",
       col = "Treatment") +
  theme(legend.position = "right",
        panel.border = element_rect(colour = "black", fill = NA),
        plot.title = element_text(size = rel(1), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  stat_ellipse(type = "t", level = 0.95) +
  guides(color = guide_legend(override.aes = list(linetype = 0, stroke = 1))) #remove lines beacause of ellipse

#>>>> Figure 5a-b-c-d ####

#plot
p <- (p1 / p2 / (p3 | p4)) +
  patchwork::plot_layout(guides = 'keep') +
  patchwork::plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 17, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#save
#ggsave("Figure5a-b-c-d.pdf", p, "pdf", "../../4_Paper/figures", width = 10, height = 10)


#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.gut.0.100, by = c("sample"="publish")) %>% 
  dplyr::select(-("id":"aquarium")) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,1:1674]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

##manhattan distance
#set.seed(2)
#dist <- vegdist(df.metabo.activ, method = "manhattan")

#euclidean distance
set.seed(1473)
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#post hoc
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")

#* Permdisp ----
test <- vegan::betadisper(dist, group = tmp$pub.treat)
permutest(test, pairwise = T)


#*********************
# 4 - LIVER d28 ####
#*********************

#* Tables ----

#metadata (common with muscle)
metadata.tissue.d28 %>% dplyr::filter(group == "liver")

#liver
liv.metabo.d28 <- read_delim("../4_analyses_metabo/R/data/liver_d28_d33_d39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  dplyr::select(c("PEPMASS", metadata.tissue.d28 %>% dplyr::filter(group == "liver") %>% .$id))

#scaling after filtering
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(liv.metabo.d28[,-1], "pareto")
rownames(X) <- liv.metabo.d28$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28


#* PCA ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.tissue.d28 %>% dplyr::filter(group == "liver") %>% .$pub.treat,
          title = "Liver metabolome", ind.names = FALSE, legend = TRUE,  ellipse = F, 
          legend.title = "Treatment",
          pch = c(18,18,18,18,18), col.per.group = col.gut.d28) #PCA samples
#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var

#cleaned plot
p1 <- ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.tissue.d28 %>% filter(group == "liver") %>% .$pub.treat),
             shape = 18, size = 5) +
  theme_bw() +
  scale_color_manual(values = col.gut.d28) +
  labs(title = "Liver metabolome (d28)", x = "PC1: 49%", y = "PC2: 19%",
       col = "Treatment") +
  theme(panel.grid = element_blank(), legend.title = element_text(size = 15),
        legend.text = element_text(size = 13), axis.title = element_text(size = 16),
        axis.text =  element_text(size = 12),legend.position = "right",
        plot.title = element_text(size = rel(1.3), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
#>>>> Sup-figure 3a ####


#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.tissue.d28 %>% filter(group == "liver"), by = c("sample"="id")) %>% 
  dplyr::select(-("group":"treatment")) %>% 
  dplyr::select(sample, pub.treat, everything()) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,-1]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

#euclidean distance
set.seed(1473)
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")


#********************
# 5 - MUSCLE d28 ####
#********************

#* Tables ----

#metadata (common with liver)
metadata.tissue.d28 %>% dplyr::filter(group == "muscle")

#muscle
mus.metabo.d28 <- read_delim("../4_analyses_metabo/R/data/muscle_d28_d33_d39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  dplyr::select(c("PEPMASS", metadata.tissue.d28 %>% dplyr::filter(group == "muscle") %>% .$id))

#scaling after filtering
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(mus.metabo.d28[,-1], "pareto")
rownames(X) <- mus.metabo.d28$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28


#* PCA ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.tissue.d28 %>% dplyr::filter(group == "muscle") %>% .$pub.treat,
          title = "Muscle metabolome", ind.names = FALSE, legend = TRUE,  ellipse = F, 
          legend.title = "Treatment",
          pch = c(18,18,18,18,18), col.per.group = col.gut.d28) #PCA samples
#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var

#cleaned plot
p2 <- ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.tissue.d28 %>% filter(group == "muscle") %>% .$pub.treat),
             shape = 18, size = 5) +
  theme_bw() +
  scale_color_manual(values = col.gut.d28) +
  labs(title = "Muscle metabolome (d28)", x = "PC1: 37%", y = "PC2: 18%",
       col = "Treatment") +
  theme(panel.grid = element_blank(), legend.title = element_text(size = 15),
        legend.text = element_text(size = 13), axis.title = element_text(size = 16),
        axis.text =  element_text(size = 12),legend.position = "right",
        plot.title = element_text(size = rel(1.3), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
#>>>> Sup-figure 3a-b -----

p <- (p1 / p2) +
  patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#save
#ggsave("Sup-figure3a-b.pdf", p, "pdf", "../../4_Paper/figures", width = 6, height = 8)

#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.tissue.d28 %>% filter(group == "muscle"), by = c("sample"="id")) %>% 
  dplyr::select(-("group":"treatment")) %>% 
  dplyr::select(sample, pub.treat, everything()) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,-1]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

#euclidean distance
set.seed(1473)
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#POST HOC (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")


#***************************
# 6 - LIVER d28 d33 d39 ####
#***************************

#* Tables ----

#metadata
metadata.tissue %>% dplyr::filter(group == "liver")

#liver
liv.metabo.d28.d33.d39 <- read_delim("../4_analyses_metabo/R/data/liver_d28_d33_d39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  dplyr::select(c("PEPMASS", metadata.tissue %>% dplyr::filter(group == "liver") %>% .$id))

#>>>> Sup_Data_5 ----
#write_csv(liv.metabo.d28.d33.d39, "../../4_Paper/soumissions/1_Nature_com/data/Supplementary_Data_5.csv")

#scaling after filtering
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(liv.metabo.d28.d33.d39[,-1], "pareto")
rownames(X) <- liv.metabo.d28.d33.d39$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28.d33.d39


#* PCA ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.tissue %>% dplyr::filter(group == "liver") %>% .$pub.treat,
          ind.names = FALSE, legend = TRUE,  ellipse = F, legend.title = "Treatment") #PCA samples
#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var

#cleaned plot
p3 <- ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.tissue %>% filter(group == "liver") %>% .$pub.treat, stroke = 1,
                           shape = metadata.tissue %>% filter(group == "liver") %>% .$pub.treat), size = 5) +
  theme_bw() +
  scale_shape_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = c(18,18,18,18,18,5,5)) +
  scale_color_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = col.gut.d28.d33.d39) +
  labs(title = "Liver metabolome (d28, d33, d39)", x = "PC1: 41%", y = "PC2: 18%",
       col = "Treatment") +
  theme(panel.grid = element_blank(), legend.title = element_text(size = 15),
        legend.text = element_text(size = 13), axis.title = element_text(size = 16),
        axis.text =  element_text(size = 12),legend.position = "right",
        plot.title = element_text(size = rel(1.3), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
#>>>> Sup-figure 4c ####


#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.tissue %>% filter(group == "liver"), by = c("sample"="id")) %>% 
  dplyr::select(-("group":"treatment")) %>% 
  dplyr::select(sample, pub.treat, everything()) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,-1]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

#euclidean distance
set.seed(1473)
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#POST HOC (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")




#****************************
# 7 - MUSCLE d28 d33 d39 ####
#****************************

#* Tables ----

#metadata (common with liver)
metadata.tissue %>% dplyr::filter(group == "muscle")

#muscle
mus.metabo.d28.d33.d39 <- read_delim("../4_analyses_metabo/R/data/muscle_d28_d33_d39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  dplyr::select(c("PEPMASS", metadata.tissue %>% dplyr::filter(group == "muscle") %>% .$id))

#>>>> Sup_Data_6 ----
#write_csv(mus.metabo.d28.d33.d39, "../../4_Paper/soumissions/1_Nature_com/data/Supplementary_Data_6.csv")

#scaling after filtering
library(MetabolAnalyze)
X <- MetabolAnalyze::scaling(mus.metabo.d28.d33.d39[,-1], "pareto")
rownames(X) <- mus.metabo.d28.d33.d39$PEPMASS
X <- t(X) #transposing

#color
col.gut.d28.d33.d39

  
#* PCA ----

detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#how many PC for the final analysis ?
tune.pca(X, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- pca(X, ncomp = 3, center = T, scale = F)
#plot
plotIndiv(pca.metabo, group = metadata.tissue %>% dplyr::filter(group == "muscle") %>% .$pub.treat,
          ind.names = FALSE, legend = TRUE,  ellipse = F, legend.title = "Treatment") #PCA samples
#export coordinates to create plot using ggplot
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
var.exp <- pca.metabo$prop_expl_var

#cleaned plot
p4 <- ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.tissue %>% filter(group == "muscle") %>% .$pub.treat, stroke = 1,
                           shape = metadata.tissue %>% filter(group == "muscle") %>% .$pub.treat), size = 5) +
  theme_bw() +
  scale_shape_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = c(18,18,18,18,18,5,5)) +
  scale_color_manual(name = "Treatment", labels = c("d28_0","d28_Z8","d28_1","d28_10","d28_100",
                                                    "d33","d39"), values = col.gut.d28.d33.d39) +
  labs(title = "Muscle metabolome (d28, d33, d39)", x = "PC1: 35%", y = "PC2: 21%",
       col = "Treatment") +
  theme(panel.grid = element_blank(), legend.title = element_text(size = 15),
        legend.text = element_text(size = 13), axis.title = element_text(size = 16),
        axis.text =  element_text(size = 12),legend.position = "right",
        plot.title = element_text(size = rel(1.3), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))
#>>>> Sup-figure 4a-b-c-d ####

#plot
p <- ((p1 | p2) / (p3 | p4)) +
  patchwork::plot_layout(guides = 'collect') +
  patchwork::plot_annotation(tag_levels = 'a') & 
  theme(plot.tag = element_text(size = 18, hjust = 0, face = "bold"),
        plot.tag.position = "topleft")
#save
#  ggsave("Sup-figure4a-b-c-d.pdf", p, "pdf", "../../4_Paper/figures", width = 12, height = 10)

 
#* Statistical analyses ----

#whole table normalized
tmp <- X %>% 
  as_tibble(.)
tmp$sample <- rownames(X)
tmp <- tmp %>%
  dplyr::select(sample, everything()) %>% 
  left_join(., metadata.tissue %>% filter(group == "muscle"), by = c("sample"="id")) %>% 
  dplyr::select(-("group":"treatment")) %>% 
  dplyr::select(sample, pub.treat, everything()) %>% 
  column_to_rownames(., "sample")
#active table (no metadata)
df.metabo.activ <- tmp[,-1]

#* Permanova ----
library(vegan)
library(RVAideMemoire)

#euclidean distance
set.seed(1473)
dist <- vegan::vegdist(df.metabo.activ, method = "euclidean")
#permanova
vegan::adonis(dist ~ pub.treat, data = tmp, permutations = 999, method = "euclidean")
#POST HOC (adjust. Bonferroni)
RVAideMemoire::pairwise.perm.manova(dist, tmp$pub.treat, nperm = 999, p.method = "bonferroni")




#****************************
# 8 - TISSUES d28 ----

#* Tables ----


df.tissue <- read_csv("../4_analyses_metabo/R/data/tissues_d28_muscle_gut_liver.csv") %>% 
  dplyr::select("PEPMASS","NAME_METABOSCAPE", metadata.tissue.d28$id) %>% #remove the fish 4 from guts and liver because muscle has only 3
  filter(rowSums(.[,-2] > 0) > 0)
#pareto transformation
df.par <- MetabolAnalyze::scaling(df.tissue[,-(1:2)], "pareto")
df.par$PEPMASS <- df.tissue$PEPMASS
df.par <- df.par %>% 
  tibble(.) %>% 
  dplyr::select(PEPMASS, everything())

#all for PCA
df.par <- df.par %>% 
  pivot_longer(., cols = !PEPMASS, names_to = "id") %>% #transpose step 1
  pivot_wider(., names_from = "PEPMASS", values_from = "value") #transpose step 2
df.par <- data.frame(df.par, row.names = "id", check.names = F)

#Colors
col.tissue

#* PCA ----
detach("package:MetabolAnalyze", unload = TRUE)
library(mixOmics)

#PCA: how many PC for the final analysis ?
mixOmics::tune.pca(df.par, ncomp = 10, center = TRUE, scale = FALSE)
#pca with 3 PC
pca.metabo <- mixOmics::pca(df.par, ncomp = 4, center = T, scale = F)
#PCA plot
mixOmics::plotIndiv(pca.metabo, group = metadata.tissue.d28$group,
                    title = "", ind.names = FALSE, legend = TRUE,  ellipse = F, 
                    legend.title = "Tissue",
                    pch = c(18,18,18,18,18), col.per.group = col.tissue)
#export coordinates
tmp <- pca.metabo$variates$X
tmp <- as.data.frame(tmp)
pca.metabo$prop_expl_var

#change plot
ggplot(tmp) +
  geom_point(mapping = aes(x = PC1, y = PC2, col = metadata.tissue.d28 %>% .$group),
             shape = 18, size = 5) +
  theme_bw() +
  scale_color_manual(values = col.tissue) +
  labs(x = "PC1: 35%", y = "PC2: 20%", col = "Tissue") +
  theme(panel.grid = element_blank(), legend.title = element_text(size = 15),
        legend.text = element_text(size = 13), axis.title = element_text(size = 16),
        axis.text =  element_text(size = 12), legend.position = "right",
        plot.title = element_text(size = rel(1.5), face = "bold.italic", colour = "gray20",
                                  hjust = 0.5)) +
  guides(color = guide_legend(ncol = 1, byrow = TRUE))



# 9 - CULTURE EXOMETABOLITES ----

#* Tables ####

#metadata 

#metadata with only common samples for the three data sets
metadata.exomet <- read_csv(file = "../4_analyses_metabo/R/data/metadata-cult-biof-gut.csv")

##data matrix (carefull with doublons in pepmass column)
exomet <- read_csv("../4_analyses_metabo/R/data/culture_gutd0_biofilmd28_pepmass.csv") %>%
  dplyr::select("PEPMASS", "NAME_METABOSCAPE", metadata.exomet$id) %>% #reorder
  rename_at(vars(`28D1b`:`0Eg1`), ~ str_replace_all(., setNames(metadata.exomet$publish, metadata.exomet$id))) %>%
  filter(rowSums(.[,-2] > 0) > 0)

#color
col.comp[c(1,3,4)]

#* Venn Diagram ####

#lists
l.gut <- exomet %>% 
  dplyr::select("PEPMASS", metadata.exomet %>% dplyr::filter(tissue == "gut") %>% .$publish) %>% 
  filter(rowSums(.[,-1] > 0) > 0) %>% 
  .$PEPMASS
length(l.gut)

l.biof <- exomet %>% 
  dplyr::select("PEPMASS", metadata.exomet %>% dplyr::filter(tissue == "biofilm") %>% .$publish) %>% 
  filter(rowSums(.[,-1] > 0) > 0) %>% 
  .$PEPMASS
length(l.biof)

l.cult <- exomet %>% 
  dplyr::select("PEPMASS", metadata.exomet %>% dplyr::filter(tissue == "culture") %>% .$publish) %>% 
  filter(rowSums(.[,-1] > 0) > 0) %>% 
  .$PEPMASS
length(l.cult)

#intersection between the 5 species
inter <- intersect(l.gut, l.biof)
inter <- intersect(inter, l.cult)
length(inter) #569 shared ASVs

#list
x <- list("Gut" = l.gut,
          "Culture" = l.cult,
          "Biofilm" = l.biof)

##VennDiagram
#VennDiagram::venn.diagram(x, cat.dist = c(0.04,0.04,0.04), cat.pos = c(-4,5,180),
#                          fill = col.comp, col = col.comp, cex = 1.8, lty = 'blank',
#                          sub.cex = 2, cat.cex = 1.8,
#                          #cat.fontface = c("plain","italic","italic","italic","italic"),
#                           #main = "12,313 total ASVs (380,236 reads)\n32 shared ASVs (34,715 reads)",
#                          filename = "figure/gut-d0_biof-d28_cult/venn.tiff") 


# Occurence in exposed guts ----

#* Tables ----

#specific metabolites in culture (n = 80)
cult <- exomet %>% 
  dplyr::filter(PEPMASS %in% setdiff(exomet$PEPMASS,union(l.gut,l.biof)))
#culture exometabolites => keep 2 values after comma
cult$PEPMASS <- round(cult$PEPMASS, 2)

#guts
exo.gut <- df.metabo
#keep 2 values after comma
exo.gut$PEPMASS <- round(as.double(exo.gut$PEPMASS), 2)
#change NA with 0
exo.gut <- exo.gut %>% 
  mutate(PEPMASS = coalesce(PEPMASS, 0))

#exposed guts
exo.gut.bloom <- df.metabo %>% 
  dplyr::select(PEPMASS, metadata.gut.d28 %>% dplyr::filter(treatment %in% c("1","10","100")) %>% .$publish)
#keep 2 values after comma
exo.gut.bloom$PEPMASS <- round(as.double(exo.gut.bloom$PEPMASS), 2)
#change NA with 0
exo.gut.bloom <- exo.gut.bloom %>% 
  mutate(PEPMASS = coalesce(PEPMASS, 0))

#livers
liv.metabo.d28
#keep 2 values after comma
liv.metabo.d28$PEPMASS <- round(liv.metabo.d28$PEPMASS, 2)

#muscles
mus.metabo.d28
#keep 2 values after comma
mus.metabo.d28$PEPMASS <- round(mus.metabo.d28$PEPMASS, 2)

#metabolites networks gut d28
metabolites <- round(as.double(metabolites), 2)
metabolites.d33.d39 <- round(as.double(metabolites.d33.d39), 2)



#* Occurence ----

#guts
intersect(exo.gut$PEPMASS, cult$PEPMASS)
#annotated?
cult %>% 
  dplyr::filter(PEPMASS %in% intersect(exo.gut$PEPMASS, cult$PEPMASS)) #non annotated
#in which exposed gut treatments?
tmp <- exo.gut %>% 
  dplyr::filter(PEPMASS %in% intersect(exo.gut$PEPMASS, cult$PEPMASS)) #non annotated


#exposed guts
intersect(exo.gut.bloom$PEPMASS, cult$PEPMASS)
#annotated?
cult %>% 
  dplyr::filter(PEPMASS %in% intersect(exo.gut.bloom$PEPMASS, cult$PEPMASS)) #non annotated
#in which exposed gut treatments?
tmp <- exo.gut.bloom %>% 
  dplyr::filter(PEPMASS %in% intersect(exo.gut.bloom$PEPMASS, cult$PEPMASS)) #non annotated


#livers
intersect(liv.metabo.d28$PEPMASS, cult$PEPMASS)
#annotated?
cult %>% 
  dplyr::filter(PEPMASS %in% intersect(liv.metabo.d28$PEPMASS, cult$PEPMASS)) #non annotated
#in which liver treatments?
tmp <- liv.metabo.d28 %>% 
  dplyr::filter(PEPMASS %in% intersect(liv.metabo.d28$PEPMASS, cult$PEPMASS)) #non annotated


#muscles
intersect(mus.metabo.d28$PEPMASS, cult$PEPMASS)
#annotated?
cult %>% 
  dplyr::filter(PEPMASS %in% intersect(mus.metabo.d28$PEPMASS, cult$PEPMASS)) #non annotated
#in which liver treatments?
tmp <- mus.metabo.d28 %>% 
  dplyr::filter(PEPMASS %in% intersect(mus.metabo.d28$PEPMASS, cult$PEPMASS)) #non annotated


#metabolites network gut d28 and gut d33 d39
intersect(metabolites, cult$PEPMASS)
intersect(metabolites.d33.d39, cult$PEPMASS)


# 10 - CULTURE METABOLITES ----

#* Tables ####

##data
pmc728.11 <- read_csv("../4_analyses_metabo/PMC728.11/Table 728.11 Alison 2021.csv") %>%
  dplyr::select("Annotation","GNPS cluster") %>% 
  dplyr::rename(cluster = "GNPS cluster") %>% 
  filter_all(any_vars(!is.na(.))) %>%
  dplyr::mutate(cluster = replace(cluster, cluster == "Unknown?", 'Unknown')) %>%
  dplyr::mutate(cluster = replace_na(cluster, "Unknown")) %>% 
  dplyr::filter(cluster != "Unknown") %>% 
  count(cluster) %>% 
  dplyr::arrange(desc(n)) %>% 
  dplyr::rename(`Metabolite annotation` = "cluster", Counts = "n")
pmc728.11$`Proportion (%)` <- pmc728.11$Counts*100/sum(pmc728.11$Counts)
pmc728.11 <- pmc728.11 %>% 
  dplyr::select(-Counts)
#>>>> Sup-table1
write_csv(pmc728.11, "../../4_Paper/soumissions/1_Nature_com/figures/Supplementary_Table_1.csv")

  ##factor
#pmc728.11$cluster <- fct_rev(fct_reorder(pmc728.11$cluster, pmc728.11$n)) %>% #reorder each phylum according to mean of their relative abundances
#  fct_drop()
##add columns
#pmc728.11$fraction <- pmc728.11$n / sum(pmc728.11$n)
#pmc728.11$ymax <- cumsum(pmc728.11$fraction)
#pmc728.11$ymin = c(0, head(pmc728.11$ymax, n=-1))
#pmc728.11$labelPosition <- (pmc728.11$ymax + pmc728.11$ymin) / 2
#pmc728.11$label <- paste0(round(pmc728.11$fraction*100, 0), "%") #check sum(pmc728.11$fraction*100)
##col
#col.pmc <- rev(c(RColorBrewer::brewer.pal(n = 7, name = "Purples"),
#             RColorBrewer::brewer.pal(n = 7, name = "Blues"),
#             RColorBrewer::brewer.pal(n = 6, name = "YlGn")))
#names(col.pmc) <- c(levels(pmc728.11$cluster))
###plot donut https://www.r-graph-gallery.com/128-ring-or-donut-plot.html
#ggplot(pmc728.11, aes(ymax=ymax, ymin=ymin, xmax=4, xmin=3, fill=cluster)) +
#  theme_void() +
#  geom_rect() +
#  coord_polar(theta="y") + # Try to remove that to understand how the chart is built initially
#  geom_label(x=3.5, aes(y=labelPosition, label=label), size=6) +
#  scale_fill_manual(values = col.pmc[levels(pmc728.11$cluster)]) +
#  xlim(c(2, 4))


# MULTI-OMICS ----
#*****************

#*****************
# 1 - GUT d28 ####
#*****************

#* Common tables ----

#metadata with only samples for the publication
metadata.gut.d28

#ASV new names for publication
asv.otu

#color
col.gut.d28

#* Tables ----

#ASV table with abundances for filtering (keep ASV if at least 1% in 1 sample)
gut.d28.abundant %>% 
  filter(rowSums(.[,-(1:10)] > 0) > 0) %>% 
  dplyr::select(ASV) -> tmp
#ASV table filtered in counts
tmp <- gut.d28 %>%
  filter(ASV %in% tmp$ASV) %>% 
  mutate(ASV = str_replace(ASV, "asv", "ASV")) #replace asv with ASV
#prepare centered log ratio transformation (CLR)
X <- t(tmp[,-1])
colnames(X) <- tmp$ASV
X <- mixOmics::logratio.transfo(as.matrix(X), logratio = 'CLR', offset = 1) #to prevent log transform later (https://mixomics.org/mixmc/mixmc-pre-processing/)
sum(which(X == 0)) #check if still zeros, if result 0 => it's ok
dim(X)
X <- as.data.frame(X[,]) #add [,] to avoid the error

#metabolites table (filter without some samples : 28A2g4; 28A3g3; 28C2g4; 28E3g4)
df.metabo <- read_delim("../4_analyses_metabo/R/data/gut_t28_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.gut.d28$id)) %>% #reorder
  rename_at(vars(`28D1g1`:`28E3g3`), ~ str_replace_all(., setNames(metadata.gut.d28$publish, metadata.gut.d28$id)))
#pareto transformation
library(MetabolAnalyze)
Y <- MetabolAnalyze::scaling(df.metabo[,-1], "pareto")
rownames(Y) <- df.metabo$PEPMASS
Y <- t(Y)
Y <- as.data.frame(Y)


#* Choose the design matrice ####
## 1- full weighted design matrice

#with this design, all blocks are connected and this full weighted design offers a trade-off between
#maximizing the correlation between datasets and the separation between classes (Y)
design <- matrix(c(0,0.1,1,
                   0.1,0,1,
                   1,1,0), ncol=3)

#* Combine matrices ####

#combine all data tables in a big list
all <- list(microbiome = X,
            metabolome = Y)

#block matrices
block.plsda.medaka1 = mixOmics::block.plsda(X = all, Y = metadata.gut.d28$pub.treat,
                                            ncomp = 3, design = design)
plotIndiv(block.plsda.medaka1)
#correlation first component
plot(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])
cor(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])


#* Network ----

#>>>> Figure 4b ####
#pdf("../../4_Paper/figures/Figure4b.pdf")
mycols <- c(colorRampPalette(c("#6587cd", "#FFFFFF"))(6),
            colorRampPalette(c("#FFFFFF", "#c85a38"))(6))
mycols = mycols[mycols != "#FFFFFF"]
p <- network(block.plsda.medaka1, cutoff = 0.695, #0.695
        color.node = c(alpha("sandybrown",0.9),alpha("white",0.7)), 
        shape.node = c("circle", "rectangle"),
        cex.node.name = 0.9,
        #color.edge = color.spectral(100),
        color.edge = mycols,
        keysize = c(0.83,0.83),
        #color.edge = color.GreenRed(n = 10, alpha = 0.4),
        lty.edge = "solid", lwd.edge = 5, #dotted instead of solid
        show.edge.labels = F,
        symkey = T)
dev.off()
#p2 <- margin(element_blank())
#
#ggarrange(p, p2,
#          widths = c(1, 1),
#          labels = c("a", "b"), font.label = list(size = 17),
#          nrow = 2, align = "hv") %>% 
#  ggexport(filename = "../../../4_Paper/figures/Figure3.pdf", width = 9, height = 12)

#export high correlated metabolites and fin occurence with cultures exometabolites (previous point 9)
metabolites <- p$gR[1]
metabolites <- names(metabolites)[-c(1:6)]
tmp <- as_tibble(metabolites, paste0(seq(1, 46)))
#write_csv(tmp, "../5_analyses_correlation_metabo-microbio/R/output/metabolites_networkd28_for_annotation.csv")


#***********************
# 2 - GUT d28 0-100 ####
#***********************

#* Common tables ----

#metadata with only samples for the publication
metadata.gut.d28.0.100
#factor
metadata.gut.d28.0.100$pub.treat <- factor(metadata.gut.d28.0.100$pub.treat)
#rename levels
levels(metadata.gut.d28.0.100$pub.treat) <- c("d28_0","d28_100")

#ASV new names for publication
asv.otu

#color
col.gut.d28[c(1,5)]

#* Tables ----

#ASV table with abundances for filtering (keep ASV if at least 1% in 1 sample)
gut.d28.abundant %>% 
  dplyr::select(c(ASV, metadata.gut.d28.0.100$publish)) %>%
  filter(rowSums(.[,-1] > 0) > 0) %>% 
  dplyr::select(ASV) -> tmp
#ASV table filtered in counts
tmp <- gut.d28 %>% 
  dplyr::select(c(ASV, metadata.gut.d28.0.100$publish)) %>%
  filter(ASV %in% tmp$ASV) %>% 
  mutate(ASV = str_replace(ASV, "asv", "ASV")) #replace asv with ASV
#prepare centered log ratio transformation (CLR)
X <- t(tmp[,-1])
colnames(X) <- tmp$ASV
X <- mixOmics::logratio.transfo(as.matrix(X), logratio = 'CLR', offset = 1) #to prevent log transform later (https://mixomics.org/mixmc/mixmc-pre-processing/)
sum(which(X == 0)) #check if still zeros, if result 0 => it's ok
dim(X)
X <- as.data.frame(X[,]) #add [,] to avoid the error

#metabolites table (filter without some samples : 28A2g4; 28A3g3; 28C2g4; 28E3g4)
df.metabo <- read_delim("../4_analyses_metabo/R/data/gut_t28_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.gut.d28.0.100$id)) %>% #reorder
  rename_at(vars(`28D1g1`:`28E3g3`), ~ str_replace_all(., setNames(metadata.gut.d28.0.100$publish, metadata.gut.d28.0.100$id)))
#pareto transformation
library(MetabolAnalyze)
Y <- MetabolAnalyze::scaling(df.metabo[,-1], "pareto")
rownames(Y) <- df.metabo$PEPMASS
Y <- t(Y)
Y <- as.data.frame(Y)


#* Choose the design matrice ####
## 1- full weighted design matrice

#with this design, all blocks are connected and this full weighted design offers a trade-off between
#maximizing the correlation between datasets and the separation between classes (Y)
design <- matrix(c(0,0.1,1,
                   0.1,0,1,
                   1,1,0), ncol=3)

#* Combine matrices ####

#combine all data tables in a big list
all <- list(microbiome = X,
            metabolome = Y)

#block matrices
block.plsda.medaka1 = mixOmics::block.plsda(X = all, Y = metadata.gut.d28.0.100$pub.treat,
                                            ncomp = 3, design = design)
#correlation first component
plot(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])
cor(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])

#* Network ----

#>>>> Figure XX ####
#pdf("../../4_Paper/figures/FigureX.pdf")
mycols <- c(colorRampPalette(c("#6587cd", "#FFFFFF"))(6),
            colorRampPalette(c("#FFFFFF", "#c85a38"))(6))
mycols = mycols[mycols != "#FFFFFF"]
network(block.plsda.medaka1, cutoff = 0.798,
        color.node = c(alpha("sandybrown",0.9),alpha("white",0.7)), 
        shape.node = c("circle", "rectangle"),
        cex.node.name = 0.9,
        #color.edge = color.spectral(100),
        color.edge = mycols,
        keysize = c(0.83,0.83),
        #color.edge = color.GreenRed(n = 10, alpha = 0.4),
        lty.edge = "solid", lwd.edge = 5, #dotted instead of solid
        show.edge.labels = F,
        symkey = T)
dev.off()




#**********************
#* 3 - GUT d33 d39 ----
#**********************


#* Common tables ----

#metadata
metadata.d33.d39 <- metadata %>% 
  dplyr::filter((time == "d33" | time == "d39") & group == "gut")
#factor
metadata.d33.d39$pub.treat <- factor(metadata.d33.d39$pub.treat)
#rename levels
levels(metadata.d33.d39$pub.treat) <- c("d33","d39")

#ASV new names for publication
asv.otu

#color
col.gut.d33.d39


#* Tables ----

#ASV table with abundances for filtering (keep ASV if at least 1% in 1 sample)
rar_df %>%
  dplyr::select(c(ASV, metadata.d33.d39$publish)) %>%
  mutate_if(is.double, abund.col) %>% 
  filter(rowSums(.[,-1] >= 1) > 0) %>% #keep lines with at least one value up to 1 (%)
  dplyr::select(ASV) -> tmp
#ASV table filtered in counts
tmp <- rar_df %>%
  dplyr::select(c(ASV, metadata.d33.d39$publish)) %>% 
  filter(ASV %in% tmp$ASV) %>% 
  mutate(ASV = str_replace(ASV, "asv", "ASV")) #replace asv with ASV
#prepare centered log ratio transformation (CLR)
X <- t(tmp[,-1])
colnames(X) <- tmp$ASV
X <- mixOmics::logratio.transfo(as.matrix(X), logratio = 'CLR', offset = 1) #to prevent log transform later (https://mixomics.org/mixmc/mixmc-pre-processing/)
sum(which(X == 0)) #check if still zeros, if result 0 => it's ok
dim(X)
X <- as.data.frame(X[,]) #add [,] to avoid the error


#metabolites table (filter without some samples : 28A2g4; 28A3g3; 28C2g4; 28E3g4)
df.metabo <- read_delim("../4_analyses_metabo/R/data/gut_t28_t33_t39_pepmass_metabo.csv", delim = ',', col_names = T) %>% 
  rename(PEPMASS = "X1") %>% 
  dplyr::select(c("PEPMASS", metadata.d33.d39$id)) %>% #reorder
  dplyr::rename(`33Eg3` = `39Eg3`, `39Eg3` = `33Eg3`) %>% #rename with correct name because of name error on spectro computer
  dplyr::relocate(`39Eg3`, .after = `39Eg2`) %>% #move samples
  dplyr::relocate(`33Eg3`, .after = `33Eg2`) %>% #move samples
  dplyr::rename(`33Dg1` = `33Dg2`, `33Dg2` = `33Dg1`) %>% #rename with correct name because of name error on spectro computer
  dplyr::relocate(`33Dg1`, .before = `33Dg2`) %>% #move samples
  rename_at(vars(`33Dg1`:`39Eg4`), ~ str_replace_all(., setNames(metadata.d33.d39$publish, metadata.d33.d39$id))) %>% #pareto transformation
  filter(rowSums(.[,-1] > 0) > 0)
library(MetabolAnalyze)
Y <- MetabolAnalyze::scaling(df.metabo[,-1], "pareto")
rownames(Y) <- df.metabo$PEPMASS
Y <- t(Y)
Y <- as.data.frame(Y)


#* Choose the design matrice ####
## 1- full weighted design matrice

#with this design, all blocks are connected and this full weighted design offers a trade-off between
#maximizing the correlation between datasets and the separation between classes (Y)
design = matrix(c(0,0.1,1,
                  0.1,0,1,
                  1,1,0), ncol = 3)

#* Combine matrices ####

#combine all data tables in a big list
all <- list(microbiome = X,
            metabolome = Y)

#block matrices: indY indicates where the outcome Y is in the list X
block.plsda.medaka1 = mixOmics::block.plsda(X = all, Y = metadata.d33.d39$pub.treat,
                                            ncomp = 3, design = design)

#correlation first component
plot(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])
cor(block.plsda.medaka1$variates$microbiome[,1],block.plsda.medaka1$variates$metabolome[,1])


#* Network ----

#>>>> Sup-figure 6b ####
#pdf("../../4_Paper/figures/Sup-figure6b.pdf")
mycols <- c(colorRampPalette(c("#6587cd", "#FFFFFF"))(6),
            colorRampPalette(c("#FFFFFF", "#c85a38"))(6))
mycols = mycols[mycols != "#FFFFFF"]
p <- network(block.plsda.medaka1, cutoff = 0.714, 
        color.node = c(alpha("sandybrown",0.9),alpha("white",0.7)), 
        shape.node = c("circle", "rectangle"),
        cex.node.name = 0.9,
        #color.edge = color.spectral(100),
        color.edge = mycols,
        keysize = c(0.83,0.83),
        #color.edge = color.GreenRed(n = 10, alpha = 0.4),
        lty.edge = "solid", lwd.edge = 5, #dotted instead of solid
        show.edge.labels = F,
        symkey = T)
dev.off()

#export high correlated metabolites and fin occurence with cultures exometabolites (previous point 9)
metabolites.d33.d39 <- p$gR[1]
metabolites.d33.d39 <- names(metabolites.d33.d39)[-c(1:4)]
tmp <- as_tibble(metabolites.d33.d39, paste0(seq(1, 46)))
#write_csv(tmp, "../5_analyses_correlation_metabo-microbio/R/output/metabolites_networkd33d39_for_annotation.csv")



