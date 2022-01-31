

#***************
# Libraries ----
#***************

#install.packages("MuMIn")
#install.packages("lmerTest")
require(MuMIn)
require(lmerTest)
library(vegan)
library(tidyverse)




#**********
# DATA ----
#**********

setwd("~/Documents/THESE/3_Projets/1_Medaka_28j/3_Analyses/11_LMMs/")

#metadata
metadata <- read_csv(file = "./data/metadata.csv") %>%
  filter(id != "28A3g3" & id != "28C2g4" & id != "food" & id != "33Dw") %>% #remove samples under rarefaction threshold
  dplyr::filter(time == "d28" & group == "gut")
metadata$condition <- fct_relevel(metadata$condition, "0", "0+Z8", "1", "10", "100")

#asv table : rarefied at 6828 reads (on QIIME2)
rar_df <- read_tsv(file = "./data/rarefied-table.txt", skip = 1) %>% 
  rename(OTUID = "#OTU ID") %>% #rename
  dplyr::select(c("OTUID", metadata$id)) %>% #reorder
  dplyr::select("28A1g1":"28E3g3") %>% t(.)

#richness
richness <- specnumber(rar_df)

#data table for lmm
df <- tibble::data_frame(varY = richness, treatment = metadata$condition, aquarium = metadata$aquarium,
                         sample = metadata$id) %>% 
  select(varY, treatment, aquarium, sample)




#**************
# FUNCTION ----
#**************


fit_lmm <- function(input_data, model_list, nm_var_Y = NULL, reml = TRUE) {
  
  res <- list()
  
  if (! is.null(nm_var_Y)) {
    res$mods <- map(model_list, 
                    ~ lmer(formula(paste0(nm_var_Y , .x)), data = input_data,
                           REML = reml)
    )
  } else {
    res$mods <- map(model_list, 
                    ~ lmer(formula(.x), data = input_data, REML = reml)
    )
  }
  
  for (i in names(res$mods)) {
    eval(parse(text = paste0(i, " <- res$mods[[i]]")))
  }
  
  res$summaries <- lapply(res$mods, summary)
  res$anova <- lapply(res$mods, anova)
  res$rdm_anova <- lapply(res$mods, function(m) ranova(m))
  res$mod_comp <- eval(parse(text = paste0("anova(",
                                           paste0(names(res$mods), 
                                                  collapse = ","), ")")))
  res$aic_r2 <- cbind(do.call(rbind, lapply(res$mods, extractAIC)),
                      data.frame(do.call(rbind, lapply(res$mods, r.squaredGLMM)),
                                 row.names = names(res$mods))) %>% 
    setNames(c("df", "AIC","r2m","r2c"))
  
  res$mod_comp <- res$mod_comp[names(model_list), ]
  res$aic_r2 <- res$aic_r2[names(model_list), ]
  
  return(res)
}



#*********
# LMM ----
#*********

#https://stats.stackexchange.com/questions/13166/rs-lmer-cheat-sheet/61466#61466

# on the right of | = indicates nested observations within treatment
model_list <- list(mod1 = " ~ treatment + (1 | treatment : aquarium)",
                   mod2 = " ~ treatment + (1 + treatment | aquarium)")

#model_list <- list(mod1 = " ~ treatment + (1 | treatment : aquarium)",
#                   mod2 = " ~ treatment + (0 + treatment | treatment : aquarium)")


fit_lmm(input_data = df, 
        model_list = model_list,
        nm_var_Y = "varY",
        reml = TRUE)






