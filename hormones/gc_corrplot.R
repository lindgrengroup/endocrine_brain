###############
###LIBRARIES###
###############

library(data.table)
install.packages("corrplot")
library(corrplot)
library(ggplot2)
install.packages("tidyverse")
library(tidyverse)

##########
###DATA###
##########

gc <- fread("/mnt/project/data/genetic_corr_update_231125.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

gc$z_region <- gc$region_h2/gc$region_se
gc$z_trait <- gc$trait_h2/gc$trait_se

gc <- gc %>% mutate(gc = ifelse(z_region<4 | z_trait<4, NA, gc))
gc <- gc %>% mutate(p = ifelse(gc=="NA", "NA", p))
gc_complete <- gc[complete.cases(gc),]
gc_p <- 0.05/length(gc_complete$p)
gc_complete_sig <- subset(gc_complete, p < gc_p)


gc <- subset(gc, trait %in% c("prog", "test", "oest", "lh", "fsh"))
gc_complete <- gc[complete.cases(gc),]

gc_p <- 0.05/length(gc_complete$p)
gc_complete_sig <- subset(gc_complete, p < gc_p)


gc_all <- subset(gc, group=="all")
gc_female <- subset(gc, group=="female")
gc_male <- subset(gc, group=="male")

gc_all_p <- gc_all[,c("region", "trait", "p")]
gc_female_p <- gc_female[,c("region", "trait", "p")]
gc_male_p <- gc_male[,c("region", "trait", "p")]

gc_all_gc <- gc_all[,c("region", "trait", "gc")]
gc_female_gc <- gc_female[,c("region", "trait", "gc")]
gc_male_gc <- gc_male[,c("region", "trait", "gc")]

gc_all_p_short <- spread(gc_all_p, "trait", "p")
gc_female_p_short <- spread(gc_female_p, "trait", "p")
gc_male_p_short <- spread(gc_male_p, "trait", "p")

gc_all_gc_short <- spread(gc_all_gc, "trait", "gc")
gc_female_gc_short <- spread(gc_female_gc, "trait", "gc")
gc_male_gc_short <- spread(gc_male_gc, "trait", "gc")

row.names(gc_all_gc_short) <- gc_all_gc_short$region
gc_all_gc_short$region <- NULL
gc_all_gc_short_filled <- gc_all_gc_short %>%
  add_column(infert1=NA) %>% 
  add_column(infert2=NA) %>% 
  add_column(infert3=NA) %>% 
  add_column(infert4=NA) %>% 
  add_column(infert5=NA) %>% 
  add_column(male_infert=NA)
gc_all_gc_short_filled <- gc_all_gc_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_all_gc_short_filled <- gc_all_gc_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_all_gc_short_filled <- as.matrix(gc_all_gc_short_filled)
colnames(gc_all_gc_short_filled) <- c("Testosterone", "FSH", "LH", "Oestrodiol", "Progesterone")
rownames(gc_all_gc_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")
row.names(gc_all_p_short) <- gc_all_p_short$region
gc_all_p_short$region <- NULL
gc_all_p_short_filled <- gc_all_p_short %>%
  add_column(infert1=NA) %>% 
  add_column(infert2=NA) %>% 
  add_column(infert3=NA) %>% 
  add_column(infert4=NA) %>% 
  add_column(infert5=NA) %>% 
  add_column(male_infert=NA)
gc_all_p_short_filled <- gc_all_p_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_all_p_short_filled <- gc_all_p_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_all_p_short_filled <- as.matrix(gc_all_p_short_filled)
colnames(gc_all_p_short_filled) <- c("Testosterone", "FSH", "LH", "Oestrodiol", "Progesterone")
rownames(gc_all_p_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")

row.names(gc_female_gc_short) <- gc_female_gc_short$region
gc_female_gc_short$region <- NULL
gc_female_gc_short_filled <- gc_female_gc_short %>%
  add_column(male_infert=NA)
gc_female_gc_short_filled <- gc_female_gc_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_female_gc_short_filled <- gc_female_gc_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_female_gc_short_filled <- as.matrix(gc_female_gc_short_filled)
rownames(gc_female_gc_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")
row.names(gc_female_p_short) <- gc_female_p_short$region
gc_female_p_short$region <- NULL
gc_female_p_short_filled <- gc_female_p_short %>%
  add_column(male_infert=NA)
gc_female_p_short_filled <- gc_female_p_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_female_p_short_filled <- gc_female_p_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_female_p_short_filled <- as.matrix(gc_female_p_short_filled)
rownames(gc_female_p_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")


row.names(gc_male_gc_short) <- gc_male_gc_short$region
gc_male_gc_short$region <- NULL
gc_male_gc_short_filled <- gc_male_gc_short %>%
  add_column(infert1=NA) %>% 
  add_column(infert2=NA) %>% 
  add_column(infert3=NA) %>% 
  add_column(infert4=NA) %>% 
  add_column(infert5=NA) %>% 
  add_column(prog=NA)
gc_male_gc_short_filled <- gc_male_gc_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_male_gc_short_filled <- gc_male_gc_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_male_gc_short_filled <- as.matrix(gc_male_gc_short_filled)
rownames(gc_male_gc_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")
row.names(gc_male_p_short) <- gc_male_p_short$region
gc_male_p_short$region <- NULL
gc_male_p_short_filled <- gc_male_p_short %>%
  add_column(infert1=NA) %>% 
  add_column(infert2=NA) %>% 
  add_column(infert3=NA) %>% 
  add_column(infert4=NA) %>% 
  add_column(infert5=NA) %>% 
  add_column(prog=NA)
gc_male_p_short_filled <- gc_male_p_short_filled[,c("test", "fsh", "lh", "oest", "prog")]
gc_male_p_short_filled <- gc_male_p_short_filled[c("HT", "PG", "OB", "HT_GM"),]
gc_male_p_short_filled <- as.matrix(gc_male_p_short_filled)
rownames(gc_male_p_short_filled) <- c("Hypothalamus", "Pituitary gland", "Olfactory bulb", "Hypothalamus\ngrey matter")




##############
###PLOTTING###
##############

p_all <- corrplot(gc_all_gc_short_filled, p.mat=gc_all_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio=0.2, method="color", col=COL2('PuOr'), na.label="X")
p_female <- corrplot(gc_female_gc_short_filled, p.mat=gc_female_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio = 0.2, method="color", col=COL2('PuOr'), na.label="X")
p_male <- corrplot(gc_male_gc_short_filled, p.mat=gc_male_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio = 0.2, method="color", col=COL2('PuOr'), na.label="X")

pdf("gc_corrplot_all_update260624.pdf", width=7, height = 5)
p_all <- corrplot(gc_all_gc_short_filled, p.mat=gc_all_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio=0.2, method="color", col=COL2('PuOr'), na.label="X", na.label.col = "grey51", insig = 'label_sig', pch.col = 'grey20')
dev.off()

pdf("gc_corrplot_female_update260624.pdf", width=7, height = 5)
p_female <- corrplot(gc_female_gc_short_filled, p.mat=gc_female_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio = 0.2, method="color", col=COL2('PuOr'), na.label="X", na.label.col = "grey51", insig = 'label_sig', pch.col = 'grey20')
dev.off()

pdf("gc_corrplot_male_update260624.pdf", width=7, height = 5)
p_male <- corrplot(gc_male_gc_short_filled, p.mat=gc_male_p_short_filled, sig.level=0.05, tl.col="black", cl.ratio = 0.2, method="color", col=COL2('PuOr'), na.label="X", na.label.col = "grey51", insig = 'label_sig', pch.col = 'grey20')
dev.off()


