###############
###LIBRARIES###
###############

library(data.table)
library(dplyr)
library(tidyr)
library(ggplot2)

###############
###FUNCTIONS###
###############

'%!in%' <- function(x,y)!('%in%'(x,y))

##########
###DATA###
##########

sex_age <- fread("/mnt/project/Phenotypes_MRI/epi_age_sex_participant.csv", data.table=FALSE)
bridge <- fread("/mnt/project/bridge_11867_8107.csv", data.table=FALSE)
pheno_file <- fread("/mnt/project/Phenotypes_MRI/allIDPs_v5_2023-07-18_norm.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

sex_age$eid <- as.character(sex_age$eid)
bridge$eid_11867 <- as.character(bridge$eid_11867)
bridge$eid_8107 <- as.character(bridge$eid_8107)
pheno_file$subjectID <- as.character(pheno_file$subjectID)
pheno_file$subjectID <- substring(pheno_file$subjectID, 2)



volumes <- pheno_file[,c("subjectID", "OB_left-volume-threshold=0.3-warpResolution=2mm", "OB_right-volume-threshold=0.3-warpResolution=2mm",
                             "HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm")]


volumes$`OB_LR-volume-threshold=0.3-warpResolution=2mm` <- rowMeans(volumes[,c("OB_left-volume-threshold=0.3-warpResolution=2mm","OB_right-volume-threshold=0.3-warpResolution=2mm")], na.rm=TRUE)

volumes$`OB_left-volume-threshold=0.3-warpResolution=2mm` <- NULL
volumes$`OB_right-volume-threshold=0.3-warpResolution=2mm` <- NULL

volumes_bridge <- merge(volumes, bridge, by.x="subjectID", by.y="eid_8107")

volume_age_sex <- merge(volumes_bridge, sex_age, by.x="eid_11867", by.y="eid")
volume_age_sex_clean <- subset(volume_age_sex, p22001 %in% c("Female", "Male"))

volume_age_sex_long <- volume_age_sex_clean %>% gather("measure", "value", -c(subjectID,eid_11867,p22001, p21003_i2))
volume_age_sex_long$p21003_i2 <- as.numeric(volume_age_sex_long$p21003_i2)
volume_age_sex_long <- volume_age_sex_long %>% mutate(age_bin=cut(p21003_i2, breaks=c(45,50,55,60,65,70,75,80,85)))
volume_age_sex_long_clean <- volume_age_sex_long[complete.cases(volume_age_sex_long),]
volume_age_sex_long_clean$measure <- as.factor(volume_age_sex_long_clean$measure)
volume_age_sex_long_clean_order <- arrange(transform(volume_age_sex_long_clean, measure=factor(measure, 
      levels=c("HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm",
               "OB_LR-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm"))), measure)

###GET VALUES FOR SUMMARY TABLE###
volume_age_sex_clean$p21003_i2 <- as.numeric(volume_age_sex_clean$p21003_i2)
volume_age_sex_clean<- volume_age_sex_clean %>% mutate(age_bin=cut(p21003_i2, breaks=c(45,55,65,75,85)))

#age_sex_mean <- aggregate(volume_age_sex_clean[,c("HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm",
#                                         "OB_LR-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm")],
#                 list(volume_age_sex_clean$p22001, volume_age_sex_clean$age_bin), FUN=function(x) c(mean=mean(x), sd=sd(x)))
age_mean <- aggregate(volume_age_sex_clean[,c("HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm",
                                              "OB_LR-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm")],
                      list(volume_age_sex_clean$age_bin), FUN=function(x) c(mean=mean(x), sd=sd(x)))
sex_mean <- aggregate(volume_age_sex_clean[,c("HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm",
                                              "OB_LR-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm")],
                      list(volume_age_sex_clean$p22001), FUN=function(x) c(mean=mean(x), sd=sd(x)))
#all_mean <- summarise(volume_age_sex_clean[,c("HT-volume-threshold=0.3-warpResolution=2mm", "PG-volume-threshold=0.3-warpResolution=2mm",
 #                                             "OB_LR-volume-threshold=0.3-warpResolution=2mm", "HT-SumGM-threshold=0.3-warpResolution=2mm")],
#                     mean=mean(), sd=sd())

volume_age_sex_clean_long <- volume_age_sex_clean %>% pivot_longer(cols=-c(eid_11867, subjectID, p22001, p21003_i2, age_bin))
volume_age_sex_clean_long_sum <- volume_age_sex_clean_long %>% group_by(name) %>% summarise(mean=mean(value), sd=sd(value))
write.csv(volume_age_sex_clean_long_sum, "epi_brain_age_sex_mean_all.csv", quote=FALSE, row.names=FALSE)
system("dx upload epi_brain_age_sex_mean_all.csv --path data/")

mean_vals <- bind_rows(age_mean, sex_mean, all_mean)
write.csv(mean_vals, "epi_brain_age_sex_mean.csv", quote=FALSE, row.names=FALSE)
system("dx upload epi_brain_age_sex_mean.csv --path data/")

##############
###PLOTTING###
##############

volume_names <- list(
  "HT-volume-threshold=0.3-warpResolution=2mm"="Hypothalamus volume",
  "PG-volume-threshold=0.3-warpResolution=2mm"="Pituitary gland volume",
  "OB_LR-volume-threshold=0.3-warpResolution=2mm"="Olfactory bulb volume",
  "HT-SumGM-threshold=0.3-warpResolution=2mm"="Hypothalamus grey matter volume"
)

volume_labeller <- function(variable,value){
  return(volume_names[value])
}

age_names <- list(
  "(45,50]"="45-50", "(50,55]"="50-55", "(55-60]"="55-60",
  "(60-65]"="60-65", "(65,70]"="65-70", "(70-75]"="70-75",
  "(75-80]"="75-80", "(80,85]"="80-85"
)
age_labeller <- function(variable,value){
  return(age_names[value])
}

 p <- ggplot(volume_age_sex_long_clean_order, aes(x=age_bin, y=value, color=p22001)) + geom_boxplot() + 
  facet_wrap(~measure, scales="free", labeller=volume_labeller) + labs(x="Age", y="Volume", colour="Sex") +
  scale_y_continuous(labels = scales::comma) + scale_x_discrete(labels=age_labeller) + scale_colour_manual(values=c("#D86E81", "#2708A0")) +
  theme_minimal() +theme(axis.text.x=element_text(angle=45))

pdf("epi_brain_age_sex.pdf", width=8, height=4.8)
p
dev.off()









