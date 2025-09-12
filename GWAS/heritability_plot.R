###############
###LIBRARIES###
###############
install.packages("readxl")
library(data.table)
library(ggplot2)
library(reporter)
library(dplyr)
library("readxl")
library(tidyr)
install.packages("ggsignif")
library(ggsignif)

###############
###FUNCTIONS###
###############

t.test2 <- function(m1,m2,s1,s2,n1,n2,m0=0,equal.variance=FALSE)
{
  if( equal.variance==FALSE ) 
  {
    se <- sqrt( (s1^2/n1) + (s2^2/n2) )
    # welch-satterthwaite df
    df <- ( (s1^2/n1 + s2^2/n2)^2 )/( (s1^2/n1)^2/(n1-1) + (s2^2/n2)^2/(n2-1) )
  } else
  {
    # pooled standard deviation, scaled by the sample sizes
    se <- sqrt( (1/n1 + 1/n2) * ((n1-1)*s1^2 + (n2-1)*s2^2)/(n1+n2-2) ) 
    df <- n1+n2-2
  }      
  t <- (m1-m2-m0)/se 
  dat <- c(m1-m2, se, t, 2*pt(-abs(t),df))    
  names(dat) <- c("Difference of means", "Std Error", "t", "p-value")
  return(dat) 
}

##########
###DATA###
##########

h2 <- fread("/mnt/project/data/h2_sex_update_231125.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

#Add z-score for genetic correlation threshold later
h2$z <- h2$h2/h2$se

#calculate sd from se as sd=se*sqrt(n)
h2$sd <- h2$se*(sqrt(h2$n))

#subset to m and f for h2 stat
h2_mf <- subset(h2, group=="male" | group=="female")


t_OB <- (0.1554 - 0.0694)/(sqrt(((32.18665^2)/1174461) + ((31.64519^2)/1174491)))
df_OB <- 1174461 + 1174491 - 2
p_OB <- pt(t_OB, df_OB)

t_HT <- t.test2(h2$h2[which(h2$group=="female" & h2$region=="HT")], h2$h2[which(h2$group=="male" & h2$region=="HT")],
                h2$sd[which(h2$group=="female" & h2$region=="HT")], h2$sd[which(h2$group=="male" & h2$region=="HT")], 
                h2$n[which(h2$group=="female" & h2$region=="HT")], h2$n[which(h2$group=="male" & h2$region=="HT")],
                equal.variance=TRUE)

t_PG <- t.test2(h2$h2[which(h2$group=="female" & h2$region=="PG")], h2$h2[which(h2$group=="male" & h2$region=="PG")],
               h2$sd[which(h2$group=="female" & h2$region=="PG")], h2$sd[which(h2$group=="male" & h2$region=="PG")], 
               h2$n[which(h2$group=="female" & h2$region=="PG")], h2$n[which(h2$group=="male" & h2$region=="PG")],
               equal.variance=TRUE)

t_OB1 <- t.test2(h2$h2[which(h2$group=="female" & h2$region=="OB")], h2$h2[which(h2$group=="male" & h2$region=="OB")],
                 h2$sd[which(h2$group=="female" & h2$region=="OB")], h2$sd[which(h2$group=="male" & h2$region=="OB")], 
                 h2$n[which(h2$group=="female" & h2$region=="OB")], h2$n[which(h2$group=="male" & h2$region=="OB")],
                 equal.variance=TRUE)

t_HT_GM <- t.test2(h2$h2[which(h2$group=="female" & h2$region=="HT_GM")], h2$h2[which(h2$group=="male" & h2$region=="HT_GM")],
                   h2$sd[which(h2$group=="female" & h2$region=="HT_GM")], h2$sd[which(h2$group=="male" & h2$region=="HT_GM")], 
                   h2$n[which(h2$group=="female" & h2$region=="HT_GM")], h2$n[which(h2$group=="male" & h2$region=="HT_GM")],
                   equal.variance=TRUE)


###########
###PLOTS###
###########

p <- ggplot(data=h2, aes(x=region, y=h2, fill=region, alpha=group)) + geom_bar(stat="identity", position=position_dodge()) +
  geom_errorbar(aes(ymin=h2-se, ymax=h2+se), position=position_dodge(.9), width=.2) +
  scale_fill_manual(values=c("#5bbda3","#9d83d6","#EEA243", "#b51a0e"), limits=c("HT", "PG", "OB", "HT_GM"), labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter"))+
  scale_x_discrete(labels=c("HT" = "Hypothalamus", "PG" = "Pituitary\nGland",
                              "OB" = "Olfactory\nBulb", "HT_GM"="Hypothalamus\nGrey Matter"), limits=c("HT","PG","OB", "HT_GM")) +
  scale_alpha_discrete(labels=c("all"="Sex\ncombined", "female"="Female\n", "male"="Male\n"), limits=c("all","female","male"), range=c(1,0.3)) +
  theme_minimal(base_size=16) +
  labs(x="Phenotype", y=bquote('Heritability' ~(h^2)), alpha="Population") +
  guides(fill="none")
  


############
###OUTPUT###
############

pdf("sex_diff_h2_paper_231125.pdf", width=8, height=6)
p
dev.off()
system("dx upload sex_diff_h2_paper_231125.pdf --path plots/")

write.csv(h2, "h2_sex_update_231125_zscore.csv", quote=FALSE, row.names=FALSE)
system("dx upload h2_sex_update_231125_zscore.csv --path data/")
