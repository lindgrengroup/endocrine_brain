###############
###LIBRARIES###
###############

library(data.table)
library(ggplot2)
install.packages("cowplot")
library(cowplot)
library(dplyr)

###############
###FUNCTIONS###
###############

qqplot <- function(pvalues, ci=0.95, is.negLog=FALSE,
                   highlight=NULL, name="", size.title=12,
                   size.text=12, raster=TRUE, point_col=NULL) {
  N  <- length(pvalues)
  if (is.negLog) {
    observed <- sort(pvalues, decreasing=TRUE)
  } else {
    observed <- -log10(sort(pvalues))
  }
  
  df <- data.frame(
    observed <- observed,
    expected <- -log10(1:N / N),
    clower   <- -log10(qbeta(ci,     1:N, N - 1:N + 1)),
    cupper   <- -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))
  )
  if (!is.null(highlight)) {
    df$highlight <- highlight
  }
  xlabel <- expression(Expected~~-log[10](italic(p)))
  ylabel <- expression(Observed~~-log[10](italic(p)))
  
  p <- ggplot2::ggplot(df)
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(x=expected, ymin=clower,
                                             ymax=cupper), fill="gray90") +
    ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=max(df$expected),
                                       yend=max(df$expected)), color="gray10") +
    ggplot2::xlim(0, max(df$expected)) +
    ggplot2::labs(x=xlabel, y=ylabel) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size=size.title),
                   axis.text=ggplot2::element_text(size=size.text)
    )
  if (!is.null(highlight)) {
    if (!raster) {
      p <- p + ggplot2::geom_point(ggplot2::aes(x=expected, y=observed,
                                                color=highlight))
    } else {
      p <- p + ggrastr::geom_point_rast(ggplot2::aes(x=expected,
                                                     y=observed,
                                                     color=highlight))
    }
    p <- p + ggplot2::scale_color_manual(values=c("#32806E","gray10"),
                                         name=name)
  } else {
    if (!raster) {
      p <- p + ggplot2::geom_point(ggplot2::aes(x=expected, y=observed),
                                   col=point_col)
    } else {
      p <- p + ggrastr::geom_point_rast(ggplot2::aes(x=expected,
                                                     y=observed),
                                        col="gray10")
    }
  }
  p
}

##########
###DATA###
##########

system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt")
system("dx download data/brain_all/genotype_process/all/regenie_step2/filtered/cojo_results/merged/meta_table_cojo_cond_edit.csv")

HT <- fread("assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG <- fread("assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB <- fread("assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM <- fread("assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

meta_cojo_cond <- fread("meta_table_cojo_cond_edit.csv", data.table=FALSE)

##############
###ANALYSIS###
##############

HT$pheno <- "HT"
PG$pheno <- "PG"
OB$pheno <- "OB"
HT_GM$pheno <- "HT_GM"

all_results <- do.call(rbind, list(HT, PG, OB, HT_GM))
#all_results <- data.frame(all_results)

all_results$char <- NA
#unsignificant
all_results[which(all_results$LOG10P < -log10(5e-8) & all_results$CHROM %% 2 == 0), "char"] <- 0
all_results[which(all_results$LOG10P < -log10(5e-8) & all_results$CHROM %% 2 == 1), "char"] <- 1
#Hypothalamus
all_results[which(all_results$LOG10P > -log10(5e-8) & all_results$pheno == "HT"), "char"] <- 2
#Pituitary gland
all_results[which(all_results$LOG10P > -log10(5e-8) & all_results$pheno == "PG"), "char"] <- 3
#Olfactory bulbs
all_results[which(all_results$LOG10P > -log10(5e-8) & all_results$pheno == "OB"), "char"] <- 4
#Hypothalamus Grey Matter
all_results[which(all_results$LOG10P > -log10(5e-8) & all_results$pheno == "HT_GM"), "char"] <- 5

all_results$char <- as.factor(all_results$char)

all_results$ALPHA <- NA
all_results[which(all_results$LOG10P < -log10(5e-8)), "ALPHA"] <- "A"
all_results[which(all_results$LOG10P > -log10(5e-8)), "ALPHA"] <- "B"

all_results$EXTRA <- NULL

###############################
###META MANAHTTAN PARAMETERS###
###############################
d=all_results
chr="CHROM"
bp="GENPOS"
p="LOG10P"
snp="ID"
title=NULL
max.y="max"
min.y="min"
is.negLog=TRUE
color=c("#C2C2C2", "#676767")
a=0.5
genomewideline=-log10(5e-8)
colorGenomewide="#676767"
linetypeGenomewide=1
pch_highlight="char"
pch_values=c("#C2C2C2", "#676767", "#5bbda3","#9d83d6","#EEA243", "#b51a0e")
size.x.title=28
size.y.title=28
size.x.text=28
size.y.text=28
size.points=3
raster=FALSE

#########################
###PLOT META MANHATTAN###
#########################

if (!(chr %in% names(d))) stop(paste("Column", chr, "not found!"))
if (!(bp %in% names(d))) stop(paste("Column", bp, "not found!"))
if (!(p %in% names(d))) stop(paste("Column", p, "not found!"))

if (!is.numeric(d[[bp]])) stop(paste(bp, "column should be numeric."))
if (!is.numeric(d[[p]])) stop(paste(p, "column should be numeric."))
if (!is.numeric(d[[chr]])) {
  stop(paste(chr, "column should be numeric. Does your [chr] column",
             "chromsomes in chr1 etc format? Are there 'X', 'Y',",
             " 'MT', etc? If so, change them to numeric encoding."))
}

names(d)[names(d) == chr] <- "CHR"
names(d)[names(d) == bp] <- "BP"
names(d)[names(d) == p] <- "P"

if (!is.null(d[[snp]])) {
  names(d)[names(d) == snp] <- "SNP"
}
if (!is.null(d[[pch_highlight]])) {
  names(d)[names(d) == pch_highlight] <- "PCH"
}

d <- na.omit(d)

if (!is.negLog) {
  if (any(d$P < 0 | d$P >= 1)) stop ("P-values have to be in range (0,1]")
  d  <- d[order(d$CHR, d$BP),]
  message("Pvalues are converted to negative log10 pvalues")
  d$logp <- -log10(d$P)
} else {
  d <- d[order(d$CHR, d$BP),]
  message("log10(p values) are used to depict results")
  d$logp <- d$P
}

d$pos <- NA
ticks <- NULL
lastbase <- 0
numchroms <- length(unique(d$CHR))

if (numchroms == 1) {
  d$pos <- d$BP
} else {
  for (i in unique(d$CHR)) {
    if (i == 1) {
      d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP
    } else {
      lastbase <- lastbase + max(subset(d, CHR==i-1)$BP)
      d[d$CHR==i, ]$pos <- d[d$CHR==i, ]$BP + lastbase
    }
    ticks <- c(ticks,
               d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
  }
  ticklim <- c(min(d$pos),max(d$pos))
}

if (max.y == "max") {
  maxy <- ceiling(max(d$logp))
} else {
  maxy <- max.y
}
if (min.y == "min") {
  miny <- floor(min(d$logp))
} else {
  miny <- min.y
}
if (maxy < 8) {
  maxy <- 8
}

mycols <- rep(color, max(d$CHR))
ylab <-expression(-log[10](italic(p)))

if (numchroms == 1) {
  p <- ggplot2::ggplot(data=d, ggplot2::aes(x=pos, y=logp))
  if (! raster) {
    p <- p + ggplot2::geom_point(size=size.points)
  } else {
    p <- p + ggrastr::geom_point_rast(size=size.points)
  }
  p <- p + ggplot2::ylab(expression(-log[10](italic(p)))) +
    ggplot2::xlab(paste("Chromosome", unique(d$CHR),"position"))
} else {
  p <- ggplot2::ggplot(data=d, ggplot2::aes(x=pos, y=logp))
  p <- p + ggplot2::ylab(expression(-log[10](italic(p))))
  p  <- p + ggplot2::scale_x_continuous(name="Chromosome", breaks=ticks,
                                        limits=ticklim, expand=c(0.01,0.01),
                                        labels=(unique(d$CHR)))
  p <- p + ggplot2::scale_y_continuous(limits = c(miny, maxy),
                                       expand=c(0.01,0.01))
}


p <- p + ggplot2::scale_colour_manual(values=color)
p <- p + ggplot2::geom_point(ggplot2::aes(color=as.factor(PCH)),
                             size=size.points)
p <- p + ggplot2::scale_colour_manual(values=pch_values, guide=FALSE)
p <- p + ggplot2::theme(legend.position="none")
p <- p + ggplot2::geom_point(ggplot2::aes(color=as.factor(PCH), alpha=ALPHA),
                             size=size.points)
p <- p + scale_alpha_manual(values=c(1,0.5), guide=F)
p  <- p + ggplot2::scale_colour_manual(name="Significant\nAssociations", values=pch_values, limits=c("0", "1", "2", "3", "4", "5"),breaks= function(x) c("2", "3", "4", "5"), drop=TRUE,
                                       labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter"), 
                                       guide=guide_legend(override.aes = list(size=28, alpha=1), keyheight=2))
if (is.null(title)) {
  p <- p + ggplot2::theme(title=title)
}
p <- p + ggplot2::theme_classic()
p <- p + ggplot2::theme(
  axis.text.x=ggplot2::element_text(size=size.x.text, colour="grey50"),
  axis.text.y=ggplot2::element_text(size=size.y.text, colour="grey50"),
  axis.title.x=ggplot2::element_text(size=size.x.title, colour="grey50"),
  axis.title.y=ggplot2::element_text(size=size.y.title, colour="grey50"),
  legend.text=(element_text(size=28)),
  legend.title=ggplot2::element_text(size=28),
  axis.ticks=ggplot2::element_blank()
)


if (genomewideline) {
  p <- p + ggplot2::geom_segment(x=min(d$pos), xend=max(d$pos),
                                 y=genomewideline, yend=genomewideline,
                                 colour=colorGenomewide,
                                 linetype=linetypeGenomewide)
}


png("meta_manhattan_paper.png", width=2400, height=1440)
p
dev.off()

pdf("meta_manhattan_paper.pdf", width=8, height=4.8)
p
dev.off()

##########################
###MAKING AN UPSET PLOT###
##########################

cond_additional <- merge(meta_cojo_cond, d, by.x="ID", by.y="SNP")
cond_additional$point_col <- NA
cond_additional$point_col[which(cond_additional$pheno=="HT" & cond_additional$HT==0)] <- 0
cond_additional$point_col[which(cond_additional$pheno=="HT" & cond_additional$HT==1)] <- 2
cond_additional$point_col[which(cond_additional$pheno=="PG" & cond_additional$PG==0)] <- 0
cond_additional$point_col[which(cond_additional$pheno=="PG" & cond_additional$PG==1)] <- 3
cond_additional$point_col[which(cond_additional$pheno=="OB" & cond_additional$OB==0)] <- 0
cond_additional$point_col[which(cond_additional$pheno=="OB" & cond_additional$OB==1)] <- 4
cond_additional$point_col[which(cond_additional$pheno=="HT_GM" & cond_additional$HT_GM==0)] <- 0
cond_additional$point_col[which(cond_additional$pheno=="HT_GM" & cond_additional$HT_GM==1)] <- 5

pch_values=c("#e8e8e8", "#676767", "#5bbda3","#9d83d6","#EEA243", "#b51a0e")
size.points=8
tick_ends <- d %>% group_by(CHR) %>% top_n(1, BP)
tick_ends <- unique(tick_ends$pos)
tick_ends <- tick_ends[-length(tick_ends)]

cond_additional$pheno <- factor(cond_additional$pheno, levels=c("HT_GM", "OB", "PG", "HT"))

p2 <- ggplot(cond_additional, aes(x=pos, y=pheno)) +
        geom_point(aes(color=factor(point_col), fill=factor(point_col)), size=size.points) +
  scale_x_continuous(name="Chromosome", breaks=ticks,
                     limits=ticklim, expand=c(0.01,0.01),
                     labels=(unique(d$CHR)))
p2 <- p2 + ggplot2::theme(legend.position="none")
p2  <- p2 + ggplot2::scale_colour_manual(name="Significant\nAssociations", values=pch_values, limits=c("0", "1", "2", "3", "4", "5"),breaks= function(x) c("2", "3", "4", "5"), drop=TRUE)
p2 <- p2 + ggplot2::scale_y_discrete(labels=c(
  "HT" = "Hypothalamus",
  "PG" = "Pituitary\nGland",
  "OB" = "Olfactory\nBulb",
  "HT_GM" = "Hypothalamus\nGrey\nMatter"
))
p2 <- p2 + ggplot2::theme_classic()
p2 <- p2 + ggplot2::theme(legend.position="none")
p2 <- p2 + ggplot2::theme(
  axis.text.x=ggplot2::element_blank(),
  axis.text.y=ggplot2::element_text(size=size.y.text, colour="grey50"),
  axis.title.x=ggplot2::element_blank(),
  axis.title.y=ggplot2::element_blank(),
  axis.line=element_blank(),
  axis.ticks = element_blank())
p2 <- p2 + geom_vline(xintercept = tick_ends, colour="#D6D6D6")


png("upset_tmp.png", width=2400, height=500)
p2
dev.off()
system("dx upload upset_tmp.png --path plots/")

p_all <- plot_grid(p, p2, ncol=1, align="v", rel_heights = c(4,1))

png("man_upset_paper.png", width=2400, height=1800)
p_all
dev.off()
system("dx upload man_upset_paper.png --path plots/")



############################
###META QQPLOT PARAMETERS###
############################

is.negLog=TRUE
raster=FALSE
ci =0.95
name=""
size.title=28
size.text=28
size.points=3
pch_values=c("#C2C2C2", "#676767", "#5bbda3","#9d83d6","#EEA243", "#b51a0e")

##################
###PLOT META QQ###
##################

ar_small <- all_results[,c("LOG10P", "char", "ALPHA")]
N <- length(ar_small$LOG10P)
ar_small_sort <- ar_small[order(ar_small$LOG10P, decreasing=TRUE),]
ar_small_sort$expected <- -log10(1:N / N)
ar_small_sort$clower   <- -log10(qbeta(ci,     1:N, N - 1:N + 1))
ar_small_sort$cupper   <- -log10(qbeta(1 - ci, 1:N, N - 1:N + 1))

xlabel <- expression(Expected~~-log[10](italic(p)))
ylabel <- expression(Observed~~-log[10](italic(p)))
  
  p <- ggplot2::ggplot(ar_small_sort)
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(x=expected, ymin=clower,
                                             ymax=cupper), fill="gray90") +
    ggplot2::geom_segment(ggplot2::aes(x=0, y=0, xend=max(ar_small_sort$expected),
                                       yend=max(ar_small_sort$expected)), color="gray10") +
    ggplot2::xlim(0, max(ar_small_sort$expected)) +
    ggplot2::labs(x=xlabel, y=ylabel) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.title=ggplot2::element_text(size=size.title),
                   axis.text=ggplot2::element_text(size=size.text)
    )

      p <- p + ggplot2::geom_point(ggplot2::aes(x=expected, y=LOG10P, color=as.factor(char), alpha=ALPHA,
                                   size=size.points))
      p <- p + scale_alpha_manual(values=c(1,0.5), guide=F)
      p  <- p + ggplot2::scale_colour_manual(name="Significant\nAssociations", values=pch_values, limits=c("0", "1", "2", "3", "4", "5"),breaks= function(x) c("2", "3", "4", "5"), drop=TRUE,
                                             labels = c("Hypothalamus", "Pituitary\nGland", "Olfactory\nBulb", "Hypothalamus\nGrey Matter"))
      p <- p + ggplot2::theme(legend.position="none")
      
      png("meta_qq_paper.png", width=1440, height=1440)
      p
      dev.off()
      
      pdf("meta_qq_paper.pdf", width=4.8, height=4.8)
      p
      dev.off()




######################################
###MAKE INDIVIDUAL COLOURED QQPLOTS###
###################################### 
      
p_HT <- qqplot(HT$LOG10P, is.negLog=TRUE, raster=FALSE, point_col="#5bbda3")
p_PG <- qqplot(PG$LOG10P, is.negLog = TRUE, raster=FALSE, point_col = "#9d83d6")
p_OB <- qqplot(OB$LOG10P, is.negLog = TRUE, raster=FALSE, point_col = "#EEA243")
p_HT_GM <- qqplot(HT_GM$LOG10P, is.negLog = TRUE, raster=FALSE, point_col= "#b51a0e")

png("HT_all_qqplot_colour.png")
p_HT
dev.off()

png("PG_all_qqplot_colour.png")
p_PG
dev.off()

png("OB_all_qqplot_colour.png")
p_OB
dev.off()

png("HT_GM_all_qqplot_colour.png")
p_HT_GM
dev.off()

system("dx upload *colour.png --path plots/")
