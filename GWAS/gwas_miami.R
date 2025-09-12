###############
###LIBRARIES###
###############

library(data.table)
library(ggplot2)
library(dplyr)
library(grid)
library(gridExtra)

###############
###FUNCTIONS###
###############


miami <- function(d, h, chr = "CHR", bp = "BP", p = "P", snp="SNP",
                      title=NULL, max.y="max", min.y="min", is.negLog=FALSE,
                      highlight=NULL, colorHighlight="green",
                      color_d=c("#b51a0e", "#F5827A"), color_h=c("#b51a0e", "#721009"),
                      a=0.5, genomewideline=-log10(5e-8), colorGenomewide="gray90",
                      linetypeGenomewide=1,
                      size.x.title=12, size.y.title=12,
                      size.x.text=12, size.y.text=12,size.points=1) {
  
  if (!(chr %in% names(d))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(d))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(d))) stop(paste("Column", p, "not found!"))
  
  if (!(chr %in% names(h))) stop(paste("Column", chr, "not found!"))
  if (!(bp %in% names(h))) stop(paste("Column", bp, "not found!"))
  if (!(p %in% names(h))) stop(paste("Column", p, "not found!"))
  
  if (!is.numeric(d[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(d[[p]])) stop(paste(p, "column should be numeric."))
  if (!is.numeric(d[[chr]])) {
    stop(paste(chr, "column should be numeric. Does your [chr] column",
               "chromsomes in chr1 etc format? Are there 'X', 'Y',",
               " 'MT', etc? If so, change them to numeric encoding."))
  }
  
  if (!is.numeric(h[[bp]])) stop(paste(bp, "column should be numeric."))
  if (!is.numeric(h[[p]])) stop(paste(p, "column should be numeric."))
  if (!is.numeric(h[[chr]])) {
    stop(paste(chr, "column should be numeric. Does your [chr] column",
               "chromsomes in chr1 etc format? Are there 'X', 'Y',",
               " 'MT', etc? If so, change them to numeric encoding."))
  }
  
  names(d)[names(d) == chr] <- "CHR"
  names(d)[names(d) == bp] <- "BP"
  names(d)[names(d) == p] <- "P"
  
  names(h)[names(h) == chr] <- "CHR"
  names(h)[names(h) == bp] <- "BP"
  names(h)[names(h) == p] <- "P"
  
  if (!is.null(d[[snp]])) {
    names(d)[names(d) == snp] <- "SNP"
    names(h)[names(h) == snp] <- "SNP"
  }
  d <- na.omit(d)
  h <- na.omit(h)
  
  if (!is.negLog) {
    if (any(d$P < 0 | d$P >= 1)) stop ("P-values have to be in range (0,1]")
    if (any(h$P < 0 | h$P >= 1)) stop ("P-values have to be in range (0,1]")
    d  <- d[order(d$CHR, d$BP),]
    h  <- h[order(h$CHR, h$BP),]
    message("Pvalues are converted to negative log10 pvalues")
    d$logp <- -log10(d$P)
    h$logp <- -log10(h$P)

  } else {
    d <- d[order(d$CHR, d$BP),]
    h <- h[order(h$CHR, h$BP),]
    message("log10(p values) are used to depict results")
    d$logp <- d$P
    h$logp <- h$P
  }
  
  d <- as_tibble(d)
  h <- as_tibble(h)
  d <- d %>% rename_at(vars(-CHR, -BP, -SNP, -ALLELE0, -ALLELE1, -TEST), funs(paste0("female:", .)))
  h <- h %>% rename_at(vars(-CHR, -BP, -SNP, -ALLELE0, -ALLELE1, -TEST), funs(paste0("male:", .)))
  dh <- full_join(d, h)
  
  dh$pos <- NA
  ticks <- NULL
  lastbase <- 0
  numchroms <- length(unique(dh$CHR))
  
  if (numchroms == 1) {
    dh$pos <- dh$BP
  } else {
    for (i in unique(dh$CHR)) {
      if (i == 1) {
        dh[dh$CHR==i, ]$pos <- dh[dh$CHR==i, ]$BP
      } else {
        lastbase <- lastbase + max(subset(dh, CHR==i-1)$BP)
        dh[dh$CHR==i, ]$pos <- dh[dh$CHR==i, ]$BP + lastbase
      }
      ticks <- c(ticks,
                 dh[dh$CHR==i, ]$pos[floor(length(dh[dh$CHR==i, ]$pos)/2)+1])
    }
    ticklim <- c(min(dh$pos),max(dh$pos))
  }
  
  if (max.y == "max") {
    maxy <- ceiling(max(na.omit(c(dh$`female:logp`, dh$`male:logp`))))
  } else {
    maxy <- max.y
  }
  if (min.y == "min") {
    miny <- floor(min(na.omit(c(dh$`female:logp`, dh$`male:logp`))))
  } else {
    miny <- min.y
  }
  if (maxy < 8) {
    maxy <- 8
  }
  
  mycols_d <- rep(color_d, max(dh$CHR))
  mycols_h <- rep(color_h, max(dh$CHR))
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
    pd <- ggplot2::ggplot(data=dh, ggplot2::aes(x=pos, y=`female:logp`))
#    pd <- pd + ggplot2::ylab(expression(-log[10](italic(p))))
    pd <- pd + ggplot2::ylab("")
    pd  <- pd + ggplot2::scale_x_continuous(name="Chromosome", breaks=ticks,
                                          limits=ticklim, expand=c(0.01,0.01),
                                          labels=(unique(dh$CHR)))
    pd <- pd + ggplot2::scale_y_continuous(limits = c(miny, maxy),
                                         expand=c(0.01,0.01))
  }
  if (genomewideline) {
    pd <- pd + ggplot2::geom_segment(x=min(dh$pos), xend=max(dh$pos),
                                     y=genomewideline, yend=genomewideline,
                                     colour=colorGenomewide,
                                     linetype=linetypeGenomewide)
  }

    pd <- pd + ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)),
                                   size=size.points)
   
    pd <- pd + ggplot2::scale_colour_manual(values=mycols_d, guide=FALSE)
    pd <- pd + ggplot2::theme(legend.position="none")
    
    
    #MALE
    ph <- ggplot2::ggplot(data=dh, ggplot2::aes(x=pos, y=`male:logp`))
    #    pd <- pd + ggplot2::ylab(expression(-log[10](italic(p))))
    ph <- ph + ggplot2::ylab("")
    ph  <- ph + ggplot2::scale_x_continuous(breaks=ticks, position="top",
                                            limits=ticklim, expand=c(0.01,0.01),
                                            labels=(unique(dh$CHR))) + labs(x=NULL)
    if (genomewideline) {
      ph <- ph + ggplot2::geom_segment(x=min(dh$pos), xend=max(dh$pos),
                                       y=-genomewideline, yend=-genomewideline,
                                       colour=colorGenomewide,
                                       linetype=linetypeGenomewide)
    }
    ph <- ph + ggplot2::scale_y_reverse(limits = c(maxy, miny),
                                          expand=c(0.01,0.01))


  ph <- ph + ggplot2::geom_point(ggplot2::aes(color=as.factor(CHR)),
                               size=size.points)


  ph <- ph + ggplot2::scale_colour_manual(values=mycols_h, guide=FALSE)
  ph <- ph + ggplot2::theme(legend.position="none")

    
  if (!is.null(highlight)) {
    if (any(!(highlight %in% as.vector(d$SNP)))) {
      warning("SNPs selected for highlighting do not exist in d")
    }
    if (any(!(highlight %in% as.vector(h$SNP)))) {
      warning("SNPs selected for highlighting do not exist in h")
    }
    dh.annotate <- dh[dh$SNP %in% highlight, ]
    pd <- pd + ggplot2::geom_point(data=dh.annotate, colour=I(colorHighlight),
                                 size=size.points)
    ph <- ph + ggplot2::geom_point(data=dh.annotate, colour=I(colorHighlight),
                                 size=size.points)
  }
  
  pd <- pd + ggplot2::theme_classic()
  pd <- pd + ggplot2::theme(
    axis.text.x=ggplot2::element_text(size=size.x.text, colour="grey50"),
    axis.text.y=ggplot2::element_text(size=size.y.text, colour="grey50"),
    axis.title.x=ggplot2::element_text(size=size.x.title, colour="grey50"),
    axis.title.y=ggplot2::element_text(size=size.y.title, colour="grey50"),
    axis.ticks=ggplot2::element_blank()
  )
  ph <- ph + ggplot2::theme_classic()
  ph <- ph + ggplot2::theme(
    axis.text.x=ggplot2::element_text(size=size.x.text, colour="grey50"),
    axis.text.y=ggplot2::element_text(size=size.y.text, colour="grey50"),
    axis.title.x=ggplot2::element_text(size=size.x.title, colour="grey50"),
    axis.title.y=ggplot2::element_text(size=size.y.title, colour="grey50"),
    axis.ticks=ggplot2::element_blank()
  )
  
#  if (genomewideline) {
#    pd <- pd + ggplot2::geom_segment(x=min(dh$pos), xend=max(dh$pos),
#                                   y=genomewideline, yend=genomewideline,
#                                   colour=colorGenomewide,
#                                   linetype=linetypeGenomewide)
#  }
#  if (genomewideline) {
#    ph <- ph + ggplot2::geom_segment(x=min(dh$pos), xend=max(dh$pos),
#                                   y=genomewideline, yend=genomewideline,
#                                   colour=colorGenomewide,
#                                   linetype=linetypeGenomewide)
#  }
  if (is.null(title)) {
    p <- p + ggplot2::theme(title=title)
  }
  p <- grid.arrange(arrangeGrob(pd + theme(legend.position="none"), 
                           ph + theme(legend.position="none"),
                           nrow = 2,
                           left = textGrob(expression(-log[10](italic(p))), rot = 90, vjust = 1)))
  if (is.null(title)) {
    p <- p + ggplot2::theme(title=title)
  }
  p
}



########################
###PARSING PARAMETERS###
########################

dir.create("female_only")
dir.create("male_only")

system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")
system("dx download data/brain_all/genotype_process/female_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o female_only/")

system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")
system("dx download data/brain_all/genotype_process/male_only/regenie_step2/filtered/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt -o male_only/")

##########
###DATA###
##########

HT_female <- fread("female_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_female <- fread("female_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_female <- fread("female_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_female <- fread("female_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

HT_male <- fread("male_only/assoc.resid_HT.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
PG_male <- fread("male_only/assoc.resid_PG.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
OB_male <- fread("male_only/assoc.resid_LR.volume.threshold.0.3.warpResolution.2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)
HT_GM_male <- fread("male_only/assoc.resid_HT-SumGM-threshold=0.3-warpResolution=2mm_norm.regenie.merged.withX.filtered.txt", data.table=FALSE)

##############
###ANALYSIS###
##############

HT_female$CHROM <- as.numeric(HT_female$CHROM)
HT_female$GENPOS <- as.numeric(HT_female$GENPOS)
HT_female$EXTRA <- NULL
PG_female$CHROM <- as.numeric(PG_female$CHROM)
PG_female$GENPOS <- as.numeric(PG_female$GENPOS)
PG_female$EXTRA <- NULL
OB_female$CHROM <- as.numeric(OB_female$CHROM)
OB_female$GENPOS <- as.numeric(OB_female$GENPOS)
OB_female$EXTRA <- NULL
HT_GM_female$CHROM <- as.numeric(HT_GM_female$CHROM)
HT_GM_female$GENPOS <- as.numeric(HT_GM_female$GENPOS)
HT_GM_female$EXTRA <- NULL

HT_male$CHROM <- as.numeric(HT_male$CHROM)
HT_male$GENPOS <- as.numeric(HT_male$GENPOS)
HT_male$EXTRA <- NULL
PG_male$CHROM <- as.numeric(PG_male$CHROM)
PG_male$GENPOS <- as.numeric(PG_male$GENPOS)
PG_male$EXTRA <- NULL
OB_male$CHROM <- as.numeric(OB_male$CHROM)
OB_male$GENPOS <- as.numeric(OB_male$GENPOS)
OB_male$EXTRA <- NULL
HT_GM_male$CHROM <- as.numeric(HT_GM_male$CHROM)
HT_GM_male$GENPOS <- as.numeric(HT_GM_male$GENPOS)
HT_GM_male$EXTRA <- NULL

#Use #5bbda3 and #B6E2D6 for females, #5bbda3 and #235849 for males
png("HT_miami_paper.png", width=800, height=550)
miami(HT_female, HT_male, chr="CHROM", bp = "GENPOS", p = "LOG10P", snp="ID", is.negLog=TRUE)
dev.off()

#Use #9d83d6 and #DBD1F0 for females, #9d83d6 and ##533399 for males
png("PG_miami_paper.png", width=800, height=550)
miami(PG_female, PG_male, chr="CHROM", bp = "GENPOS", p = "LOG10P", snp="ID", is.negLog=TRUE)
dev.off()

#Use #EEA243 and #F6D0A2 for females, #EEA243 and #BB6E11 for males
png("OB_miami_paper.png", width=800, height=550)
miami(OB_female, OB_male, chr="CHROM", bp = "GENPOS", p = "LOG10P", snp="ID", is.negLog=TRUE)
dev.off()

#Use #b51a0e and #F5827A for females, #b51a0e and #721009 for males
png("HT_GM_miami_paper.png", width=800, height=550)
miami(HT_GM_female, HT_GM_male, chr="CHROM", bp = "GENPOS", p = "LOG10P", snp="ID", is.negLog=TRUE)
dev.off()

system("dx upload *miami_paper.png --path plots/")
