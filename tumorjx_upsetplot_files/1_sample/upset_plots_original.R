#!/usr/bin/env/R

# Generate upset plots from data
# Code author: SKM
# Contact: maden@ohsu.edu

#-------------------
# Load data objects
#-------------------
load("jx-datlist-mainfig_final.rda")
load("jx-datlist-suppfig_final.rda")

#------------------
# Main figure plot
#------------------
lro <- lro.main

# filter data for display
wholeset <- lro$sizes_data_scaled$set_sizes
subset <- lro$sizes_data_scaled$subset_sizes
wholeset <- wholeset[,!grepl("gencode|paired|gtex|other", colnames(wholeset))]
subset <- subset[,!grepl("gencode|paired|gtex|other", colnames(subset))]
wholeset <- wholeset[,rev(order(apply(wholeset,2,mean)))]
colnames(wholeset)

# make ggplot2 graph objects
# note: setsize.grplabels should match colnames of wholeset
# make the axis labels
setseq <- rev(seq(-1, 2, 1))
setsizeplot.ylab <- c()
for(s in 1:length(setseq)){
  ns <- paste0(eval(parse(text=paste0("1e",setseq[s]))),"%")
  setsizeplot.ylab <- c(setsizeplot.ylab, ns)
}
subsetseq <- rev(seq(-3, 2, 1))
subsetplot.ylab <- c()
for(s in 1:length(subsetseq)){
  ns <- paste0(eval(parse(text=paste0("format(1e",subsetseq[s],",scientific=FALSE)"))),"%")
  subsetplot.ylab <- c(subsetplot.ylab, ns)
}
lr <- make.ggset.final(setsizes.data = wholeset, subsetsizes.data = subset, 
                       setsize.grplabels = c("Unexplained", "Other Adult", "Developmental", "Stem Cells", "Placenta"),
                       hl.name = "unexplained", ugrp = "", 
                       setsizeplot.dims = unit(c(1.4, -2.6, 0.005, 0.1), 'cm'),
                       setsizeplot.ybreaks = setseq,
                       setsizeplot.ylab = setsizeplot.ylab,
                       subsetplot.dims = unit(c(0.1, 2.8, -1.65, 0.7), 'cm'),
                       subsetplot.ybreaks = subsetseq,
                       subsetplot.ylab =  subsetplot.ylab,
                       upsetplot.dims = unit(c(1.5, 1, 1.1, 1), 'cm'),
                       subtextsize = 15, sstextsize = 15, ssxaxistextsize = 12)
save(lr, file="gglist_mainfig.rda")

plotname = "upset_mainfig_inctxt.pdf"
pdf(plotname, 10, 6.5)
grid.arrange(blank,
             lr$sets_plot + labs(y = "Mean Size (% all junctions)"), 
             lr$subsets_plot + labs(y = "Intersect (% all junctions)"), 
             lr$upset_plot, 
             layout_matrix=matrix(c(1,2,rep(c(3,4),3)),nrow=2))
dev.off()

# store the new plot
plotname = "upset_mainfig.jpg"
blank <- ggplot() + theme_void()
jpeg(plotname,10,6.5,units="in",res=400)
grid.arrange(blank,
             lr$sets_plot + labs(y = "Mean Size (% all junctions)"), 
             lr$subsets_plot + labs(y = "Intersect (% all junctions)"), 
             lr$upset_plot, 
             layout_matrix=matrix(c(1,2,rep(c(3,4),3)),nrow=2))
dev.off()

#------------------------
# Supplement figure plot
#------------------------
lro <- lro.supp
# filter data for display
wholeset <- lro$sizes_data_scaled$set_sizes
subset <- lro$sizes_data_scaled$subset_sizes
set.exclude <- c(colnames(wholeset)[!colnames(wholeset) %in% unique(unlist(strsplit(colnames(subset),";")))])
wholeset <- wholeset[,!grepl(paste0(c("gencode","paired","gtex",set.exclude),collapse="|"), colnames(wholeset))]
subset <- subset[,!grepl(paste0(c("gencode","paired","gtex",set.exclude),collapse="|"), colnames(subset))]
wholeset <- wholeset[,order(apply(wholeset,2,mean))]
colnames(wholeset)

# make ggplot2 graph objects
# make the axis labels
{
  setseq <- rev(seq(-1, 2, 1))
  setsizeplot.ylab <- c()
  for(s in 1:length(setseq)){
    ns <- paste0(eval(parse(text=paste0("1e",setseq[s]))),"%")
    setsizeplot.ylab <- c(setsizeplot.ylab, ns)
  }
  subsetseq <- rev(seq(-4, 2, 1))
  subsetplot.ylab <- c()
  for(s in 1:length(subsetseq)){
    ns <- paste0(eval(parse(text=paste0("format(1e",subsetseq[s],",scientific=FALSE)"))),"%")
    subsetplot.ylab <- c(subsetplot.ylab, ns)
  }
  
  cp <- make.ggset.final(setsizes.data = wholeset, 
                         subsetsizes.data = subset,
                         hl.name = "unexplained",
                         setsize.grplabels = c("Oocyte", "Ectoderm", "Placenta", "Zygote",
                                               "Myoblast", "Neonate Fetal", "Stem Cells", "Other Embryo",
                                               "Developmental", "Other Adult", "Unexplained"),
                         ugrp="", use.nsubfilt=TRUE, nsub.filt = 40,
                         setsizeplot.dims = unit(c(0.9, -2.62, 0.35, 0.3), 'cm'), 
                         subsetplot.dims = unit(c(0.2, 1.1, -1.2, 0.85), 'cm'),
                         upsetplot.dims = unit(c(1, 0.5, 1.2, 2), 'cm'),
                         setsizeplot.ybreaks = setseq, 
                         setsizeplot.ylab = setsizeplot.ylab,
                         subsetplot.ybreaks = subsetseq,
                         subsetplot.ylab = subsetplot.ylab)
  
  save(cp, file="gglist_suppfig.rda")
  
  # store the new plot
  lr <- cp
  plotname = "upset_supplfig.jpg"
  blank <- ggplot() + theme_void()
  jpeg(plotname,17,8,units="in",res=600)
  grid.arrange(blank,
               lr$sets_plot + labs(y = "Mean Size (% all junctions)"), 
               lr$subsets_plot + labs(y = "Intersect (% all junctions)"), 
               lr$upset_plot, 
               layout_matrix=matrix(c(1,2,rep(c(3,4),3)),nrow=2))
  dev.off()
  
}
