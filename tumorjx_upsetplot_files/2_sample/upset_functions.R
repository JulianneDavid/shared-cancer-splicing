#!/usr/bin/env/R

# Upset plot utility functions
# Code author: SKM
# Contact: maden@ohsu.edu

#-------------------
# Define Functions
#-------------------
# 1. get.datlist: function to read in cancer jx datasets
get.datlist <- function(fdpath = "August6_19_finalsetmemberships", 
                        save.dflist = TRUE, fl = list.files(fdpath),
                        dflist.name = "listdf_allcx_jxset-augfinal.rda"){
  # get.datlist
  # Description: Generate list of data frames for set analysis.
  # Arguments:
  #   fdpath (String): Path to data table csv's to read
  #   save.dflist (Boolean): Whether to save the final data list object
  #   fl (list): List of files to read from target dir
  #   dflist.name (String): Filename of data list object to save
  # Returns:
  #   cxdl: List of data frame tables for sets analysis.
  require(data.table)
  cxdl <- list()
  message("reading in data tables...")
  for(i in 1:length(fl)){
    fni = fl[i]
    ti = as.data.frame(fread(paste0(fdpath,"/", fni), sep=',', header = T))
    cxtype = gsub("_pie.*","",fni)
    cxdl[[cxtype]] <- ti
    message("finished with table (index): ",i)
  }
  message("finished reading in data tables.")
  
  if(save.dflist & !dflist.name == ""){
    message("saving the df list...")
    save(cxdl, file = dflist.name)
  }
  
  return(cxdl)
}

# 2. jx.calcsets: function to tabulate whole-set and subset tallies across junctions
jx.calcsets <- function(cxdl, cxfn = names(cxdl), read.from.disk = F, readdn = "cxdat", scale = "log10", hl.name = "novel", denom.grp="other",
                        dtall = c("gencode","paired","gtex","sra_stemcells","sra_developmental","sra_adult"),
                        ugrp = c("gtex", "paired"),
                        save.lro = TRUE, lro.name = "jx-datlist.rda"){
  # jx.calcsets
  # Description: Calculates set and overlap set sizes from cancer junction identity tables.
  # Arguments:
  #   cxdl (list) : Either of: 1. List of tables from which to calculate sets and overlaps; 2. List of table 
  #     names to read from disk
  #   read.from.disk (bool) : Whether to read in tables from disk. If T, uses cxdl as a list of filenames, 
  #     and 'readdn' as dirname containing tables.
  #   readdn (char) : Dirname containing tables to read from disk (only used if read.from.disk==T)
  #   scale (string): Specify scale for cross-table count abundances (currently only supports "log10").
  #   hl.name (string): Name of group to highlight in plots.
  #   denom.grp (string): Name of group to use as denominator in fraction calculations (either 'all' or 'other').
  #   dtall (vector of strings): Names of set categories to target (should corresond to colnames in tables).
  #   ugrp (vector of strings): Names of sets (column names in cxdl tables) to exclude/pre-filter before calculating 
  #     overlaps and subsets.
  #   save.lro (boolean): Whether to save the returned list of data tables.
  #   lro.name (string): Name of new file for returned list of data tables.
  # Returns:
  # lro : A list of data tables summarizing set sizes and subsets (scaled and unscaled).
  library(data.table)
  
  # get subsets of interest for union/filter
  setlist <- list()
  for(u in 1:length(ugrp)){
    setlist[[length(setlist)+1]] <- ugrp[u]
  }
  message("Calculting observable subsets from data...")
  dtfilt <- dtall[!dtall %in% ugrp] # filtered set names
  message("Generated dtfilt columns: ", paste0(dtfilt, collapse=" "))
  str.grp <- c()
  # populate str.grp
  # this tracks all possible identity sets across jxs
  for(c in 1:length(cxfn)){
    message("Working on table: ", cxfn[c],"...")
    if(read.from.disk==T){
      message("Reading table from disk, using cxdl names...")
      fpathc = paste0(readdn, "/", cxfn[c], ".csv")
      mx <- fread(fpathc, sep = ',', data.table = F, header = T)
    } else{
      message("Loading table from cxdl...")
      mx <- as.matrix(cxdl[[c]])
    }
    
    message("Identifying set indices in cx table...")
    mx.set <- mx[, colnames(mx) %in% dtall] # cx set members indices
    mx.set <- mx.set[, order(match(colnames(mx.set), c(dtall)))]
    # exclude union/filter grps
    message("Excluding jxs on union/filter group membership...")
    for(i in 1:length(ugrp)){
      mx.setfilt1 <- mx.set[!mx.set[, ugrp[i]]==1,]
    }
    message("Filtering data table on sets of interest in 'dtfilt'...")
    mx.setfilt2 <- mx.setfilt1[, colnames(mx.setfilt1) %in% dtfilt] # cx data with all sets of interest 
    # get filtered set identity binary codes
    message("Making filtered set id codes...")
    str.grp <- unique(c(str.grp, unique(apply(mx.setfilt2, 1, paste0, collapse=";"))))
    message("Finished with table (index): ",c)
    save(str.grp, file="str_grp.rda")
  }
  # get subset names as 'setlist'
  message("Defining sets and subsets...")
  for(s in 1:length(str.grp)){
    ss <- str.grp[s]
    ssid <- as.numeric(unlist(strsplit(ss,";")))
    setlist[[length(setlist)+1]] <- dtfilt[ifelse(ssid==1, TRUE, FALSE)]
  }
  message("Identifying 'hl.name' jx group...")
  setlist[[which(setlist == "character(0)")]] <- hl.name
  message("Calculating data tables...")
  
  #------------------------------ 
  # Calculate set and subset data
  #------------------------------
  # dim of subsets data:
  subset.cnames <- c("cx", ugrp, unlist(lapply(setlist, paste, collapse=";")),
                     "other", "all")
  ncol.subset <- length(subset.cnames) # n subsets plus 3 additional cols
  bpdf <- matrix(nrow = 0, ncol = ncol.subset) # subset sizes graph
  # whole set data just has cols:
  wholeset.cnames <- c("cx", ugrp, dtfilt, hl.name, "other", "all")
  wholeset.ncol <- length(wholeset.cnames)
  setdf <- matrix(nrow = 0, ncol = wholeset.ncol) # whole set sizes graph
  # append cx data by row to each graph dataset
  for(c in 1:length(cxfn)){
    message("Working on table ",cxfn[c],"...")
    cxtype <- cxfn[c]
    cxdat <- cxsetdat <- c(cxtype)
    if(read.from.disk==T){
      message("Reading table from disk, using cxdl names...")
      fpathc = paste0(readdn, "/", cxfn[c], ".csv")
      mx.all <- fread(fpathc, sep = ',', header = T)
    } else{
      message("Loading table from cxdl...")
      mx.all <- as.matrix(cxdl[[c]])
    }
    nallc <- nrow(mx.all) # all jx count
    mx.all <- as.data.frame(mx.all)
    mx.filt <- mx.all[, colnames(mx.all) %in% dtall]
    mx.filt <- mx.filt[, order(match(colnames(mx.filt), dtall))]
    # Apply prefilter on union groups (order matters!)
    message("Excluding union/filter sets...")
    for(u in 1:length(ugrp)){
      # subset sizes graph data
      cxdat <- cxsetdat <- c(cxdat, nrow(mx.filt[mx.filt[, ugrp[u]]==1,])) # identity check for ugrp[u]
      mx.filt <- mx.filt[!mx.filt[, ugrp[u]]==1,] # condition next groups excluding ugrp[u]
    }
    notherc <- nrow(mx.filt) # count of 'other' jxs
    mx.filt <- mx.filt[, colnames(mx.filt) %in% c(dtfilt)] # mx.filt to only have non-ugrp sets
    # append additional set sizes from dtfilt, after ugrp filt
    for(s in 1:length(dtfilt)){
      cxsetdat <- c(cxsetdat, nrow(mx.filt[mx.filt[, colnames(mx.filt)==dtfilt[s]]==1,])) # whole set sizes data
    }
    
    message("Calculting set overlaps...") # calculate overlaps using length of shared identities
    mx.str <- apply(mx.filt, 1, paste0, collapse=";")
    cnn <- colnames(mx.filt)
    # get subset group sizes
    # note: sizes exclusive, so 'setA' implies 'setA==1' & 'setB==0' etc.
    for(s in 1:length(setlist)){
      # mxff <- mx.filt
      fo <- setlist[[s]]
      if(fo==hl.name){
        cns <- paste0(rep(0, length(cnn)), collapse=";")
        cxsetdat <- c(cxsetdat, length(mx.str[mx.str==cns])) # append hl.name count to whole set data
      }else{
        cns <- paste0(ifelse(cnn %in% fo, 1, 0), collapse=";")
      }
      cxdat <- c(cxdat, length(mx.str[mx.str==cns]))
      #message(s)
    }
    
    # append total "other" and "all" jx counts to both sets
    cxdat <- c(cxdat, notherc, nallc) # subset sizes graph data
    cxsetdat <- c(cxsetdat, notherc, nallc) # sets sizes graph data
    # finally, grow return dfs
    bpdf <- rbind(bpdf, cxdat) # subsets sizes graph data
    setdf <- rbind(setdf, cxsetdat) # sets sizes graph data
    message("Finished cx (index): ", c)
    # message("setdf dim:",dim(setdf))
  }
  message("Making set column names...") # make colnames
  colnames(setdf) <- wholeset.cnames
  colnames(bpdf) <- subset.cnames
  cn <- c(ugrp, unlist(lapply(setlist[c((length(ugrp)+1):length(setlist))], paste0, collapse=";")))
  
  #-------------------
  # rescale jx counts
  #-------------------
  # get % all jx count for both set types
  bpdff <- as.data.frame(bpdf, stringsAsFactors = F) # group sizes
  setdff <- as.data.frame(setdf, stringsAsFactors = F) # subsets
  ndf <- as.data.frame(bpdf[,c(2:(ncol(bpdf)-2))],stringsAsFactors = F)
  gdf <- as.data.frame(setdf[,c(2:(ncol(setdf)-1))], stringsAsFactors = F)
  message("Calculating cross-table jx set size percentages...")
  for(c in 1:ncol(ndf)){ndf[,c] <- as.numeric(ndf[,c])}
  for(c in 1:ncol(gdf)){gdf[,c] <- as.numeric(gdf[,c])}
  if(denom.grp==""){
    message("Error: no denominator group specified! returning...")
  } else{
    if(denom.grp=="all"){
      message("using group 'all' as denominator for fractions...")
      for(i in 1:nrow(ndf)){ndf[i,] <- 100*(ndf[i,]/as.numeric(bpdff$all[i]))}
      for(i in 1:nrow(gdf)){gdf[i,] <- 100*(gdf[i,]/as.numeric(setdff$all[i]))}
    }
    if(denom.grp=="other"){
      message("using group 'other' as denominator for fractions...")
      for(i in 1:nrow(ndf)){ndf[i,] <- 100*(ndf[i,]/as.numeric(bpdff$other[i]))}
      for(i in 1:nrow(gdf)){gdf[i,] <- 100*(gdf[i,]/as.numeric(setdff$other[i]))}
    }
  }
  
  dfplot1 <- ndf # subset graph data
  dfplot2 <- gdf # whole set graph data
  # calculate scale
  message("Rescaling according to 'scale'...")
  if(scale=="log10"){
    dfplot1.scale <- log10(dfplot1+1e-10) # adds offset to avoid Inf values
    dfplot2.scale <- log10(dfplot2+1e-10) # adds offset to avoid Inf values
  }
  
  #-------------
  # return data
  #-------------
  lro <- list(sizes_data_unscaled = list(by_set = setdff, by_subset = bpdff),
              sizes_data_scaled = list(set_sizes = dfplot2.scale,
                                       subset_sizes = dfplot1.scale,
                                       scale_info = paste0("percent;",scale)))
  if(save.lro & !lro.name==""){
    save(lro, file=lro.name)
  }
  return(lro)
}

# 3. make.ggset.final: function to plot jx set data in an upset plot with error bars.
make.ggset.final <- function(setsizes.data, subsetsizes.data, 
                             setsize.grplabels = c("Placenta", "Unexplained", "Other Adult", "Stem Cells", "Developmental"),
                             use.nsubfilt = TRUE, nsub.filt = 20,
                             ugrp = c("paired","gtex"),
                             hl.name = "novel",
                             setsizeplot.dims = unit(c(1, 0, 0, 0), 'cm'),
                             subsetplot.dims = unit(c(0.3, 1.8, -0.6, 3), 'cm'),
                             upsetplot.dims = unit(c(0.5, 1, 0.8, 0), 'cm'),
                             setsizeplot.ybreaks = c(5, 4, 3, 2, 1, 0, -1, -2, -3, -4, -5),
                             setsizeplot.ylab = c("1e5%", "1e4%", "1e3%", "1e2%", "1e1%", "1%", "1e-1%", "1e-2%", "1e-3%", "1e-4%", 
                                                  "1e-5%"),
                             subsetplot.ybreaks = c(2,1,0,-1), 
                             subsetplot.ylab = c("1e2%","1e1%","1%","1e-1%"),
                             subset.barcolor = "limegreen"){
  # make.ggset.final
  # Description: Define and return ggplot2 plot objects for jx sets
  # Arguments:
  #   setsizes (data.frame): Set size dataframe object
  #   subsetsizes (data.frame): Subset sizes dataframe object
  #   setsize.grplabels (vector of strings): Group labels for set sizes barplots (order should match columns in setsizes.data).
  #   use.nsubfilt (boolean): Whether to filter total subsets displayed, ordered from largest to smallest.
  #   nsub.filt (integer): If use.nsubfilt==TRUE, what number of subsets are to be displayed, ordered from largest to smallest.
  #   ugrp (vector of strings): What sets/columns are to be pre-filtered/removed before generating plots.
  #   hl.name (string): Name of set/column data to highlight in plot (currently supports red color highlight).
  #   setsizeplot.dims (unit vector object): Plot margin dimensions of set size barplot (order of values: top, right, bottom, left, unit).
  #   subsetplot.dims (unit vector object): Plot margin dimensions of subset size barplot (order of values: top, right, bottom, left, unit).
  #   upsetplot.dims (unit vector object): Plot margin dimensions of upset plot (order of values: top, right, bottom, left, unit).
  #   setsizeplot.ybreaks (vector of numerics or integers): Units and limits for x-axis on set sizes barplot.
  #   setsizeplot.ylab (vector of strings): Custom labels for x-axis on set sizes barplot.
  #   subsetplot.ybreaks (vector of numerics or integers): Units and limits for x-axis on subset sizes barplot.
  #   subsetplot.ylab (vector of strings): Custom labels for x-axis on subset sizes barplot.
  #   subset.barcolor (string): Specify color for non-highlighted barplots of subset barplot.
  # Returns:
  #   List of ggplot2 graph objects
  
  # dependencies
  require("ggplot2")
  require("gridExtra")
  
  #----------------------
  # Subset Sizes Barplot
  #----------------------
  {
    message("Making subsets barplot...")
    dfi <- subsetsizes.data
    cmeans <- apply(dfi,2,mean) # reorder by colmeans
    dfi <- dfi[,rev(order(cmeans))]
    dfig <- matrix(nrow=0,ncol=3)
    for(i in 1:ncol(dfi)){
      dati <- c()
      dati <- c(colnames(dfi)[i],
                mean(dfi[,i]),
                sd(dfi[,i]))
      dfig <- rbind(dfig, dati)
    }
    dfig <- as.data.frame(dfig, stringsAsFactors = F)
    dfig[,1] <- as.character(dfig[,1])
    dfig[,2] <- as.numeric(dfig[,2])
    dfig[,3] <- as.numeric(dfig[,3])
    colnames(dfig) <- c("group","len","sd")
    
    if(use.nsubfilt){
      message("Applying subset number filter...")
      if(nrow(dfig) > nsub.filt){
        dfig <- dfig[c(1:nsub.filt),] # apply filter for subset barplot
        dfi <- dfi[,colnames(dfi) %in% dfig$group] # apply filter for upset plot
      }
    } else{
      message("Continuing without subset number filter...")
    }
    
    # barplot color scheme
    bpcol <- rep(subset.barcolor, nrow(dfig)) # bar colors
    bpcol[which(dfig$group==hl.name)] <- "red" # 'novel' set bar color
    
    message("Forming ggplot2 plot object...")
    bpgg2.setol <- ggplot(dfig, aes(x = reorder(group, -len), y = len)) + 
      geom_bar(stat="identity", fill=c(bpcol), 
               position=position_dodge(),colour="black") +
      geom_errorbar(aes(ymin = len-sd, ymax = len+sd), width=.2,
                    position=position_dodge(.9)) +
      theme(strip.text.x = element_blank(),
            axis.title.x = element_blank(),
            axis.text.x = element_blank(),
            plot.margin = subsetplot.dims,
            panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      scale_y_continuous(breaks = subsetplot.ybreaks,
                         labels = subsetplot.ylab,
                         limits = c(min(subsetplot.ybreaks), max(subsetplot.ybreaks)))
    
  }
  
  #-------------------
  # Set sizes barplot
  #-------------------
  {
    message("Making set sizes barplot...")
    dfi2 <- setsizes.data
    dfi2 <- dfi2[,!grepl("other", colnames(dfi2))]
    dfig2 <- matrix(nrow=0,ncol=3)
    for(i in 1:ncol(dfi2)){
      dati2<- c()
      dati2 <- c(colnames(dfi2)[i],
                 mean(dfi2[,i]),
                 sd(dfi2[,i]))
      dfig2 <- rbind(dfig2, dati2)
    }
    dfig2 <- as.data.frame(dfig2, stringsAsFactors = F)
    dfig2[,1] <- as.character(dfig2[,1])
    dfig2[,2] <- as.numeric(dfig2[,2])
    dfig2[,3] <- as.numeric(dfig2[,3])
    colnames(dfig2) <- c("group","len","sd")
    # reassign union grp colnames
    if(!ugrp==""){
      for(u in 1:length(ugrp)){
        dfig2[grepl(paste0("^", ugrp[u]), 
                    dfig2[,1]),1] <- colnames(dfi)[grepl(paste0("^",ugrp[u]), 
                                                         colnames(dfi))]
      }
    }
    if(length(setsize.grplabels)==nrow(dfig2)){
      dfig2$group.label <- setsize.grplabels # set the group labels
    } else{
      dfig2$group.label <- dfig2$group
    }
    lvl.order <- rev(order(dfig2$len)) # group order for factors
    dfig2 <- dfig2[lvl.order,]
    dfig2$group <- factor(dfig2$group, levels = dfig2$group) # make group a factor var
    dfig2$group.label <- factor(dfig2$group.label, levels = dfig2$group.label) # make group a factor var
    # barplot color scheme
    dfig2$bp.color <- rep("lightblue", nrow(dfig2))
    dfig2$bp.color[which(dfig2$group==hl.name)] <- "red"
    
    message("Forming ggplot2 plot object...")
    bpgg2.grpsize <- ggplot(dfig2, aes(x = group.label, y = len)) + 
      geom_bar(stat = "identity", fill = c(dfig2$bp.color), 
               position = position_dodge(), colour = "black") +
      scale_y_reverse(breaks = setsizeplot.ybreaks,
                      labels = setsizeplot.ylab,
                      limit = rev(c(min(setsizeplot.ybreaks), 
                                    max(setsizeplot.ybreaks)))) +
      coord_flip() +
      geom_errorbar(aes(ymin = len - sd, ymax = len + sd), 
                    width=.2, position = position_dodge(.9)) +
      theme(strip.text.y = element_blank(),
            axis.title.y = element_blank(),
            plot.margin = setsizeplot.dims,
            panel.grid.major = element_blank(), 
            panel.grid.minor = element_blank(),
            panel.background = element_blank(), 
            axis.line = element_line(colour = "black"))
    
  }
  
  #------------
  # Upset Plot
  #------------
  {
    message("Making upset plot...")
    dff <- matrix(nrow=0,ncol=3)
    cfilt <- colnames(dfi) 
    # use optional ugrp filter
    if(!ugrp==""){
      for(u in 1:length(ugrp)){
        usc <- colnames(dfi)[grepl(paste0("^",ugrp[u]), colnames(dfi))]
        dff <- rbind(dff, c(usc, usc, 1))
        cfilt <- cfilt[!grepl(paste0("^",ugrp[u]), cfilt)]
      }
    }
    # form tall df of all subsets
    for(c in 1:length(cfilt)){
      scl <- unlist(strsplit(cfilt[c],";"))
      for(s in 1:length(scl)){
        dff <- rbind(dff, c(scl[s], cfilt[c], 1))
      }
    }
    dff <- as.data.frame(dff,stringsAsFactors=F)
    colnames(dff) <- c("set","group","proportion")
    dff$proportion <- as.numeric(dff$proportion)
    dff$group <- factor(dff$group, levels=colnames(dfi)) # set subset order
    dff$set <- factor(dff$set, levels=rev(levels(dfig2$group))) # rev group levels from set sizes
    
    # make the coloration var
    classes <- rep("other",nrow(dff))
    classes[which(dff$set==hl.name)] <- "unexplained"
    
    p <- ggplot(dff, aes(x="", y=proportion, fill=classes)) +
      geom_bar(width = 1, stat = 'identity') +
      theme(axis.text = element_blank(),
            axis.ticks = element_blank(),
            panel.grid  = element_blank())
    pie <- p + coord_polar('y', start = 0) + 
      scale_fill_manual(values=c("black","red"))
    
    message("Forming upset plot...")
    upset <- pie +
      facet_grid(rows = vars(dff$set), 
                 cols = vars(dff$group), 
                 switch='y') +
      theme(axis.text = element_blank(),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            strip.text.x = element_blank(),
            strip.text.y = element_blank(),
            strip.background = element_rect(colour = "black", fill = "white",
                                            size = 0.4, linetype = "solid"),
            panel.border = element_rect(color = "black", fill = NA, size = 0.5),
            panel.spacing = unit(0, "lines"),
            legend.position = 'none',
            plot.margin = upsetplot.dims,
            plot.background = element_blank())
  }
  
  message("Plot formation complete. Returning...")
  lr <- list(subsets_plot = bpgg2.setol,
             sets_plot = bpgg2.grpsize,
             upset_plot = upset)
  return(lr)
}

#----------------
# Save functions
#----------------
save(get.datlist, file="get_datlist.rda")
save(jx.calcsets, file="jx_calcsets_function.rda")
save(make.ggset.final, file="make_ggsetplot_function.rda")
