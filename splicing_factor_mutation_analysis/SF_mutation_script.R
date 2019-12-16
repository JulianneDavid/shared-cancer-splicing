#!/usr/bin/env/R

# SF_mutation_script.R
# Identifies patients with splicing factor mutations to compare jxn count and sharedness
# uses TCGA mutation files from GDAC firehose, stored in dir ./TCGA_mut_files/
# uses non core junction files produced by running jx_indexer in “experiment” mode, stored in dir ./non_core_jxns/
# uses non core junction annotation files produced by running set_membership_annotation.py on the output of jx_indexer in “experiment” mode, stored in dir ./non_core_jxn_annotations/
# Code author: Ben Weeder
# Contact: weeder@ohsu.edu

# load required libraries
require(data.table)
require(ggplot2)
require(gridExtra)
require(RColorBrewer)
require(tidyverse)
require(dplyr)

# pull in splicing factor genes from uniprot
SF_genes <- fread("uniprot-SF-genes.tsv", select = c("Entry name", "Organism", "Gene names"))
SF_genes <- subset(SF_genes, SF_genes$Organism == "Homo sapiens (Human)") # only needed if human filter not used on dl
SF_genes$`Entry name` <- gsub("_HUMAN", "", SF_genes$`Entry name`) # remove _human so that gene names match with other files

# pull in TCGA mutation calls
TCGA_mut_dirs <- list.files("./TCGA_mut_files")
TCGA_mut_dir_paths <- paste("./TCGA_mut_files/", TCGA_mut_dirs, sep = "")
TCGA_labels <- unlist(lapply(TCGA_mut_dirs, function(x){strsplit(x, split = "_")[[1]][2]}))
TCGA_labels <- unlist(lapply(TCGA_labels, function(x){strsplit(x, split = "[.]")[[1]][1]}))
TCGA_labels

# pull in junction calls
jxn_files <- list.files("./non_core_jxns")
jxn_file_paths <- paste("./non_core_jxns/", jxn_files, sep = "")
jxn_labels <- unlist(lapply(jxn_files, function(x){strsplit(x, split = "_")[[1]][2]}))
jxn_labels <- jxn_labels[!(jxn_labels %in% "MESO")] # no mutation calls here?

# pull in gene maps for junctions
# note that jx annotation file names were changed to match TCGA cancer codes
jxn_map_files <- list.files("./non_core_jxn_annotations")
jxn_map_file_paths <- paste("./non_core_jxn_annotations/", jxn_map_files, sep = "")
jxn_map_labels <- unlist(lapply(jxn_map_files, function(x){strsplit(x, split = "_")[[1]][2]}))
jxn_map_labels <- jxn_map_labels[!(jxn_labels %in% "MESO")] # ignore as above

# compile list of mutations for all patients with annotations in TCGA
Compiled_mut_list <- list()

# loop through mut_file_paths for each disease
for (d in 1:length(TCGA_mut_dirs)){
  pt_files <- list.files(TCGA_mut_dir_paths[d])[which(!(list.files(TCGA_mut_dir_paths[d]) %in% c("MANIFEST.txt")))]
  pt_ids <- unlist(lapply(pt_files, function(x){strsplit(x, split = "[.]")[[1]][1]}))
  pt_file_paths <- paste(TCGA_mut_dir_paths[d],"/", pt_files, sep = "")
  
  pt_mut_list <- list()
  for(f in 1:length(pt_files)){
    pt_mut_list[[f]] <- fread(file = pt_file_paths[f], sep = "\t", select = c("Hugo_Symbol","Variant_Classification"))
    npts <- dim(pt_mut_list[[f]])[1]
    pt_mut_list[[f]]$Disease <- TCGA_labels[d]
    pt_mut_list[[f]]$pt_id <- rep(pt_ids[f], npts)
    merged_pt_muts <- rbindlist(pt_mut_list, fill = T)
    Compiled_mut_list[[d]] <- merged_pt_muts
    names(Compiled_mut_list)[d] <- TCGA_labels[d]
  }
}

# Compile across all cancers
full_dat <- rbindlist(Compiled_mut_list, fill = T)

# trim to patient specific, not sample specific id's
full_dat$tcga_id <- strtrim(full_dat$pt_id, 12)

# make list for number of occurences per gene mutation
SF_mut_sums <- list()
SF_gene_list <- c(SF_genes$`Entry name`, "U2AF1", "SF3B1", "TADA1", "PPP2R1A", "IDH1") # all in uniprot + "new" identified by Kahles et al.
SF_gene_list <- unique(SF_gene_list)
# SF_gene_list <- SF_genes$`Entry name`
for(sf in 1:length(SF_gene_list)){
  SF_mut_sums[[sf]] <- data.table(SF = SF_gene_list[sf],
                                  nEvent = uniqueN(full_dat$pt_id[full_dat$Hugo_Symbol == SF_gene_list[sf] & full_dat$Variant_Classification != "Silent"]),
                                  nDisease = uniqueN(full_dat$Disease[full_dat$Hugo_Symbol == SF_gene_list[sf] & full_dat$Variant_Classification != "Silent"]))
}

mut_sums <- rbindlist(SF_mut_sums)
# filter out those with no ocurrences
mut_sums <- subset(mut_sums, mut_sums$nEvent > 0)
# order by most common
mut_sums <- mut_sums[order(mut_sums$nEvent, decreasing = T),]

# flag those specifically from Kahles et al.
prev_list <-  c("U2AF1", "SF3B1", "TADA1", "PPP2R1A", "IDH1")
mut_sums$prev_referenced <- rep("no", dim(mut_sums)[1])
mut_sums$prev_referenced[mut_sums$SF %in% prev_list] <- "yes"

# compile mutation prevalence across cancer types
SF_percentages <- rep(NA, length(TCGA_labels))
prev_SF_percentages <- rep(NA, length(TCGA_labels))

for (dis in 1:length(TCGA_labels)){
  tmp_dat <- subset(full_dat, full_dat$Disease == TCGA_labels[dis])
  n_patients <- uniqueN(tmp_dat$pt_id)
  n_SFs <- uniqueN(tmp_dat$pt_id[tmp_dat$Variant_Classification != "silent" & !duplicated(tmp_dat[,c("Hugo_Symbol", "pt_id")]) & (tmp_dat$Hugo_Symbol %in% SF_gene_list)])
  n_prevSFs <- uniqueN(tmp_dat$pt_id[tmp_dat$Variant_Classification != "silent" & !duplicated(tmp_dat[,c("Hugo_Symbol", "pt_id")]) & (tmp_dat$Hugo_Symbol %in% c("U2AF1", "SF3B1", "TADA1", "PPP2R1A", "SRSF2"))])
  SF_percentages[dis] <- n_SFs/n_patients
  prev_SF_percentages[dis] <- n_prevSFs/n_patients
}

# generate table mapping cancer type and SF prevalence
SF_by_cancer <- data.frame(cancer = TCGA_labels, SF_percentages = SF_percentages, established_SF_percentages = prev_SF_percentages)

# load in junction data 
compiled_jxn_list <- list()
for (f in 1:length(jxn_file_paths)){
  tmp_dat <- fread(jxn_file_paths[f])
  tmp_label <- jxn_labels[f]
  compiled_jxn_list[[f]] <- data.frame(tmp_dat, Disease = tmp_label)
}

# combine across cancers and patients
full_jx_by_pt <- rbindlist(compiled_jxn_list)

# compile map of junctions to genes affected
compiled_jxn_map_list <- list()
for(f in 1:length(jxn_map_file_paths)){
  compiled_jxn_map_list[[f]] <- fread(jxn_map_file_paths[f], select = c("jx", "coding_regions"))
}

# merge
compiled_jxn_map <- rbindlist(compiled_jxn_map_list)
# remove duplicates
compiled_jxn_map <- unique(compiled_jxn_map)

# merge jx by pt and jx to gene map... for only coding jx's code is commented out, but can be swapped
compiled_jxn_df <- merge(full_jx_by_pt, compiled_jxn_map, all.x = T)

jxn_aggregates_by_pt <- compiled_jxn_df %>% group_by(tcga_id, Disease) %>% summarise(total_jxns = length(jx))
jxn_aggregates_by_pt <- jxn_aggregates_by_pt[order(jxn_aggregates_by_pt$Disease, jxn_aggregates_by_pt$total_jxns),]
jxn_aggregates_by_pt$pt_order <- 1:dim(jxn_aggregates_by_pt)[1]

# filter out patients with no RNAseq (no jxn calls)
full_dat <- subset(full_dat, full_dat$tcga_id %in% unique(jxn_aggregates_by_pt$tcga_id))
# filter out patients with no mutation calls
jxn_aggregates_by_pt <- subset(jxn_aggregates_by_pt, jxn_aggregates_by_pt$tcga_id %in% unique(full_dat$tcga_id))

# note that some TCGA here have "no" for cancer type... some issue there... looks like all are UVM in the full_dat
jxn_aggregates_by_pt$Disease <- factor(jxn_aggregates_by_pt$Disease, levels = levels(SF_by_cancer$cancer))
jxn_aggregates_by_pt$tcga_id <- factor(jxn_aggregates_by_pt$tcga_id, levels = jxn_aggregates_by_pt$tcga_id[1:dim(jxn_aggregates_by_pt)[1]])

# add SF status for each patient
jxn_aggregates_by_pt$SF_status <- rep(NA, dim(jxn_aggregates_by_pt)[1])
for(pt in 1:dim(jxn_aggregates_by_pt)[1]){
  if(pt%%1000 == 0){print(pt)}
  if (TRUE %in% (full_dat$tcga_id == jxn_aggregates_by_pt$tcga_id[pt] & 
                 full_dat$Variant_Classification != "Silent" & 
                 full_dat$Hugo_Symbol %in% SF_gene_list)){
    jxn_aggregates_by_pt$SF_status[pt] <- "yes"
  }
  else{jxn_aggregates_by_pt$SF_status[pt] <- "no"}
}

# add est_SF_status for each patient (based on genes from Kahles et al.)
jxn_aggregates_by_pt$est_SF_status <- rep(NA, dim(jxn_aggregates_by_pt)[1])
for(pt in 1:dim(jxn_aggregates_by_pt)[1]){
  if(pt%%1000 == 0){print(pt)}
  if (TRUE %in% (full_dat$tcga_id == jxn_aggregates_by_pt$tcga_id[pt] & 
                 full_dat$Variant_Classification != "Silent" & 
                 full_dat$Hugo_Symbol %in% prev_list)){
    jxn_aggregates_by_pt$est_SF_status[pt] <- "yes"
  }
  
  else{jxn_aggregates_by_pt$est_SF_status[pt] <- "no"}
}

# change SF status to yes/no, make graph grouping
jxn_aggregates_by_pt$SF_status <- factor(jxn_aggregates_by_pt$SF_status, levels = c("yes", "no"))
jxn_aggregates_by_pt$graph_grp <- paste(jxn_aggregates_by_pt$Disease, jxn_aggregates_by_pt$SF_status, sep = "_")

# change est_SF status to yes/no, make graph grouping
jxn_aggregates_by_pt$est_SF_status <- factor(jxn_aggregates_by_pt$est_SF_status, levels = c("yes", "no"))
jxn_aggregates_by_pt$est_graph_grp <- paste(jxn_aggregates_by_pt$Disease, jxn_aggregates_by_pt$est_SF_status, sep = "_")

# order diseases by median junctions
disease_order <- jxn_aggregates_by_pt %>% group_by(Disease) %>% summarise(median_jxns = median(total_jxns))
disease_order <- disease_order[order(disease_order$median_jxns, decreasing = T),]

# reorder data for proper graphing
jxn_aggregates_by_pt$Disease <- factor(jxn_aggregates_by_pt$Disease, levels = disease_order$Disease)
jxn_aggregates_by_pt <- jxn_aggregates_by_pt[order(jxn_aggregates_by_pt$Disease, jxn_aggregates_by_pt$SF_status),]
jxn_aggregates_by_pt$graph_grp <- factor(jxn_aggregates_by_pt$graph_grp, levels = unique(jxn_aggregates_by_pt$graph_grp[1:length(jxn_aggregates_by_pt$graph_grp)]))
jxn_aggregates_by_pt <- jxn_aggregates_by_pt[order(jxn_aggregates_by_pt$Disease, jxn_aggregates_by_pt$est_SF_status),]
jxn_aggregates_by_pt$est_graph_grp <- factor(jxn_aggregates_by_pt$est_graph_grp, levels = unique(jxn_aggregates_by_pt$est_graph_grp[1:length(jxn_aggregates_by_pt$est_graph_grp)]))

# create color pallete for plotting TCGA data
color_map <- data.frame(Disease = TCGA_labels, color_mapping=rep(NA, length(TCGA_labels)))
color_map$color_mapping <- c("#b79400",
                             "#fff4f2",
                             "#cb0162",
                             "#ffa756",
                             "#01386a",
                             "#95d0fc",
                             "#1e488f",
                             "#107ab0",
                             "#9e43a2",
                             "#90e4c1",
                             "#ff474c",
                             "#ffc5cb",
                             "#f7879a",
                             "#ff7855",
                             "#c292a1",
                             "#b7c9e2",
                             "#e4cbff",
                             "#966ebd",
                             "#d8863b",
                             "#607c8e",
                             "#dbb40c",
                             "#840000",
                             "#d5ffff",
                             "#029386",
                             "#b0dd16",
                             "#7bc8f6",
                             "#be0119",
                             "#fdee73",
                             "#b9a281",
                             "#ffd8b1",
                             "#ff9408",
                             "#048243")

# order colors to match cancer types
color_map$Disease <- factor(color_map$Disease, levels = disease_order$Disease)
color_map <- color_map[order(color_map$Disease),]

# order SF by cancer to match jxn_aggregates_by_pt
SF_by_cancer$cancer <- factor(SF_by_cancer$cancer, levels = disease_order$Disease)

# plot junctions per patient by cancer
p_jxns_by_cancer_sep <- ggplot(data = jxn_aggregates_by_pt, aes(x=graph_grp, y=log(total_jxns)))+
  stat_boxplot(geom="errorbar", width = .35)+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width=.07), aes(fill=Disease, color=Disease), pch=21)+
  scale_x_discrete(labels = NULL)+
  scale_fill_manual(values =  c(color_map$color_mapping))+
  scale_color_manual(values = c(color_map$color_mapping))+
  guides(fill=F, color=F)+
  theme_classic()+
  xlab("")+
  ylab("Novel Junctions per Patient (log scale) For Those \n With v. Without Splicing Associated Mutation")
# p_jxns_by_cancer_sep

# plot SF mutation prevalence
p_muts_by_cancer <- ggplot(SF_by_cancer, aes(x = cancer, y=SF_percentages*100))+
  geom_bar(stat = "identity", color = "black", fill = color_map$color_mapping)+
  theme_classic()+
  scale_x_discrete()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size = 16))+
  xlab("TCGA Cancer Type")+
  ylab("Percentage of Patients With at Least \n One Splicing Associated Mutation")
# p_muts_by_cancer

# combine the two graphs
fig_s1f <- grid.arrange(p_jxns_by_cancer_sep,p_muts_by_cancer+scale_y_reverse(limits=c(100,0)), ncol = 1)
ggsave(filename = "fig_s1e.jpg", plot = fig_s1e, device = "jpeg", width = 8, height = 4, units = "in")

# look at statistical significant interaction terms
adj_pvals <- p.adjust(summary(lm(log(total_jxns)~Disease*SF_status-1, data = jxn_aggregates_by_pt))$coefficients[,"Pr(>|t|)"], method = "BH")
adj_pvals[adj_pvals<=.05] # only sig. interaction term is BRCA

# plot above based on kahles et al. SF's
p_jxns_by_cancer_est_sep <- ggplot(data = jxn_aggregates_by_pt, aes(x=est_graph_grp, y=log(total_jxns)))+
  stat_boxplot(geom="errorbar", width = .35)+
  geom_boxplot(outlier.alpha = 0)+
  geom_point(position = position_jitter(width=.07), aes(fill=Disease, color=Disease), pch=21)+
  scale_x_discrete(labels = NULL)+
  scale_fill_manual(values =  c(color_map$color_mapping))+
  scale_color_manual(values = c(color_map$color_mapping))+
  guides(fill=F, color=F)+
  theme_classic()+
  xlab("")+
  ylab("Novel Junctions per Patient (log scale) For Those \n With v. Without Pre-Identified Splicing Related Mutation")
# p_jxns_by_cancer_est_sep

# plot est_SF prevalence by cancer
p_est_muts_by_cancer <- ggplot(SF_by_cancer, aes(x = cancer, y=established_SF_percentages*100))+
  geom_bar(stat = "identity", color = "black", fill = color_map$color_mapping)+
  theme_classic()+
  scale_x_discrete()+
  theme(axis.text.x = element_text(angle = 60, hjust = 1),
        axis.title.x = element_text(size = 16))+
  xlab("TCGA Cancer Type")+
  ylab("Percentage of Patients With at Least \n One Pre-Identified Splicing Related Mutation")
# p_est_muts_by_cancer

# combine the two graphs
fig_s1g <- grid.arrange(p_jxns_by_cancer_est_sep,p_est_muts_by_cancer+scale_y_reverse(limits=c(100,0)), ncol = 1)
ggsave(filename = "fig_s1f.jpg", plot = fig_s1f, device = "jpeg", width = 8, height = 4, units = "in")

# look at statistical significance
adj_est_pvals <- p.adjust(summary(lm(log(total_jxns)~Disease*est_SF_status-1, data = jxn_aggregates_by_pt))$coefficients[,"Pr(>|t|)"], method = "BH")
adj_est_pvals[adj_est_pvals<=.05] # no sig. interaction terms

#### reformat data at pt. level without junction aggregation (used for assessing sharedness)
# generate data frame with patient SF status
pt_status <- data.frame(tcga_id = unique(full_dat$tcga_id), Disease = full_dat$Disease[which(!duplicated(full_dat$tcga_id))])
pt_status$SF_mut <- rep(NA, nrow(pt_status))
pt_status$est_SF_mut <- rep(NA, nrow(pt_status))

# for each patient loop through and identify status for each criteria set
for (i in 1:nrow(pt_status)){
  tmp_genes <- full_dat$Hugo_Symbol[which(full_dat$Variant_Classification != "Silent" & full_dat$tcga_id == pt_status$tcga_id[i])]
  if (length(intersect(tmp_genes, SF_gene_list)) > 0){
    pt_status$SF_mut[i] <- "yes"
  }
  else {
    pt_status$SF_mut[i] <- "no"
  }
  
  if (length(intersect(tmp_genes, prev_list)) > 0){
    pt_status$est_SF_mut[i] <- "yes"
  }
  else {
    pt_status$est_SF_mut[i] <- "no"
  }
  if (i %% 1000 == 0){print(i)}
}

# create cancer + status var for plotting
pt_status$p_group <- paste(pt_status$Disease, pt_status$SF_mut, sep = "_")
pt_status$p_est_group <- paste(pt_status$Disease, pt_status$est_SF_mut, sep = "_")

# pull only shared junctions in the set for downstream analysis
shared_jxns <- subset(full_jx_by_pt,full_jx_by_pt$jx %in% full_jx_by_pt$jx[duplicated(full_jx_by_pt$jx)])
shared_jxns$plt_group <- rep(NA, nrow(shared_jxns))
shared_jxns$plt_est_group <- rep(NA, nrow(shared_jxns))
unique_pts <- unique(shared_jxns$tcga_id)

# for each unique patient ID, loop
for (i in 1:length(unique_pts)){
  # if data available for patient...
  if (unique_pts[i] %in% pt_status$tcga_id){
    plt_group <- pt_status$p_group[which(pt_status$tcga_id == unique_pts[i])]
    est_plt_group <- pt_status$p_est_group[which(pt_status$tcga_id == unique_pts[i])]
  }
  else{
    plt_group <- NA
    est_plt_group <- NA
  }
  
  # populate shared jxn table
  shared_jxns$plt_group[which(shared_jxns$tcga_id == unique_pts[i])] <- plt_group
  shared_jxns$plt_est_group[which(shared_jxns$tcga_id == unique_pts[i])] <- est_plt_group
  
  if (i %% 100 == 0){print(i)}
}

#
write.table(shared_jxns, file = "shared_jxn_df.txt", sep = "\t", quote = F, row.names = F)
# shared_jxns <- read.table("shared_jxn_df.txt", sep = "\t", header = T, stringsAsFactors = F)

