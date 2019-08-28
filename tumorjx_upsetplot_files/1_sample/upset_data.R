#!/usr/bin/env/R

# Calculate upset data from jx db
# Code author: SKM
# Contact: maden@ohsu.edu

#-----------------------
# Load Utility Functions
#-----------------------
load("get_datlist.rda")
load("jx_calcsets_function.rda")
load("make_ggsetplot_function.rda")

#---------------------
# Read in Cancer Data
#---------------------

dn = "1-sample"
fl = list.files(dn)
fl.filt = fl[grepl(".*\\.csv", fl)]
cxdl <- get.datlist(dn, fl = fl.filt)

# define developmental category
for(c in 1:length(cxdl)){
  cxx <- cxdl[[c]]
  cxx$developmental <- ifelse(cxx$sra_developmental==1 & cxx$sra_placenta==0, 1, 0)
  cxdl[[c]] <- cxx
  message(c)
}

save(cxdl, file = "listdf_allcx-devadd-augfinal.rda")

# save new tables in dir
library(data.table)
ddn = "cxdat"
for(c in 1:length(cxdl)){
  # write.csv(cxdl[[c]], file = paste0(ddn, "/", names(cxdl)[c],".csv"))
  fwrite(cxdl[[c]], file = paste0(ddn, "/", names(cxdl)[c],".csv"), sep = ",", row.names = F)
  message("finished writing table ",c,", name: ", names(cxdl)[c])
}

#--------------------
# Get Upset Plot Data
#--------------------
# 1. Expanded data (for supplement figure)
cxm.lf = list.files("cxdat")
dtall.sf <- c("paired", "gtex", "sra_stemcells", "sra_adult",
                "sra_embryo_ectoderm", "sra_embryo_embryo", "sra_embryo_lateembryo", 
                "sra_embryo_mesenchyme", "sra_embryo_myoblast", 
                "sra_neonate_fetal", "sra_zygote", "sra_oocyte", "sra_placenta", 
              "developmental")
cxfndat = gsub("\\.csv", "", cxm.lf)

lro.supp <- jx.calcsets(cxfn = cxfndat, read.from.disk = T, readdn = "cxdat", 
                        scale = "log10", dtall = dtall.sf, 
                        hl.name = "unexplained",
                         ugrp = c("paired", "gtex"), 
                        denom.grp = "other",
                        save.lro = TRUE, 
                         lro.name = "jx-datlist-suppfig_final.rda")

# 2. Consolidated set data (for main figure)
# define developmental (non-placenta) column
# generate data
# cxfn <- list.files("cxdat")
# cxfn <- gsub("\\.csv", "", cxfn)
cxfn = cxfndat
dtall.figfinal <- c("paired", "gtex", "sra_stemcells", "sra_adult", "sra_placenta", "developmental")
lro.main <- jx.calcsets(cxfn = cxfn, readdn = "cxdat", read.from.disk = T, 
                        dtall = dtall.figfinal, ugrp = c("paired", "gtex"),
                         hl.name = "unexplained", denom.grp = "other", save.lro = TRUE, 
                         lro.name = "jx-datlist-mainfig_final.rda")
# correct unexplained data from expanded data
lro.main$sizes_data_unscaled$by_set$unexplained <- lro.supp$sizes_data_unscaled$by_set$unexplained
lro.main$sizes_data_unscaled$by_subset$unexplained <- lro.supp$sizes_data_unscaled$by_subset$unexplained
lro.main$sizes_data_scaled$set_sizes$unexplained <- lro.supp$sizes_data_scaled$set_sizes$unexplained
lro.main$sizes_data_scaled$subset_sizes$unexplained <- lro.supp$sizes_data_scaled$subset_sizes$unexplained
save(lro.supp, file="jx-datlist-suppfig_final.rda")

#------------------
# Write data tables
#------------------
write.csv(lro.supp$sizes_data_unscaled$by_set, "supptable-supplfig-dat_wholesetsizes.csv",  row.names=F)
write.csv(lro.supp$sizes_data_unscaled$by_subset, "supptable-supplfig-dat_subsetsizes.csv",  row.names=F)
write.csv(lro.main$sizes_data_unscaled$by_set, "supptable-mainfig-dat_wholesetsizes.csv",  row.names=F)
write.csv(lro.main$sizes_data_unscaled$by_subset, "supptable-mainfig-dat_subsetsizes.csv",  row.names=F)


