#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(ggplot2)
library(survival)
library(survminer)
library(gridExtra)
library(grid)

setwd(args[1])
dat <- read.csv("OV_mslnjx_survival_data.csv", header = T)
surv_model <- survfit(Surv(time, death)~as.factor(with_MSLN_jx), data=dat)
surv_plot <- ggsurvplot(surv_model,
                        conf.int = F,
                        risk.table = T,
                        legend.title = "TCGA ovarian\ncancer samples",
                        legend.labs = c("without jx","with target MSLN jx"),
                        palette = c("#0072B2","#D55E00"),
                        pval = F,
                        break.time.by = 500,
                        xlab="Time (days from diagnosis)")

# adjust size as needed, also accepts PDF
ggsave("OV_msln_survival_plot.jpg",
       plot = print(surv_plot),
       device = "jpeg",
       units = "in",
       height = 6,
       width = 10)
