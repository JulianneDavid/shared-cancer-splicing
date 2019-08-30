# shared-cancer-splicing
Code for reproducing analyses and figures for shared alternative cancer splicing paper

# Requirements

software:
Anaconda3
Mmh3
sqlite3
intervaltree
seaborn

Clone snaptron-experiments in directory UTILITIES_DIRECTORY from ____github.


data files:
From https://jhubiostatistics.shinyapps.io/recount/
Under the TCGA tab:
Click “jx_bed” in the junction column to download the TCGA_JUNCTION_BED file, TCGA.junction_id_with_transcripts.bed
Click “jx_cov” in the junction column to download the TCGA_JUNCTION_COVERAGE file, TCGA.junction_coverage.tsv
Click the “link” in the phenotype column to download the TCGA_PHEN tab-separated phenotype file, TCGA.tsv

GTEx tab:
Click “jx_bed” in the junction column to download the GTEx_JUNCTION_BED file, SRP012682.junction_id_with_transcripts.bed
Click “jx_cov” in the junction column to download the GTEx_JUNCTION_COVERAGE file, SRP012682.junction_coverage.tsv
Click “link” in the phenotype column to download the GTEx_PHEN tab-separated phenotype file, SRP012682.tsv

Documentation tab:
Click sample_ids.tsv under the “junction raw coverage file” entry for the RECOUNT_SAMPLE_ID_FILE

GENCODE_ANNOTATION_GTF: download  gencode.v28.annotation.gtf from https://www.gencodegenes.org/human/release_28.html

<<fasta files>>??

CANCER_GENE_CENSUS COSMIC cancer_gene_census.csv from https://cancer.sanger.ac.uk/census

ONCOKB_GENES cancerGeneList.tsv from https://oncokb.org/cancerGenes

MetaSRA run files, processed to remove potential “tumor” samples -- upload to directory SRA_junction_download/MetaSRA_queries

Patient somatic mutation call files:
“Mutation_Packager_Oncotated_Calls” tar.gz for each cancer type from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/


UniProt splicing-associated gene file (tab-separated): keyword:"mRNA splicing [KW-0508]" & Reviewed:yes & organism:"Homo sapiens (Human) [9606]" from https://www.uniprot.org/uniprot/

——
Collect TCGA junctions

1. Specify a directory DB_DIRECTORY to hold the junction database

2. create junction index database by running jx_indexer in “index” mode

jx_indexer.py -d DB_DIRECTORY index -c GTEx_JUNCTION_COVERAGE -C TCGA_JUNCTION_COVERAGE -b GTEx_JUNCTION_BED -B TCGA_JUNCTION_BED -p GTEx_PHEN -P TCGA_PHEN -s RECOUNT_SAMPLE_ID_FILE -g GENCODE_ANNOTATION_GTF

3. run “experiment” mode on junction index database to collect junction files for analysis.
jx_indexer.py -d DB_DIRECTORY experiment -o OUTPUT_DIRECTORY

Within OUTPUT_DIRECTORY, generates:
ALL_JX_DIRECTORY for use in set membership annotation.

NON_CORE_NORMAL_DIRECTORY, and NON_TISSUE_MATCHED_NORMAL_DIRECTORY for use in set membership annotation and Figures S3A, S3B, S4A, and S4B.

COUNTS_PER_SAMPLE_DIRECTORY (junctions not found in core normals) to generate Figure 1B

PREVALENCE_FILE_DIRECTORY (tcga cancer junctions not found in core normals with cancer-type prevalence, per cancer type) to generate Figures 2A, 2B, 2C, and S5

FILTERED_NTM_JX_PER_SAMPLE_DIRECTORY (coverage- and annotation-filtered junctions not found in GTEx or TCGA tissue-matched normal samples) and FILTERED_NCN_JX_PER_SAMPLE_DIRECTORY (coverage- and annotation-filtered junctions not found in core normals) to generate Figure S1A

ALL_JXS_PER_SAMPLE_DIRECTORY containing all junctions for TCGA SKCM normal samples, TCGA SKCM tumor samples, and GTEx bulk skin normal samples to generate Figure S7.

——
Collect SRA junctions

All MetaSRA query files are in this repo under SRA_junction_download/MetaSRA_queries. These files should be downloaded to directory METASRA_QUERY_FILES

python create_snaptron_queries.py -s METASRA_QUERY_FILES -i RECOUNT_SAMPLE_ID_FILE -o LOGGING_OUTPUT_DIRECTORY -p SRA_JUNCTION_OUTPUT_DIRECTORY -u UTILITIES_DIRECTORY

Where UTILITIES_DIRECTORY is where snaptron-experiments was previously cloned; SRA_JUNCTION_OUTPUT_DIRECTORY is where the results will be stored.

Run each SAMPLE_TYPE_indiv_chrom_batch_query_script.sh to collect SRA junctions.

Output will be as follows:
Non-cancer junction results will be saved in 
SRA_JUNCTION_OUTPUT_DIRECTORY/SRA_noncancer_rawresults, hereafter referred to as SNAPTRON_NONCANCER_DIRECTORY
Non-cancer experiment lists will be saved in
SRA_JUNCTION_OUTPUT_DIRECTORY/SRA_noncancer_exptlists, hereafter referred to as SNAPTRON_NONCANCER_EXPTLIST_DIRECTORY 

Cancer junction results will be saved in 
SRA_JUNCTION_OUTPUT_DIRECTORY/SRA_cancer_rawresults, hereafter referred to as SNAPTRON_CANCER_DIRECTORY
Cancer experiment lists will be saved in
SRA_JUNCTION_OUTPUT_DIRECTORY/SRA_cancer_exptlists, hereafter referred to as SNAPTRON_CANCER_EXPTLIST_DIRECTORY 


——
Data preparation
Run 1-read TCGA junction collection 
python collect_1-read_TCGA_jxs.py -C TCGA_JUNCTION_COVERAGE -B TCGA_JUNCTION_BED -o JSON_DIRECTORY

Generates file single-read-tcga-jxs_json.txt in JSON_DIRECTORY for use in set membership annotation runs below.

Add set membership/piechart annotation to full junction files:
For primary data and figures:
python set_membership_annotation.py --db-path DB_DIRECTORY --snaptron-results SNAPTRON_NONCANCER_DIRECTORY -d ALL_JX_DIRECTORY -g NON_CORE_NORMAL_DIRECTORY -p NON_TISSUE_MATCHED_NORMAL_DIRECTORY --gtf-file GENCODE_ANNOTATION_GTF --single-read-jx-json JSON_DIRECTORY/single-read-tcga-jxs_json.txt --cancer-sra-directory SNAPTRON_CANCER_DIRECTORY --cancer-gene-census CANCER_GENE_CENSUS --oncokb-cancer-genes ONCOKB_GENES --min-overall-set-count 1 -o FULL_PIECHART_DIRECTORY

Note: this also creates subdirectories FULL_PIECHART_DIRECTORY/unexplained and FULL_PIECHART_DIRECTORY/developmental containing the relevant subsets of the full files.

Run 2-sample minimum requirement for category membership for Figures S8B and S8C:
python set_membership_annotation.py --db-path DB_DIRECTORY --snaptron-results SNAPTRON_NONCANCER_DIRECTORY -d ALL_JX_DIRECTORY -g NON_CORE_NORMAL_DIRECTORY -p NON_TISSUE_MATCHED_NORMAL_DIRECTORY --gtf-file GENCODE_ANNOTATION_GTF --single-read-jx-json 1_READ_TCGA_JX_JSON --cancer-sra-directory SNAPTRON_CANCER_DIRECTORY --cancer-gene-census CANCER_GENE_CENSUS --oncokb-cancer-genes ONCOKB_GENES --min-overall-set-count 2 -o 2-SAMPLE_PIECHART_DIRECTORY

Note: this also creates subdirectories 2-SAMPLE_PIECHART_DIRECTORY/unexplained and 2-SAMPLE_PIECHART_DIRECTORY/developmental containing the relevant subsets of the full files.

——
Data analyses:

Analyze set memberships
python set_membership_analysis.py -o OUTPUT_DIRECTORY -s FULL_PIECHART_DIRECTORY

Note: creates directory FULL_PIECHART_DIRECTORY/true_TCGA_prevalence_files populated with non-core-normal junctions for use in Figures S1B and S1C

count SRA experiments:
python count_unique_SRA_expts.py -c SNAPTRON_CANCER_EXPTLIST_DIRECTORY -n SNAPTRON_NONCANCER_EXPTLIST_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

——
Figure generation:

Fig 1A:
python fig1A_overall_set_barplots.py -s FULL_PIECHART_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 1B:
python fig1B_ncn_jx_counts_per_sample.py -j COUNTS_PER_SAMPLE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 1C:
python fig1C_overall_set_prevalence_boxplot.py --set-memberships FULL_PIECHART_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 2A:
python fig2A_TCGA_heatmap.py -d PREVALENCE_FILE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 2B:
python fig2B_TCGA_subtypes_heatmap.py -d PREVALENCE_FILE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 2C:
python fig2C_cell_of_origin_heatmap.py  -d PREVALENCE_FILE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY --snaptron-results SNAPTRON_NONCANCER_DIRECTORY -e SNAPTRON_NONCANCER_EXPTLIST_DIRECTORY

Fig 3A:--already added by sean on github

Fig 3B:
python fig3B_antisense_boxplot.py -o FIGURE_OUTPUT_DIRECTORY -s FULL_PIECHART_DIRECTORY

Fig S1A:
python fig1B_ncn_jx_counts_per_sample.py -j FILTERED_NCN_JX_PER_SAMPLE_DIRECTORY -p FILTERED_NTM_JX_PER_SAMPLE_DIRECTORY -d -g THYM CESC UVM DLBC --prepared-sort-order -o FIGURE_OUTPUT_DIRECTORY

Fig S1B and S1C:
python figS1BC_S9BC_junction_sharedness.py -d FULL_PIECHART_DIRECTORY/true_TCGA_prevalence_files populated -o FIGURE_OUTPUT_DIRECTORY

Note: requires contents of “true_TCGA_prevalence_files” directory created by set_membership_analysis.py

Fig S2:
python figS2_set_prevalences_per_cancer.py -s FULL_PIECHART_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig S3A & S3B:
R SF_mutation_script.R

Fig S4A & S4B:
R SF_sharedness_plots.R
Ben’s instructions:
To generate splicing factor plots, ensure you have the following dependencies installed for R:
Data.table
ggplot2
gridExtra
RColorBrewer
tidyverse
dplyr

1. download the TCGA mutation files from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ for each cancer type and move into the sub-directory ./TCGA_mut_files/, where ./ is the location of "SF_mutation_script.R" and "SF_sharedness_plots.R"

2. Copy files from NON_CORE_NORMAL_DIRECTORY to a new directory within ./ called "non_core_jxns/". Rename files from full cancer names to match TCGA cancer codes.

3. Copy files from [PIECHART DIRECTORY] to a new directory within ./ called "non_core_jxn_annotations/". Rename files from full cancer names to match TCGA cancer codes.

4. run SF_mutation_script.R

Generates "fig_s3a.jpg", "fig_s3b.jpg" and shared_jxn_df.txt. shared_jxn_df.txt is used downstream for SF_sharedness_plots.R

If run interactively, also presents statistics for jxn comparisons between patients with and without splicing factor mutations across cancers.

5. run SF_sharedness_plots.R

Generates "fig_s4a.jpg" and "fig_s4b.jpg"

If run interactively, also presents statistics for jxn comparisons between patients with and without splicing factor mutations across cancers.
---

Fig S5:
python figS5_S9A_SRA_cancer_shared_prevalence.py --snaptron-results SNAPTRON_CANCER_DIRECTORY -e SNAPTRON_CANCER_EXPTLIST_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY -d PREVALENCE_FILE_DIRECTORY

Fig S6:
python figS6_full_TCGA_SRA_heatmap.py  -d PREVALENCE_FILE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY --snaptron-results SNAPTRON_NONCANCER_DIRECTORY -e SNAPTRON_NONCANCER_EXPTLIST_DIRECTORY

Fig S7:
python figS7_SKCM_jx_similarity.py --snaptron-results SNAPTRON_NONCANCER_DIRECTORY -d ALL_JXS_PER_SAMPLE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig S8A:
Fig S8B:
Fig S8C: -- already added by sean on github

Fig S9A:
python figS5_S9A_SRA_cancer_shared_prevalence.py --snaptron-results SNAPTRON_CANCER_DIRECTORY -e SNAPTRON_CANCER_EXPTLIST_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY -d UNEXPLAINED_PIECHART_DIRECTORY  --unexplained-junctions

Fig S9B and S9C:
python figS1BC_S9BC_junction_sharedness.py -d UNEXPLAINED_PIECHART_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY


# Generate Junction Upset Plot Figures
To generate junction upset plots, ensure you have the following dependencies installed for R:

1. [ggplot2](https://cran.r-project.org/web/packages/ggplot2/index.html)
2. [data.table](https://cran.r-project.org/web/packages/data.table/index.html)
3. [gridExtra](https://cran.r-project.org/web/packages/gridExtra/index.html)

Next, download the cancer-specific junction tables from the data resource to a subdirectory in your working directory.

Instantiate the functions in `upset_functions.R`, and make a new subdirectory (default name: `cxdat`), and run `get.datlist()`, specifying the name of the subdirectory containing the junction data tables. Once this has run, the working directory should contain a large list object containing the junction data tables, as well as a newly populated subdirectory (default: `cxdat`). Depending on your session memory limit, you may need to remove the object `cxdl` from your environment before proceeding.

Follow the remaining steps in `upset_data.R` to: Define the new `developmental` junction category; Run `jx.calcsets()` for the supplemental data, then for the main figure data; Replace the `unexplained` category in the main figure data with the supplemental figure data; Generate the supplemental data tables. To save memory, tables are read processively and processed.

Next, follow the steps in `upset_plots_original.R` to: Specify the data table containing the whole/total group sets and subset group sets, respectively; Exclude the pre-filtered set columns from visualization; Properly order whole set columns by mean value; Automatically generate the axis limits and labels; Specify the group labels (should match whole set column labels after ordering on means); Generate the plot data objects; Generate the final plot images.

Following these steps should generate both a main and a supplemental upset plot with annotation set and subset abundances. It may be necessary to experiment with different plot and image dimensions in `upset_plots_original.R` and `make.ggset.final()` to refine the figure properties. For additional information, refer to the function docstrings and script comments.




