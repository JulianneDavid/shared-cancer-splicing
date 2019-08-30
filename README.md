# shared-cancer-splicing
Code here reproduces analyses and figures from "Putatively cancer-specific alternative splicing is shared across patients and present in developmental and other non-cancer cells," by David et al.

## Requirements

### Software
We ran Python scripts with Python 3.6.6, which we installed as part of [Anaconda](https://www.anaconda.com/) 5.2.0. We made use of the third-party modules `intervaltree` 2.1.0 and `seaborn` 0.9.0.

We ran R scripts with R 3.6.1 and made use of `Data.table` 1.12.2, `ggplot2` 3.2.1, `gridExtra` 2.3, `RColorBrewer` 1.1-2, and `dplyr` 0.8.3.

We further used [Snaptron](http://snaptron.cs.jhu.edu/) via its `qs` utility. To obtain it, run
```cd UTILITIES_DIRECTORY
git clone https://github.com/ChristopherWilks/snaptron-experiments
```

### Data
We downloaded exon-exon junction BEDs for GTEx and TCGA and accompanying metadata from [recount2](https://jhubiostatistics.shinyapps.io/recount/):
* [GTEx junction IDs](http://duffel.rail.bio/recount/SRP012682/SRP012682.junction_id_with_transcripts.bed.gz) (`GTEX_JUNCTION_BED`)
* [GTEx samples with junctions](http://duffel.rail.bio/recount/SRP012682/SRP012682.junction_coverage.tsv.gz) (`GTEX_JUNCTION_COVERAGE`)
* [GTEx phenotype file](http://duffel.rail.bio/recount/SRP012682/SRP012682.tsv) (`GTEX_PHEN`)
* [TCGA junction IDs](http://duffel.rail.bio/recount/TCGA/TCGA.junction_id_with_transcripts.bed.gz) (`TCGA_JUNCTION_BED`)
* [TCGA samples with junctions](http://duffel.rail.bio/recount/TCGA/TCGA.junction_coverage.tsv.gz) (`TCGA_JUNCTION_COVERAGE`)
* [TCGA phenotype file](http://duffel.rail.bio/recount/TCGA/TCGA.tsv) (`TCGA_PHEN`)
* [recount sample IDs](https://jhubiostatistics.shinyapps.io/recount/sample_ids.tsv) (`RECOUNT_SAMPLE_IDS`)

We also used:
* [GENCODE v28](ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_28/gencode.v28.annotation.gtf.gz) (`GENCODE_GTF`)
* [COSMIC's cancer gene census](https://cancer.sanger.ac.uk/cosmic/login?d=1&r_url=%2Fcosmic%2Fcensus%2Fall%3Fhome%3Dy%26name%3Dall%26tier%3D%26sEcho%3D1%26iColumns%3D20%26sColumns%3D%26iDisplayStart%3D0%26iDisplayLength%3D25%26mDataProp_0%3D0%26sSearch_0%3D%26bRegex_0%3Dfalse%26bSearchable_0%3Dtrue%26bSortable_0%3Dtrue%26mDataProp_1%3D1%26sSearch_1%3D%26bRegex_1%3Dfalse%26bSearchable_1%3Dtrue%26bSortable_1%3Dtrue%26mDataProp_2%3D2%26sSearch_2%3D%26bRegex_2%3Dfalse%26bSearchable_2%3Dtrue%26bSortable_2%3Dtrue%26mDataProp_3%3D3%26sSearch_3%3D%26bRegex_3%3Dfalse%26bSearchable_3%3Dtrue%26bSortable_3%3Dtrue%26mDataProp_4%3D4%26sSearch_4%3D%26bRegex_4%3Dfalse%26bSearchable_4%3Dtrue%26bSortable_4%3Dtrue%26mDataProp_5%3D5%26sSearch_5%3D%26bRegex_5%3Dfalse%26bSearchable_5%3Dtrue%26bSortable_5%3Dtrue%26mDataProp_6%3D6%26sSearch_6%3D%26bRegex_6%3Dfalse%26bSearchable_6%3Dtrue%26bSortable_6%3Dtrue%26mDataProp_7%3D7%26sSearch_7%3D%26bRegex_7%3Dfalse%26bSearchable_7%3Dtrue%26bSortable_7%3Dtrue%26mDataProp_8%3D8%26sSearch_8%3D%26bRegex_8%3Dfalse%26bSearchable_8%3Dtrue%26bSortable_8%3Dtrue%26mDataProp_9%3D9%26sSearch_9%3D%26bRegex_9%3Dfalse%26bSearchable_9%3Dtrue%26bSortable_9%3Dtrue%26mDataProp_10%3D10%26sSearch_10%3D%26bRegex_10%3Dfalse%26bSearchable_10%3Dtrue%26bSortable_10%3Dtrue%26mDataProp_11%3D11%26sSearch_11%3D%26bRegex_11%3Dfalse%26bSearchable_11%3Dtrue%26bSortable_11%3Dtrue%26mDataProp_12%3D12%26sSearch_12%3D%26bRegex_12%3Dfalse%26bSearchable_12%3Dtrue%26bSortable_12%3Dtrue%26mDataProp_13%3D13%26sSearch_13%3D%26bRegex_13%3Dfalse%26bSearchable_13%3Dtrue%26bSortable_13%3Dtrue%26mDataProp_14%3D14%26sSearch_14%3D%26bRegex_14%3Dfalse%26bSearchable_14%3Dtrue%26bSortable_14%3Dtrue%26mDataProp_15%3D15%26sSearch_15%3D%26bRegex_15%3Dfalse%26bSearchable_15%3Dtrue%26bSortable_15%3Dtrue%26mDataProp_16%3D16%26sSearch_16%3D%26bRegex_16%3Dfalse%26bSearchable_16%3Dtrue%26bSortable_16%3Dtrue%26mDataProp_17%3D17%26sSearch_17%3D%26bRegex_17%3Dfalse%26bSearchable_17%3Dtrue%26bSortable_17%3Dtrue%26mDataProp_18%3D18%26sSearch_18%3D%26bRegex_18%3Dfalse%26bSearchable_18%3Dtrue%26bSortable_18%3Dtrue%26mDataProp_19%3D19%26sSearch_19%3D%26bRegex_19%3Dfalse%26bSearchable_19%3Dtrue%26bSortable_19%3Dtrue%26sSearch%3D%26bRegex%3Dfalse%26iSortCol_0%3D0%26sSortDir_0%3Dasc%26iSortingCols%3D1%26export%3Dcsv) (after login, this should download as a CSV; `CANCER_GENE_CENSUS`)
* [OncoKB cancer gene list](https://oncokb.org/api/v1/utils/cancerGeneList.txt) (`ONCOKB_GENES`)
* Patient somatic mutation call files: “Mutation_Packager_Oncotated_Calls” tar.gz for each cancer type from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ (!!)
* UniProt splicing-associated gene file (TSV): keyword:"mRNA splicing [KW-0508]" & Reviewed:yes & organism:"Homo sapiens (Human) [9606]" from https://www.uniprot.org/uniprot/

Note we also made use of [metaSRA](http://metasra.biostat.wisc.edu/)'s ontology terms obtained via its API and available in https://github.com/JulianneDavid/shared-cancer-splicing/tree/master/SRA_junction_download/MetaSRA_run_files .

## Execution

1. Create TCGA and GTEx junction index by running `jx_indexer.py` in `index` mode:

        python3 jx_indexer.py -d DB_DIR index -c GTEx_JUNCTION_COVERAGE -C TCGA_JUNCTION_COVERAGE -b GTEx_JUNCTION_BED -B TCGA_JUNCTION_BED -p GTEx_PHEN -P TCGA_PHEN -s RECOUNT_SAMPLE_IDS -g GENCODE_ANNOTATION_GTF

2. Run `jx_indexer.py` in `experiment` mode on junction index to collect junction files for analysis:

        python3 jx_indexer.py -d DB_DIR experiment -o OUTPUT_DIR

`OUTPUT_DIR` will now contain:
* `all_jxs` (`ALL_JX_DIR`) for use in set membership annotation;
* `non-core-normal_jxs_per_sample` (`NON_CORE_NORMAL_DIR`) and `non-tissue-matched_jxs_per_sample` (`NON_TISSUE_MATCHED_NORMAL_DIR`) for use in set membership annotation and Figures S3A, S3B, S4A, and S4B;
* `non-core-normal_counts_per_sample` (`COUNTS_PER_SAMPLE_DIR`) (junctions not found in core normals) to generate Figure 1B;
* `non-core-normal_jx_prevalences` (`PREVALENCE_FILE_DIR`) (TCGA cancer junctions not found in core normals with cancer-type prevalence, per cancer type) to generate Figures 2A, 2B, 2C, and S5;
* `FigS1a_non-core-normal_counts_per_sample` (`FILTERED_NTM_JX_PER_SAMPLE_DIR`) (coverage- and annotation-filtered junctions not found in GTEx or TCGA tissue-matched normal samples) and `FigS1a_non-tissue-matched_counts_per_sample` (`FILTERED_NCN_JX_PER_SAMPLE_DIR`) (coverage- and annotation-filtered junctions not found in core normals) to generate Figure S1A;
* (`all_jxs_per_sample_paired_normals`) `ALL_JXS_PER_SAMPLE_DIR` containing all junctions for TCGA SKCM normal samples, TCGA SKCM tumor samples, and GTEx bulk skin normal samples to generate Figure S7.

3. Call [this directory](SRA_junction_download/MetaSRA_run_files) `METASRA_QUERY_FILES`. Run:

        python3 create_snaptron_queries.py -s METASRA_QUERY_FILES -i RECOUNT_SAMPLE_IDS -o LOGGING_OUTPUT_DIR -p SRA_JUNCTION_OUTPUT_DIR -u UTILITIES_DIR

where `UTILITIES_DIRECTORY` is where `snaptron-experiments` was previously cloned (see Requirements). `SRA_JUNCTION_OUTPUT_DIR` is where the results will appear.

4. Run each `SAMPLE_TYPE_indiv_chrom_batch_query_script.sh` to collect SRA junctions.

`SRA_JUNCTION_OUTPUT_DIR` will now contain:
* Non-cancer junction results will be saved in 
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




