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

### Preliminaries
1. Create TCGA and GTEx junction index by running `jx_indexer.py` in `index` mode:

        python3 jx_indexer.py -d DB_DIR index -c GTEx_JUNCTION_COVERAGE -C TCGA_JUNCTION_COVERAGE -b GTEx_JUNCTION_BED -B TCGA_JUNCTION_BED -p GTEx_PHEN -P TCGA_PHEN -s RECOUNT_SAMPLE_IDS -g GENCODE_ANNOTATION_GTF

2. Run `jx_indexer.py` in `experiment` mode on junction index to collect junction files for analysis:

        python3 jx_indexer.py -d DB_DIR experiment -o OUTPUT_DIR

`OUTPUT_DIR` will now contain:
* `all_jxs` (`ALL_JX_DIR`) for use in set membership annotation;
* `non-core-normal_jxs_per_sample` (`NON_CORE_NORMAL_DIR`) and `non-tissue-matched_jxs_per_sample` (`NON_TISSUE_MATCHED_NORMAL_DIR`) for use in set membership annotation;
* `non-core-normal_counts_per_sample` (`COUNTS_PER_SAMPLE_DIR`) (junctions not found in core normals) to generate Figure 1B;
* `non-core-normal_jx_prevalences` (`PREVALENCE_FILE_DIR`) (TCGA cancer junctions not found in core normals with cancer-type prevalence, per cancer type) to generate Figures 2A, 2B, 2C, and S2A;
* `FigS1a_non-core-normal_counts_per_sample` (`FILTERED_NCN_JX_PER_SAMPLE_DIR`) (coverage- and annotation-filtered junctions not found in GTEx or TCGA tissue-matched normal samples) and `FigS1a_non-tissue-matched_counts_per_sample` (`FILTERED_NTM_JX_PER_SAMPLE_DIR`) (coverage- and annotation-filtered junctions not found in core normals) to generate Figure S1A;
* `jx_non-core-normal_jxs_per_sample` (`NON_CORE_NORMAL_JX_COORD_DIR`) for use in Figures S1E, S1F, S1G, and S1H;
* (`all_jxs_per_sample_paired_normals`) `ALL_JXS_PER_SAMPLE_DIR` containing all junctions for TCGA SKCM normal samples, TCGA SKCM tumor samples, and GTEx bulk skin normal samples to generate Figure S2B.

3. Call [this directory](SRA_junction_download/MetaSRA_run_files) `METASRA_QUERY_FILES`. Run:

        python3 create_snaptron_queries.py -s METASRA_QUERY_FILES -i RECOUNT_SAMPLE_IDS -o LOGGING_OUTPUT_DIR -p SRA_JUNCTION_OUTPUT_DIR -u UTILITIES_DIR

where `UTILITIES_DIRECTORY` is where `snaptron-experiments` was previously cloned (see Requirements). `SRA_JUNCTION_OUTPUT_DIR` is where the results will appear.

4. Run each `SAMPLE_TYPE_indiv_chrom_batch_query_script.sh` in `SRA_JUNCTION_OUTPUT_DIR` to collect SRA junctions.

`SRA_JUNCTION_OUTPUT_DIR` will now contain:
* `SRA_noncancer_rawresults` (`SNAPTRON_NONCANCER_DIR`) for non-cancer junction results
* `SRA_noncancer_exptlists` (`SNAPTRON_NONCANCER_EXPTLIST_DIR`) for non-cancer experiment lists
* `SRA_cancer_rawresults` (`SNAPTRON_CANCER_DIR`) for cancer junction results
* `SRA_cancer_exptlists` (`SNAPTRON_CANCER_EXPTLIST_DIR`) for cancer experiment lists 
        
5. Run `collect_1-read_TCGA_jxs.py` to collect single-read TCGA junction collection:

        python3 collect_1-read_TCGA_jxs.py -C TCGA_JUNCTION_COVERAGE -B TCGA_JUNCTION_BED -o JSON_DIR

`JSON_DIR` will now contain:
* `single-read-tcga-jxs_json.txt` containing single-read junctions for use in set membership annotation runs below, using `set_membership_annotation.py`, as follows:

        python3 set_membership_annotation.py --db-path DB_DIR --snaptron-results SNAPTRON_NONCANCER_DIR -d ALL_JX_DIR -g NON_CORE_NORMAL_DIR -p NON_TISSUE_MATCHED_NORMAL_DIR --gtf-file GENCODE_ANNOTATION_GTF --single-read-jx-json JSON_DIR/single-read-tcga-jxs_json.txt --cancer-sra-directory SNAPTRON_CANCER_DIR --cancer-gene-census CANCER_GENE_CENSUS --oncokb-cancer-genes ONCOKB_GENES --min-overall-set-count 1 -o FULL_PIECHART_DIR

`FULL_PIECHART_DIR` will now contain:
* `unexplained` (`UNEXPLAINED_DIR`) containing unexplained subset of all single-read junctions
* `developmental` (`DEVELOPMENTAL_DIR`) containing developmentally-occurring subset of all single-read junctions

To generate two-sample minimum junction sets (used in Figures S3E and S3F), run:

        python3 set_membership_annotation.py --db-path DB_DIR --snaptron-results SNAPTRON_NONCANCER_DIR -d ALL_JX_DIR -g NON_CORE_NORMAL_DIR -p NON_TISSUE_MATCHED_NORMAL_DIR --gtf-file GENCODE_ANNOTATION_GTF --single-read-jx-json 1_READ_TCGA_JX_JSON --cancer-sra-directory SNAPTRON_CANCER_DIR --cancer-gene-census CANCER_GENE_CENSUS --oncokb-cancer-genes ONCOKB_GENES --min-overall-set-count 2 -o 2-SAMPLE_PIECHART_DIR

`2-SAMPLE_PIECHART_DIR` will now contain:
* `unexplained` (`UNEXPLAINED_2-SAMPLE_DIR`) containing unexplained subset of all two-sample junctions
* `developmental` (`DEVELOPMENTAL_2-SAMPLE_DIR`) containing developmentally-occurring subset of all two-sample junctions

6. To perform additional set membership analyses using `set_membership_analysis.py` (for use in Figures S1B and S1C), run:

        python3 set_membership_analysis.py -o OUTPUT_DIR -s FULL_PIECHART_DIR

`FULL_PIECHART_DIR` will now contain:
* `unexplained` (`UNEXPLAINED_DIR`) containing unexplained subset of all single-read junctions
* `developmental` (`DEVELOPMENTAL_DIR`) containing developmentally-occurring subset of all single-read junctions
* `true_TCGA_prevalence_files` (`TCGA_PREVALENCE_DIR`) containing non-core-normal junctions

To tally SRA experiments using `count_unique_SRA_expts.py`, run:

        python3 count_unique_SRA_expts.py -c SNAPTRON_CANCER_EXPTLIST_DIR -n SNAPTRON_NONCANCER_EXPTLIST_DIR -o FIGURE_OUTPUT_DIR

### Figure generation

Fig 1A:
        
        python3 fig1A_overall_set_barplots.py -s FULL_PIECHART_DIR -o FIGURE_OUTPUT_DIR

Fig 1B:

        python3 fig1B_ncn_jx_counts_per_sample.py -j COUNTS_PER_SAMPLE_DIRECTORY -o FIGURE_OUTPUT_DIRECTORY

Fig 1C:

        python3 fig1C_overall_set_prevalence_boxplot.py --set-memberships FULL_PIECHART_DIR -o FIGURE_OUTPUT_DIR

Fig 2A:

        python3 fig2A_TCGA_heatmap.py -d PREVALENCE_FILE_DIR -o FIGURE_OUTPUT_DIR

Fig 2B:

        python3 fig2B_TCGA_subtypes_heatmap.py -d PREVALENCE_FILE_DIR -o FIGURE_OUTPUT_DIR

Fig 2C:

        python3 fig2C_cell_of_origin_heatmap.py  -d PREVALENCE_FILE_DIR -o FIGURE_OUTPUT_DIR --snaptron-results SNAPTRON_NONCANCER_DIR -e SNAPTRON_NONCANCER_EXPTLIST_DIR

Figs 3A, S3A, S3E, and S3F:

(Contact on instructions for regenerating this figure: https://github.com/metamaden)

Instantiate the functions in `upset_functions.R`, and make a new subdirectory (default name: `cxdat`), and run `get.datlist()`, specifying the name of the subdirectory containing the junction data tables. Once this has run, the working directory should contain a large list object containing the junction data tables, as well as a newly populated subdirectory (default: `cxdat`). Depending on your session memory limit, you may need to remove the object `cxdl` from your environment before proceeding.

Follow the remaining steps in `upset_data.R` to: Define the new `developmental` junction category; Run `jx.calcsets()` for the supplemental data, then for the main figure data; Replace the `unexplained` category in the main figure data with the supplemental figure data; Generate the supplemental data tables. To save memory, tables are read processively and processed.

Next, follow the steps in `upset_plots_original.R` to: Specify the data table containing the whole/total group sets and subset group sets, respectively; Exclude the pre-filtered set columns from visualization; Properly order whole set columns by mean value; Automatically generate the axis limits and labels; Specify the group labels (should match whole set column labels after ordering on means); Generate the plot data objects; Generate the final plot images.

Following these steps should generate both a main and a supplemental upset plot with annotation set and subset abundances. It may be necessary to experiment with different plot and image dimensions in `upset_plots_original.R` and `make.ggset.final()` to refine the figure properties. For additional information, refer to the function docstrings and script comments.)

Fig 3B:

        python3 fig3B_antisense_boxplot.py -o FIGURE_OUTPUT_DIR -s FULL_PIECHART_DIR

Fig S1A:

        python3 fig1B_S1A_ncn_jx_counts_per_sample.py -j FILTERED_NCN_JX_PER_SAMPLE_DIR -p FILTERED_NTM_JX_PER_SAMPLE_DIR -d -g THYM CESC UVM DLBC --prepared-sort-order -o FIGURE_OUTPUT_DIR

Figs S1B and S1C:

        python3 figS1BC_S3CD_junction_sharedness.py -d TCGA_PREVALENCE_DIR -o FIGURE_OUTPUT_DIR

Note: requires contents of `TCGA_PREVALENCE_DIR` created by running set_membership_analysis.py in step 6 above.

Fig S1D:

        python3 figS1D_set_prevalences_per_cancer.py -s FULL_PIECHART_DIR -o FIGURE_OUTPUT_DIR

Figs S1E, S1F, S1G, and S1H:

(Contact for instructions on regenerating this figure: https://github.com/weederb23)

1. Download the TCGA mutation files from http://gdac.broadinstitute.org/runs/stddata__2016_01_28/data/ for each cancer type and move into the sub-directory ./TCGA_mut_files/, where ./ is the location of "SF_mutation_script.R" and "SF_sharedness_plots.R"

2. Copy files from NON_CORE_NORMAL_JX_COORD_DIR to a new directory within ./ called "non_core_jxns/". Rename files from full cancer names to match TCGA cancer codes.

3. Copy files from [PIECHART DIRECTORY] to a new directory within ./ called "non_core_jxn_annotations/". Rename files from full cancer names to match TCGA cancer codes.

4. Run `Rscript SF_mutation_script.R`.

This generates `fig_s1e.jpg`, `fig_s1f.jpg` and `shared_jxn_df.txt`. `shared_jxn_df.txt` is used downstream for `SF_sharedness_plots.R`.

If run interactively, also presents statistics for junction comparisons between patients with and without splicing factor mutations across cancers.

5. Run `SF_sharedness_plots.R`.

This generates `fig_s1g.jpg` and `fig_s1h.jpg`. If run interactively, also presents statistics for jxn comparisons between patients with and without splicing factor mutations across cancers.

Fig S1I:

        python3 figS1I_S3B_SRA_cancer_shared_prevalence.py --snaptron-results SNAPTRON_CANCER_DIR -e SNAPTRON_CANCER_EXPTLIST_DIR -o FIGURE_OUTPUT_DIR -d PREVALENCE_FILE_DIR

Fig S2A:

        python3 figS2A_full_TCGA_SRA_heatmap.py  -d PREVALENCE_FILE_DIR -o FIGURE_OUTPUT_DIR --snaptron-results SNAPTRON_NONCANCER_DIR -e SNAPTRON_NONCANCER_EXPTLIST_DIR

Fig S2B:

        python3 figS2B_SKCM_jx_similarity.py --snaptron-results SNAPTRON_NONCANCER_DIR -d ALL_JXS_PER_SAMPLE_DIR -o FIGURE_OUTPUT_DIR

Fig S3B:

        python3 figS1I_S3B_SRA_cancer_shared_prevalence.py --snaptron-results SNAPTRON_CANCER_DIR -e SNAPTRON_CANCER_EXPTLIST_DIR -o FIGURE_OUTPUT_DIR -d UNEXPLAINED_PIECHART_DIR --unexplained-junctions

Fig S3C and S3D:

        python3 figS1BC_S3CD_junction_sharedness.py -d UNEXPLAINED_PIECHART_DIR -o FIGURE_OUTPUT_DIR
