# shared-cancer-splicing
Code for reproducing analyses and figures for shared alternative cancer splicing paper

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




