#' ---
#' title: "Targeted Lipid ExploreR QC Report"
#' author: ANPC
#' output: html_document
#' 
#' ---
#'
#' ***
#' #### Project summary
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$project_summary)
#' 
#' 
#' ***
#' ### Data import summary
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$transposed_summary)
#'
#'
#' ***
#' #### PCA: Raw skyline imports
#' PCA displaying raw data imported from skyline. 
#' * There has been no outlier removal or data processing at this point.
#'
#'
#'
#'
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
subplot(master_list$pca_output$data_sorted$sample_qc$plot_scores, master_list$pca_output$data_sorted$plate$plot_scores)
#'
#'
#' ***
#' 
#' ### Process step: Missing value filter
#' Missing value filter: 
#' * Step 1: Remove all samples that have > 50% missing values (removes any mis-injections etc that may be present in the data)
#' * Step 2: Remove all metabolite features that have > 50% missing values (zero, NA, NaN etc)
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$missing_value_filter_summary)
#'
#'
#' The following internal standards failed the missing value filter and were removed from the project. 
#' Note - metabolite targets that use these internal standards for calculation of response ratio and concentrations are also removed.
#'
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_lists$missing_value_filter$failed_SIL %>% as_tibble() %>% rename(internal_standard = value))
#' 
#' ***
#' 
#' #### PCA: Post-missing value filter
#' PCA displaying data that has undergone:
#' * missing value filtering
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
subplot(master_list$pca_output$missing_value_filter$sample_qc$plot_scores, master_list$pca_output$missing_value_filter$plate$plot_scores)
#' 
#' 
#' ***
#' 
#' ### Process step: Imputation 
#' * Imputation of the remaining zero value and missing data 
#' * Imputation is completed using x/2, where x is minimum intensity of that feature in the batch 
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$impute_table)
#'
#'
#' ***
#' 
#' ### Process step: Response ratio and concentration value 
#' Two step process:
#' 
#' * Calculation of target metabolite/stable isotope labelled (SIL) internal standard ratio, using predefined target metabolite/internal standard pairs
#' * Conversion of response ratio to concentration values using single point calibration
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$concentration_summary)
#' 
#' 
#' ***
#' 
#' ##### PCA: Post-calculation of internal standard ratio's and conversion to concentration values
#' PCA displaying data that has undergone:
#' * missing value filtering
#' * imputation of remaining missing values
#' *calculation of response ratios (target analyte peak area/internal standard peak area) 
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
subplot(master_list$pca_output$concentration$sample_qc$plot_scores, master_list$pca_output$concentration$plate$plot_scores)
#'
#'
#' ***
#' 
#' ### Process step: Sample outlier filter
#' Filter to remove all outlier samples with excessive principal component (PC) variation
#' * Step 1 : Create PCA scores for PC 1:3
#' * Step 2: Find samples with PC score > 1.5 standard deviation of median PC
#' 
#' note - filter is in development and under evaluation
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$pc_filter_summary)
#'
#' ***
#' #### PCA filter plot: Run order vs PC score
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pc_filter$PC1
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pc_filter$PC2
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pc_filter$PC3
#'
#' ***
#' 
#' #### PCA: Post-principal component filter
#' PCA displaying data that has undergone:
#' * missing value filtering
#' * imputation of remaining missing values
#' * calculation of response ratios (target analyte peak area/internal standard peak area) 
#' * principal component filtration
#' 
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
subplot(master_list$pca_output$pc_filter$sample_qc$plot_scores, master_list$pca_output$pc_filter$plate$plot_scores)
#'
#'
#' ***
#'
#' ### Process step: Signal drift and batch correct the data (per project)
#' * Data from each individual batch undergoes signal drift correction using statTarget package (https://stattarget.github.io/)
#' * This is performed within individual batches at this point to evaluate the performance of each batch
#' 
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
knitr::kable(master_list$summary_tables$batch_correction_overview)
#' 
#'
#' ***
#' 
#' #' #### Run order plot: pre-statTarget correction
#' Plot of run order vs principle components
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pre_statTarget$PC1$plotly
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pre_statTarget$PC2$plotly
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$pre_statTarget$PC3$plotly
#' 
#' ***
#' 
#' #### Run order plot: post-statTarget correction
#' Plot of run order vs principle components
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$post_statTarget$PC1$plotly
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$post_statTarget$PC2$plotly
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
master_list$pc_run_plot$post_statTarget$PC3$plotly
#'
#' ***
#' 
#' #### PCA: Post-signal drift/batch correction
#' PCA displaying data that has undergone:
#' * missing value filtering
#' * imputation of remaining missing values
#' * calculation of response ratios (target analyte peak area/internal standard peak area) 
#' * principal component filtration
#' * signal drift/batch correction using statTarget package
#' 
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
subplot(master_list$pca_output$statTarget_corrected$sample_qc$plot_scores, master_list$pca_output$statTarget_corrected$plate$plot_scores)
#'
#' 
#' ***
#' #### Environment summary
#' R version
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
print(master_list$environment$r_version)
#' Base packages
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
print(master_list$environment$base_packages)
#' User packages
#+ echo=FALSE, message=FALSE, fig.width=10, fig.height=5
print(master_list$environment$user_packages)
#' ***





