#impute function

# FUNC_data = tibble containing data
# FUNC_metabolite_list = character array of metabolites (to extract from data tibble for impute)
# FUNC_option_impute_missing_data = TRUE/FALSE option to impute

lgw_lipid_impute <- function(FUNC_data,
                             FUNC_metabolite_list,
                             FUNC_option_impute_missing_data)
  {

  FUNC_list <- list()
  FUNC_list$data <- FUNC_data
  FUNC_list$metabolite_data <- FUNC_list$data %>% select(all_of(FUNC_metabolite_list))
#step 1 = replace all NA with 0 (missing value) prior to impute
  FUNC_list$metabolite_data[is.na(FUNC_list$metabolite_data)] <- 0

# #find missing value % in current data frame
# temp_frame <- master_list$impute_data[[idx_data]] %>% select(!contains("sample")) 

# total_values <- dim(temp_frame)[1] * dim(temp_frame)[2]

total_missing_values <- length(which(FUNC_list$metabolite_data==0))
  
if(total_missing_values != 0 & FUNC_option_impute_missing_data == TRUE){
  #find idx of columns containing a 0 value:
  zero_value_column_idx <- which(
    colSums(
      FUNC_list$metabolite_data == 0) > 0) %>% 
    as.matrix() %>% 
    c()

  #min value impute
  for(idx_num_impute in zero_value_column_idx){
    temp_impute_data <- FUNC_list$metabolite_data[,idx_num_impute]
    non_zero_values_idx <- which(temp_impute_data > 0) 
    zero_values_idx <- which(temp_impute_data == 0)
    if(length(non_zero_values_idx) > 0){
    impute_value <- min(temp_impute_data[non_zero_values_idx,])/2
    FUNC_list$metabolite_data[zero_values_idx,idx_num_impute] <- impute_value
    }
  }
  
  }

FUNC_list$data_out <- bind_cols(
  FUNC_list$data %>% select(contains("sample")),
  FUNC_list$metabolite_data
)

FUNC_list$data_out

}