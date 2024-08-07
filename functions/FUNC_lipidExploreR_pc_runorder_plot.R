#ANPC PCA quality control visualisation

# FUNC_data = data containing individual lipid data - MUST CONTAIN sampleID column 
# FUNC_colour_by = how to colour the plot (e.g. sample class, is_ltr or cohort)
# FUNC_plot label = what to label the scores plot with (e.g. sampleID)
# FUNC_scaling = UV or Pareto

lgw_pc_run_plot <- function(FUNC_data, 
                    FUNC_metabolite_list, 
                    FUNC_colour_by, 
                    FUNC_plot_label, 
                    FUNC_scaling,
                    FUNC_title,
                    FUNC_project_colours,
                    FUNC_option_point_size,
                    FUNC_option_plot_qc
                    ){
  
  pca_output <- list()
  pca_output$plots <- list()
  
  title_text <- FUNC_title
  
  qc_idx <- which(FUNC_data[["sample_type"]] == "qc")
  
  if(FUNC_option_plot_qc == FALSE){
    FUNC_data <- FUNC_data %>% 
      filter(sample_type != "qc")
  }
  

  #create data matrix for PCA
  pca_x <- FUNC_data %>%  select(all_of(FUNC_metabolite_list)) %>% as.matrix() 
  pca_x[pca_x == 0] <- NA #remove all 0 values (above adds 1 to all values therefore anything that = 1 was a 0)
  pca_x[is.infinite(pca_x)] <- NA #remove all infinite values
  min_value <- min(pca_x, na.rm = TRUE) # find the lowest value in the matrix
  pca_x[is.na(pca_x)] <- min_value/100 # replace all NA, Inf, and 0 values with the lowest value in the matrix/100 to represent a value below limit of detection
  
  pca_x <- log(pca_x+1) #log values for plotting
  
  #create PCA model
  pca_output$pca_model <- pca(pca_x, pc=3, scale = paste(FUNC_scaling), center = TRUE)
  
  # extract score values for plotting in plot_ly
  PC1 <- as.numeric(as.matrix(pca_output$pca_model@t[,1]))
  PC2 <- as.numeric(as.matrix(pca_output$pca_model@t[,2]))
  PC3 <- as.numeric(as.matrix(pca_output$pca_model@t[,3]))

  plot_Val <- as_tibble(cbind(PC1, PC2, PC3, FUNC_data$sample_idx, FUNC_data$sample_name,  FUNC_data$sample_type, FUNC_data$sample_plate_id)) %>% 
    setNames(c("PC1", "PC2", "PC3", "sample_idx", "sample_name", "sample_type", "sample_plate_id"))
  #add file_name to sample_idx for plotly plot label
  plot_Val$sample_idx <- c(1:nrow(plot_Val)) %>% as.numeric()
  #create factor to control plot order
  plot_Val$sample_idx <- paste0(plot_Val$sample_idx, "_", plot_Val$sample_name) %>%
    factor(ordered = TRUE, levels = paste0(plot_Val$sample_idx, "_", plot_Val$sample_name))
  plot_Val$PC1 <- plot_Val$PC1 %>% as.numeric()
  plot_Val$PC2 <- plot_Val$PC2 %>% as.numeric()
  plot_Val$PC3 <- plot_Val$PC3 %>% as.numeric()

  #set object to store
  pca_output$PC1 <- list()
 
  #set object to store
  pca_output$PC2 <- list()
  
  #set object to store
  pca_output$PC3 <- list()
  
  # set plot attributes (controlled by FUNC_colour_by and FUNC_plot_label)
  pca_colour <- list()
  pca_colour <- FUNC_data %>% select(all_of(FUNC_colour_by)) #%>% as.matrix()
  colnames(pca_colour) <- "pca_colour" 
 
  pca_plot_colour <- pca_colour$pca_colour
  #pca_plot_colour[is.na(pca_plot_colour)] <- "none"
  
  #set colours
  plot_colours <- c(FUNC_project_colours)

 
  ##################  ##################  ##################  ##################
  #produce run order plot
  ##################  ##################  ##################  ##################
  
  for(idx_PC in plot_Val %>% select(-contains("sample")) %>% names()){
  
  bp <- ggplot(data=plot_Val,
               aes(x=sample_idx,
                   y=get(idx_PC))
  )
  
  bp <- bp + geom_point(aes(fill = pca_plot_colour#,
  ),
  shape = 21,
  size = FUNC_option_point_size
  )
  
  bp <- bp + scale_fill_manual(values = c(plot_colours))
  
  bp <- bp + labs(x = paste("Sample order"),
                  y = paste0(idx_PC))
  bp <- bp + ggtitle(paste0(idx_PC))
  bp <- bp + theme_cowplot() 
  bp <- bp + theme(
    plot.title = element_text(hjust = 0.5, size=14),
    axis.text.y = element_text(size = 12, margin = margin(t = 0, r = 0, b = 0, l = 2)),
    #axis.text.x = element_text(size = 12, angle = 45, vjust = 1, hjust = 1),
    axis.text.x = element_blank(),
    axis.title = element_text(size = 14),
    #legend.text=element_text(size=12),
    legend.title=element_blank(),
    axis.ticks.x = element_blank(),
    
  )
  
  #create vertical lines to separate classes on plot
  FUNC_data_batches <- plot_Val$sample_plate_id %>% unique()
  batch_idx <- NULL
  
  if(length(FUNC_data_batches) > 1){
   for(idx_batch in FUNC_data_batches[2:length(FUNC_data_batches)]){
     batch_idx <- c(batch_idx, min(which(plot_Val$sample_plate_id == idx_batch)))
   }
    
    bp <- bp + geom_vline(xintercept=c(batch_idx),color="grey")
    }

  pca_output$plots[[idx_PC]]$plotly <- bp %>% ggplotly() %>% layout(legend = list(orientation = "h",   # show entries horizontally
                                                                                      xanchor = "center",  # use center of legend as anchor
                                                                                      x = 0.5,
                                                                                      y = -0.2,
                                                                                      title=list(text="Group comparison:"))
  )
     
  
 
  
  
  
  
  
  
   }
 
  pca_output$plots
}
