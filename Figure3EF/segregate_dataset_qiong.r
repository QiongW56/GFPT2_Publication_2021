# A function for extracting the top/bottom n% of a 
# gene expression matrix based on the expression of 'list_of_genes'

segregate_dataset <- function(dataset, # gene expression datasets (e.g. RNAseq)
                              genes, # gene list or one gene
                              nth_percentage, # 0.01 or 0.25 etc.
                              cor_data2, # survical data of patients
                              type_of_effect, # the effects ("Beneficial", "Harmful" or "No effect")
                              print_plot = TRUE){ # if it is TRUE, then effects cannot be empty
  
  if (missing(print_plot))
  {
    print_plot <- TRUE
  }
  
  # If you want to plot the KM plot (print_plot == TRUE), "type_of_effect" should not be empty
  if (print_plot == TRUE){ 
    id_no_effect <- which(type_of_effect == 'No effect!')
  } else {
    id_no_effect <- c()
  }
  
  # If you want to plot, "type_of_effect" should not be empty, Delete the "No effect!" from 
  # the "type_of_effect".
  if (class(genes) != 'data.frame'){ # set up the targeted gene list into data frame
    
    genes <- data.frame(genes)
    
    if(length(id_no_effect) != 0){
      
      type_of_effect <- type_of_effect[-id_no_effect]
      genes <- data.frame(genes[-id_no_effect,]) # Delete the "No effect" genes
      colnames(genes) <- 'Group1'
      
    } else {
      
      type_of_effect <- type_of_effect
      genes <- genes
      colnames(genes) <- 'Group1'
      
    }
    
  } else {
   
    if(length(id_no_effect) != 0){
      
      type_of_effect <- type_of_effect[-id_no_effect]
      genes <- genes[-id_no_effect, ] # Delete the "No effect" genes
      colnames(genes) <- 'Group1'
      
    } else {
      
      type_of_effect <- type_of_effect
      genes <- genes
      colnames(genes) <- 'Group1'
      
    }
    
  }
  
  # at least some of the targeted list of genes should be within the gene expression dataset
  # print(genes[,1])
  id_used_3 <- which(rownames(dataset) %in% genes[, 1])
  if (length(id_used_3) == 0){
    print('No gene with any significant effect on survival, this will result in an error..')
  }
  
  # select gene expression dataset with targeted list of genes
  # print(id_used_3)
  dataset2 <- dataset[id_used_3, ]
  
  # z-score the gene expression data (package: "siggitRausti")
  #print(dataset2)
  zdata <- standardize(dataset2)
  
  # ADDED 090519 as a test-------------------
  if (print_plot == TRUE){
    for (i in 1:dim(dataset2)[1]){
     
      id_gene <- which(genes$Group1 %in% rownames(dataset2)[i])
      effect_type <- type_of_effect[id_gene]
      
      if (length(effect_type == 'Beneficial!') == 0){
        
        zdata[i, ] = zdata[i, ] # "Harmful"
        
      } else if (length(effect_type == 'Beneficial!') == 1){
        
        if (effect_type == 'Beneficial!'){
          zdata[i, ] = -as.numeric(as.character(zdata[i, ])) # "Beneficial"
        } else {
          zdata[i, ] = zdata[i, ]
        }
        
      }
    }
  }
  
  # test if all values are the same (no sd = NA values in standardized data)
  if (all(is.na(c(NA, zdata)))){ # dataframe to list
    
    if (print_plot == T){
      
      p3 <- 'Standard deviation is 0, no standardization possible!'
      
    } else if (print_plot == F){
      
      p3 <- 1
      
    }
    
  } else {
    # This could use some adjustment, 1-2 genes could be carrying the rest.
    sum_genes = colSums(zdata) 
    
    # top
    thresh_val_top = quantile(sum_genes,
                              1-nth_percentage)
    patients_top <- data.frame(names(sum_genes)[which(sum_genes >= thresh_val_top)])
    colnames(patients_top) = 'patients'
    
    # bottom
    thresh_val_bottom = quantile(sum_genes,
                                 nth_percentage)
    patients_bottom <- data.frame(names(sum_genes)[which(sum_genes <= thresh_val_bottom)])
    colnames(patients_bottom) = 'patients'
    
    dataset_top <- dataset[,which(colnames(dataset) %in% patients_top$patients)]
    dataset_bottom <- dataset[,which(colnames(dataset) %in% patients_bottom$patients)]
    
    # Here I need to assert whether cor_data2 contains a 'patient','survival' and 'status' vectors. Otherwise, it will
    # result in a boring error that is hard to track unless going directly through all lines in the function...
   
    # top
    df_patients_top <- cor_data2[which(cor_data2$patient %in% colnames(dataset_top)), ]
    df_patients_top <- cbind(df_patients_top,rep('top',
                                                 nrow(df_patients_top))) # add column "top"
    colnames(df_patients_top)[ncol(df_patients_top)] = 'top_or_bottom' # change column name
    
    # bottom
    df_patients_bottom <- cor_data2[which(cor_data2$patient %in% colnames(dataset_bottom)), ]
    df_patients_bottom <- cbind(df_patients_bottom,rep('bottom', 
                                                       nrow(df_patients_bottom))) # add column "bottom"
    colnames(df_patients_bottom)[ncol(df_patients_bottom)] = 'top_or_bottom' # change column name
   
    # row combine
    df_fin = rbind(df_patients_top, 
                    df_patients_bottom)
    
    # Now do the survival curve:
    # package: "survival"
    fit_chosen <- survfit(Surv(survival, status) ~ top_or_bottom, data = df_fin) 
    
    if (print_plot == TRUE){
      
      if (dim(dataset2)[1] == 0) { # Checking if the dataframe is empty (if checking one gene's effect on survival)
        print('No genes left in dataframe!')
      }
      
      # package: "survminer"
      p3 <- ggsurvplot(fit_chosen, 
                       data = df_fin, 
                       risk.table = FALSE,
                       pval = TRUE, 
                       conf.int = TRUE, # confident intervial
                       xlim = c(0, 60), # 5 years
                       xlab = 'Time in months', # Choose length of follow-up time (120 = 10 years, 60 = 5 years)
                       break.time.by = 12, # breaks for x axis (every year)
                       legend.title = "Effective Gene Expression Profile",
                       legend.labs = c('Beneficial Profile','Harmful Profile'), # change legends
                       legend =  c(0.78, 0.4), # "right"
                       # palette = c("gray72", "gray24"),
                       palette = c('steelblue','red'), # change line colors
                       ggtheme = theme_classic2(base_size = 24),
                       surv.scale="percent")
      
      # package: "ggpubr"
      p3 <- ggpar(p3,
                  font.x = c(20, "bold", "black"), # x axis title
                  font.y = c(20, "bold", "black"), # y axis title
                  font.tickslab = c(18, "bold", "black"),
                  font.legend = c(12, "bold", "black")) # x and y axis
      
    } else {
      
      if (dim(dataset2)[1] == 0) { # Checking if the dataframe is empty (if checking one gene's effect on survival)
        
        print('No genes left in dataframe!')
        p3 <- NA
        
      } else {
        # package: "survival"
        res.cox_cand <- coxph(Surv(survival, status) ~ top_or_bottom, data = df_fin)
        hehe <- summary(res.cox_cand)
        p3 <- hehe$sctest[3] # This is the p-value - return that!
        
      }
      
    }
    
  }
  
  return(p3)
  
}
