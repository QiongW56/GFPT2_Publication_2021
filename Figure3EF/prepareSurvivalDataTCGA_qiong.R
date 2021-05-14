function (dataset, # gene expression dataset
          genes,  # list of genes
          quantile_used, # one element from the quantile vector
          percentage_patients)  # e.g. 0.6
{
  # colnames of the gene list
  if (colnames(genes) == "Group1") {
    
    id_used_3 <- which(rownames(dataset) %in% genes$Group1)
    
  } else if (colnames(genes) == "Group2") {
    
    id_used_3 <- which(rownames(dataset) %in% genes$Group2)
    
  } else {
    
    colnames(genes) = "Group1"
    id_used_3 <- which(rownames(dataset) %in% genes$Group1)
    
  }
  
  # if no genes in the gene list are in the dataset 
  if (length(id_used_3) == 0) {
    
    patients_assignment_vector <- rep(0, ncol(dataset))
    print(paste0("Gene(s) not found: ", genes$Group1))
    
    return(patients_assignment_vector)
    
  } else { # if there are genes...
    
    if (colnames(dataset)[1] == "Hugo_Symbol") {
      dataset <- dataset[, 3:ncol(dataset)]
    } else {
      dataset <- dataset
    }
    
    dataset2 <- dataset[id_used_3, ] # select patient data for genes in the gene list
    
    CL3_z <- standardize(dataset2) # z-score
    
    patients_assignment_vector <- rep(0, ncol(CL3_z)) # patient number 
    q_vals <- rep(0, nrow(CL3_z)) # gene number 
    
    for (j in 1:nrow(CL3_z)) { # row by row, reurn the number at 20% or 50% etc quantile
      q_vals[j] <- quantile(as.numeric(CL3_z[j, 1:ncol(CL3_z)]), 
                            quantile_used, 
                            na.rm = T)
    }
    
    for (i in 1:ncol(CL3_z)) {
      
      q_comp = rep(0, nrow(CL3_z))
      
      for (j in 1:nrow(CL3_z)) {
        q_comp[j] <- CL3_z[j, i] > q_vals[j] # q_comp contains "TRUE" or "FALSE"
      }
      
      q_comp <- as.numeric(q_comp) # q_comp contains "1" or "0"
      
      # Percentage of genes in one column 
      if (length(which(q_comp == 1))/nrow(CL3_z) > percentage_patients) {
        patients_assignment_vector[i] = 1
      }
      else {
        patients_assignment_vector[i] = 0
      }
      
    }
    
    return(patients_assignment_vector)
    
  }
}
