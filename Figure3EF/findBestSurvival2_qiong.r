findBestSurvival2 <- function(dataset, # gene expression data
                              genelist, # gene list (one by one)
                              quantile_list, # seq(0.3, 0.7, by = 0.05)
                              cor_data, # survival data
                              print_plot = FALSE, 
                              effect_type = TRUE){
  # change gene to data.frame
  if (class(genelist) != 'data.frame'){
    genelist <- data.frame(genelist)
  }
  
  # change colname to "Group1"
  if (colnames(genelist) != 'Group1'){
    colnames(genelist) = 'Group1'
  }
  
  if (missing(effect_type))
  {
    effect_type <- FALSE
  }
  
  if (missing(print_plot))
  {
    print_plot <- FALSE
  } else if (effect_type == TRUE){ # if there are effects, no plotting
    print_plot <- FALSE
  }
  
  p_comp = Inf # start with the highest...
  lowest_quantile = NULL # Start with the lowest quantile
  q_vector <- quantile_list
  cor_data2 <- cor_data
  
  if (length(q_vector) > 1){
    # try all the elements in the q_vector and find the smallest p value (p_comp), 
    # the coresponding q_vector (lowest_quantile) 
    # and the "survival data + (0/1) with this element in q_vector" (cor_data_chosen)
    
    for (i in 1:length(q_vector)){
      # Package: siggitRausti (return "0" or "1" for each patient)
      # source('prepareSurvivalDataTCGA_qiong.R')
      patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,
                                                            genelist,
                                                            q_vector[i],
                                                            0.6) # percentage_patients This 0.75 value needs optimization!!
      
      if (all(patients_assignment_vector == 0)){ # all "0"
        
        lowest_quantile = q_vector[i]
        break
        
      } else {
        cor_data4 <- cbind(cor_data2, as.factor(patients_assignment_vector))
        colnames(cor_data4)[ncol(cor_data4)] <- c('Cl3_high_low')
       
        # Package "survival"
        res.cox2 <- coxph(Surv(survival, status) ~ Cl3_high_low, data = cor_data4)
        sum_coxph <- summary(res.cox2)
        pval_to_test <- sum_coxph$sctest[3]
        print(pval_to_test)
        
        # find the minimal p value, the coresponding vector in the quantile 
        # and new "suvival data + 0/1"
        if (pval_to_test < p_comp){
          
          p_comp <- pval_to_test
          lowest_quantile <- q_vector[i]
          cor_data_chosen <- cor_data4
          
        } else if (pval_to_test >= p_comp) {
          
          p_comp <- p_comp
          lowest_quantile <- lowest_quantile
          cor_data_chosen <- cor_data_chosen
          # lowest_quantile <- lowest_quantile
          
        }
        
      }
      
    }
    
  } else {
    lowest_quantile = q_vector
  }
  
  # Now output the correct graph:
  # Package "survival", generate 0/1 as before, but with lowest_quantile
  patients_assignment_vector <- prepareSurvivalDataTCGA(dataset,
                                                        genelist,
                                                        lowest_quantile,
                                                        0.6)
  cor_data_chosen <- cbind(cor_data2,
                           as.factor(patients_assignment_vector))
  colnames(cor_data_chosen)[ncol(cor_data_chosen)] <- c('Cl3_high_low')
  
  if (all(patients_assignment_vector == 0)){
    
    type_of_effect <- 0
    
    if (print_plot == TRUE){
      
      p3 <- 'No plot will be plotted, since there is no significant segregation'
      print(p3)
      
    } else if (print_plot == FALSE){
      
      if (effect_type == FALSE){
        
        p3 <- 1
        
      } else {
        
        if (type_of_effect == 0){
          p3 <- 'No effect!'
        }
        
      }
      
    }
    
  } else {
    # package "survival"
    res.cox_cand <- coxph(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
    
    # return -1, 0 or +1
    type_of_effect <- sign(res.cox_cand$coefficients)
    
    # package "survival"
    fit_chosen <- survfit(Surv(survival, status) ~ Cl3_high_low, data = cor_data_chosen)
    
    if(print_plot == TRUE){
      # plotting survival curve
      # package: "survminer"
      p3 <- ggsurvplot(fit_chosen, 
                       data = cor_data_chosen, 
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
      
      print(paste0("Q-value chosen is ", lowest_quantile))
      
    } else if (print_plot == FALSE){
      
      if (effect_type == FALSE){
        
        hehe <- summary(res.cox_cand)
        p3 <- hehe$sctest[3] # This is the p-value - return that!
        
      } else {
        hehe <- summary(res.cox_cand)
        p_value2 <- hehe$sctest[3]
        
        # res.cox_cand$coefficients (type_of_effect)
        if (type_of_effect == 1 & p_value2 < 0.05){
          p3 <- 'Beneficial!'
        } else if (type_of_effect == -1 & p_value2 < 0.05){
          p3 <- 'Harmful!'
        } else {
          p3 <- 'No effect!'
        }
      }
      
    }
  }
  
  return(p3)
  
}
