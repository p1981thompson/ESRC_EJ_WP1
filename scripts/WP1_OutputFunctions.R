
# lpa_enum_table: summarise output from class enumeration process ----
### function takes single argument of type mplus.model.list (created by readModels function)
### must also have MplusAutomation package loaded, and scientific notation turned off
lpa_enum_table <- function(output = NA){
  
  # Extract initial summary table from model output
  mm_summaries <- mixtureSummaryTable(output, keepCols = c("Title", "Classes", "Parameters", "Observations", 
                                                           "LL", "AIC", "BIC", "aBIC", "T11_VLMR_PValue", "T11_LMR_PValue", "BLRT_PValue")) %>% 
    
    # Exclude duplicates (e.g., if re-run to extract individual level data)
    distinct(Title, .keep_all = TRUE) %>% 
    
    # Order by number of classes
    arrange(Classes) %>% 
    
    # Create information indices
    mutate(CAIC = -2*LL + Parameters*(log(Observations) + 1),
           AWE = -2*LL + Parameters*(log(Observations) + 1.5),
           SIC = -0.5*BIC,
           BF = round(exp(SIC-lead(SIC)), 3)) ########################################### PT PLEASE CHECK                
  
  # Use SIC to compute approximate correct model probability out of k-class models
    max_sic <- max(mm_summaries$SIC, na.rm = TRUE) 
    
    # Step 1 of cmP computation (top half of equation)
    mm_summaries <- mm_summaries %>% 
      mutate(cmp_top = exp(SIC - max_sic))
    
    cmp_bottom <- sum(mm_summaries$cmp_top, na.rm = TRUE)
  
  # Add cmP_k to table output, and format for output
  mm_summaries <- mm_summaries %>% 
    mutate(cmP_k = round(cmp_top/cmp_bottom, 2)) %>% ###################################### PT PLEASE CHECK  cmP
    select(Title, Classes, LL, Parameters, BIC, CAIC, AWE, T11_VLMR_PValue, 
           T11_LMR_PValue, BLRT_PValue, BF, cmP_k) %>% 
    rename(Specification = Title, 
           VLMR_p = T11_VLMR_PValue,
           LMR_p = T11_LMR_PValue,
           bLRT_p = BLRT_PValue) %>% 
    mutate(across(where(is.numeric), round, 2)) %>% 
    mutate(VLMR_p = pvalue(VLMR_p, accuracy = 0.01), 
           LMR_p = pvalue(LMR_p, accuracy = 0.01),
           bLRT_p = pvalue(bLRT_p, accuracy = 0.01),
           BF = ifelse(BF == 0.00, "<0.01",
                       ifelse(BF > 100, ">100", BF)))
  
  mm_summaries
}



# fmm_enum_table: summarise output from class enumeration process ----
### function takes single argument of type mplus.model.list (created by readModels function)
### must also have MplusAutomation package loaded, and scientific notation turned off
fmm_enum_table <- function(output = NA){
  
  # Extract initial summary table from model output
  mm_summaries <- mixtureSummaryTable(output, keepCols = c("Title", "Classes", "Parameters", "Observations", 
                                                           "LL", "AIC", "BIC", "aBIC", "T11_VLMR_PValue", "T11_LMR_PValue", "BLRT_PValue")) %>% 
    
    # Exclude duplicates (e.g., if re-run to extract individual level data)
    distinct(Title, .keep_all = TRUE) %>% 
    
    # Order by number of classes
    arrange(Classes) %>% 
    
    # Create information indices
    mutate(SIC = -0.5*BIC,
           BF = round(exp(SIC-lead(SIC)), 3)) ########################################### PT PLEASE CHECK                
  
  # Use SIC to compute approximate correct model probability out of k-class models
  max_sic <- max(mm_summaries$SIC, na.rm = TRUE) 
  
  # Step 1 of cmP computation (top half of equation)
  mm_summaries <- mm_summaries %>% 
    mutate(cmp_top = exp(SIC - max_sic))
  
  cmp_bottom <- sum(mm_summaries$cmp_top, na.rm = TRUE)
  
  # Add cmP_k to table output, and format for output
  mm_summaries <- mm_summaries %>% 
    mutate(cmP_k = round(cmp_top/cmp_bottom, 2)) %>% ###################################### PT PLEASE CHECK  cmP
    select(Title, Classes, LL, Parameters, AIC, BIC, aBIC, T11_VLMR_PValue, 
           T11_LMR_PValue, BLRT_PValue) %>% 
    rename(Specification = Title, 
           VLMR_p = T11_VLMR_PValue,
           LMR_p = T11_LMR_PValue,
           bLRT_p = BLRT_PValue) %>% 
    mutate(across(where(is.numeric), round, 2)) %>% 
    mutate(VLMR_p = pvalue(VLMR_p, accuracy = 0.01), 
           LMR_p = pvalue(LMR_p, accuracy = 0.01),
           bLRT_p = pvalue(bLRT_p, accuracy = 0.01))
  
  mm_summaries
}

# lpa_enum_elbow: plot output from class enumeration process ----
### function takes one argument for model output summary table, created by the above function
### also one optional argument for benchmark statistics
lpa_enum_elbow <- function(mm_summaries = NA, 
                          benchmark_stats = NA){
  elbow_fig <- mm_summaries %>% 
    select(Classes, LL, BIC, CAIC, AWE) %>%  
    pivot_longer(LL:AWE, names_to = "statistic", values_to = "value") %>%
    mutate(statistic_relevel = factor(statistic, levels = c("LL", "BIC", "CAIC", "AWE"))) %>% 
    ggplot(aes(x = as.factor(Classes), y = value)) +
    geom_point() +
    geom_line(group = 1) +
    xlab("n classes") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    facet_wrap(statistic_relevel ~., scales = "free") 

  if (!is.na(benchmark_stats)){
    elbow_fig <- elbow_fig +
      geom_hline(data = benchmark_stats, aes(yintercept = benchmark), linetype = "dashed", colour = "red")
  }
  
  elbow_fig
}

# fmm_enum_elbow: plot output from class enumeration process ----
### function takes one argument for model output summary table, created by the above function
### also one optional argument for benchmark statistics
fmm_enum_elbow <- function(mm_summaries = NA, 
                           benchmark_stats = NA){
  elbow_fig <- mm_summaries %>% 
    select(Classes, LL, AIC, BIC, aBIC) %>%  
    pivot_longer(LL:aBIC, names_to = "statistic", values_to = "value") %>%
    mutate(statistic_relevel = factor(statistic, levels = c("LL", "AIC", "BIC", "aBIC"))) %>% 
    ggplot(aes(x = as.factor(Classes), y = value)) +
    geom_point() +
    geom_line(group = 1) +
    xlab("n classes") +
    theme_bw() +
    theme(panel.grid.major.x = element_blank(),
          panel.grid.minor.y = element_blank()) +
    facet_wrap(statistic_relevel ~., scales = "free") 
  
  if (!is.na(benchmark_stats)){
    elbow_fig <- elbow_fig +
      geom_hline(data = benchmark_stats, aes(yintercept = benchmark), linetype = "dashed", colour = "red")
  }
  
  elbow_fig
}

# mm_extract_data: refit subset of models to save out participant classifications ----
### this currently requires the original model list to still be in the environment - could be improved
mm_extract_data <- function(orig_mods = NA,        # list of original models (in environment)
                            candidate_mods = NA,   # candidate number of classes
                            filepath = NA,         # folder to save models to
                            analysis_id = "sv",
                            rerun = TRUE,          # whether to re-run models to extract data, set to false if just loading 
                            one_fit = TRUE) {      # whether a one-class model was fitted (if not, adjusts selection from list)
  
  # Adjust index of model if one-class model not fitted 
  if (one_fit == FALSE){
    candidate_mods = candidate_mods-1
  }
  
  # Extract file name
  model_name <- str_remove(deparse(substitute(orig_mods)), "m_")
  
  # Update scripts 
  if (rerun == TRUE){
    for (model in candidate_mods){
      
      n_classes <- ifelse(one_fit == TRUE, model, model+1) # adjust name to reflect class n
      
      body <- update(orig_mods[[model]],
                     VARIABLE = ~ . + "IDVARIABLE IS yp_no;",
                     SAVEDATA = as.formula(sprintf(" ~ 'FILE IS sv_%s_%dclass.dat; SAVE = cprobabilities;'", model_name, n_classes)))
      
      mplusModeler(body, sprintf("%s/%s_%s_rerun_%dclass.dat", filepath, analysis, model_name, n_classes), run = TRUE)
    }
  }
  filefilter = paste0(model_name, "_rerun")
  fmm_classification <- readModels(target = filepath, filefilter = filefilter)
  fmm_classification
  
}


# plotMixtures_simpView: plot standardised class means for task variables ----
### currently makes use of standard plotMixtures function, but could be extract manually for more flexibility (link saved in slack)
plotMixtures_simpView <- function(output = NA){
  plotMixtures(output, variables = c("Nonwacc", "Wordacc","Naraacc", "Naracomp", "Woldcomp"),
               coefficients = "stdyx.standardized") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_x_discrete(limits = c("Nonwacc", "Wordacc","Naraacc", "Naracomp", "Woldcomp"))
}

# plotMixtures_compProf: plot standardised class means for task variables ----
### currently makes use of standard plotMixtures function, but could be extract manually for more flexibility
plotMixtures_compProf <- function(output = NA, select_model = NA){
  
  # Specify all models if no model selected
  if (is.na(select_model)){
    select_model <- 1:length(output)
  }
  
  # For each model...
  for (model in select_model){
    
    # Make separate plots of tasks for each factor
    
    acc <- plotMixtures(output[[model]], 
                        variables = c("Nonwacc", "Wordacc","Naraacc"),
                        coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Nonwacc", "Wordacc","Naraacc")) + 
      ggtitle("accuracy")
    
    comp <- plotMixtures(output[[model]], 
                         variables = c("Naracomp", "Woldcomp"),
                         coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Naracomp", "Woldcomp")) + 
      ggtitle("comprehension")
    
    rate <- plotMixtures(output[[model]], 
                         variables = c("Nararate"),
                         coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Nararate")) + 
      ggtitle("fluency")
    
    vocab <- plotMixtures(output[[model]], 
                          variables = c("Wiscvcb", "Woldvcb"),
                          coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Wiscvcb", "Woldvcb")) + 
      ggtitle("vocabulary")
    
    perfiq <- plotMixtures(output[[model]], 
                           variables = c("Wiscpcmp", "Wisccode", "Wiscparr", "Wiscbloc", "Wiscobja"),
                           coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Wiscpcmp", "Wisccode", "Wiscparr", "Wiscbloc", "Wiscobja")) + 
      ggtitle("performance IQ")
    
    execfun <- plotMixtures(output[[model]], 
                            variables = c("Wiscbwsp", "Cntspan", "Attnsel", "Attndiv", "Attnopp", "Sdqhyp"),
                            coefficients = "stdyx.standardized") +
      scale_shape_manual(values = c(18, 18, 18, 18)) + 
      scale_size_manual(values = c(5, 5, 5, 5)) +
      scale_linetype_manual(values = c("solid", "solid", "solid", "solid")) +
      theme(axis.text.x = element_text(angle = 45, vjust = 0.9, hjust = 1),
            axis.title.x = element_blank()) +
      scale_x_discrete(limits = c("Wiscbwsp", "Cntspan", "Attnsel", "Attndiv", "Attnopp", "Sdqhyp")) + 
      ggtitle("executive functions")
    
    # Create title and arrange plots together
    plot_title <- output[[model]][["input"]][["title"]] 
    fig <- ggarrange(acc, comp, rate, vocab, perfiq, execfun, rows = 3, common.legend = TRUE)
    fig <- annotate_figure(fig, top = text_grob(plot_title, color = "black", face = "bold", size = 14))
    print(fig)
    
    # Save out for clearer inspection
    plot_title <- str_replace_all(plot_title, "[^[:alnum:]]", "_")
    filename <- paste0("../output/figures/", plot_title, ".tiff")
    print(filename)
    ggsave(filename, plot = fig, dpi = 600, height = 9, width = 10)
  }
  
}

# classification_diagnostics ---- 
### not formatted in easy-to-read way, could maybe be improved
class_diag <- function(output = NA){
  {
    for (model in 1:length(output)) {
      
      # Print model title 
      print(paste0("CLASSIFICATION DIAGNOSTICS FOR ", toupper(output[[model]]$input$title)))
      
      # Entropy
      print(paste0("Entropy: ", output[[model]]$summaries$Entropy))
      cat("\n")
      
      # Average posterior class probability
      print("Average posterior class probability:") 
      ave_pp <- (diag(output[[model]][["class_counts"]][["classificationProbs.mostLikely"]]))
      print(ave_pp)
      
      ### Print additional warning if do not meet rule of thumb
      if(any(ave_pp < 0.7)){
        print("Class assignment accuracy potentially inadequate (Nagin, 2005)")
      }
      
      cat("\n")
      
      # Odds of correct classification ratio 
      est_prop <- output[[model]][["class_counts"]][["modelEstimated"]][["proportion"]]
      occ <- round(ave_pp/(1-ave_pp)/(est_prop/(1-est_prop)))
      
      print("Odds of correct classification ratio:")
      print(occ)
              
      ### Print additional interpretation
      if(all(occ > 5)){
        print("Good class separation and assignment accuracy (Nagin, 2005)")
      }
      
      cat("\n")
      
      # Modal class assignment proportion
      print("Model estimated proportion for each class:")
      print(est_prop)
      
      cat("\n")
      
      print("Proportion modally assigned to each class:")
      mcap <- output[[model]][["class_counts"]][["mostLikely"]][["proportion"]]
      print(mcap)
      
      cat("\n")
      print("--------------------------------------------")
      cat("\n")
    }}
  
}


# extract_classes: 
extract_classes <- function(output = NULL, type = NULL){
  
  # for each model in the output
  for (model in 1:length(output)){
    
    # extract class assignment for each participant
    class_assign <- as.data.frame(output[[model]]$savedata$C)
    
    # name output with number of classes in model
    n_classes <- max(unique(class_assign))
    var_name <- paste0(type, "_C_k", n_classes)
    assign(var_name, class_assign)
  }
  
  # combine all class assignments into single dataframe
  classes_name <- paste0(type, "_classes")
  list_c_vars <- ls(pattern = paste0(type, "_C_"))
  all_classes <- do.call(cbind, mget(list_c_vars))
  
  # extract comprehension ability for plotting too
  ability <- as.data.frame(output[[1]]$savedata[,1:5])
  ability_names <- names(output[[1]]$savedata[1:5])
  all_classes <- cbind(all_classes, ability)
  
  # rename and return to global environment
  names(all_classes) <- c(list_c_vars, ability_names)
  assign(classes_name, all_classes, env = .GlobalEnv)
  
  #PLOTTING
  #summarise transitions
  class_trans <- all_classes %>%
    mutate_at(1:length(output), as.factor) %>% 
    group_by(.[[1]], .[[2]], .[[3]]) %>%
    summarise(n = n(), mean_comp = mean(NARACOMP))
  
  # if two
  if (length(output) == 2){
    plot_base <- all_classes %>%
      mutate_at(1:length(output), as.factor) %>% 
      group_by(.[[1]], .[[2]]) %>%
      summarise(n = n(), mean_comp = mean(NARACOMP)) %>% 
      ggplot(aes(y = n, axis1 = .[[1]], axis2 = .[[2]], fill = mean_comp))
  }
  
  # if three
  if (length(output) == 3){
    plot_base <- all_classes %>%
      mutate_at(1:length(output), as.factor) %>% 
      group_by(.[[1]], .[[2]], .[[3]]) %>%
      summarise(n = n(), mean_comp = mean(NARACOMP)) %>% 
      ggplot(aes(y = n, axis1 = .[[1]], axis2 = .[[2]], axis3 = .[[3]], fill = mean_comp))
    
  }
  
  plot_title <- paste0("Group transitions for Model ", toupper(substr(type, nchar(type), nchar(type))), " Candidate Models")
  
  plot <- plot_base +
    geom_alluvium(aes(fill = mean_comp), aes.bind = TRUE, width = 1/12) +
    scale_fill_gradient(low = "blue", high = "red", name = "Comprehension \nscore") +
    geom_stratum(fill = "white", colour = "black", width = 1/4) +
    geom_text(stat = "stratum", aes(label = after_stat(stratum))) +
    scale_x_discrete(expand = c(.05, .05)) +
    labs(y = "Participants") +
    theme_minimal() +
    theme(legend.position = "right",
          axis.text = element_text(colour = "black")) +
    ggtitle(plot_title)
  
  print(plot)
}








################ GRAVEYARD ----


# ave_pp: compute the average posterior class probability (classification uncertainty) ----
### takes list of models as input
ave_pp <- function(output = NA){
  
  for (model in 1:length(output)) {
    
    
    print(paste0("Classification diagnostics - Average posterior class probability for:", output[[model]]$input$title))
    ave_pp <- (diag(output[[model]][["class_counts"]][["classificationProbs.mostLikely"]]))
    print(ave_pp)
    
    # Print additional warning if do not meet rule of thumb
    if(any(ave_pp < 0.7)){
      print("Class assignment accuracy potentially inadequate (Nagin, 2005)")
    }
    cat("\n")
  }}


# classification diagnostics
class_diagn2 <- function(output = NA){
  
  # Create a table from model output including title and entropy measure 
  diagnostics <- mixtureSummaryTable(output, keepCols = c("Title", "Classes", "Entropy"))
  
  # Extract average posterior class probability
  for (model in 1:length(output)) {

    Title <- output[[model]]$input$title
    ave_pp <- (diag(output[[model]][["class_counts"]][["classificationProbs.mostLikely"]]))

    
    diag_df <- data.frame(Title, avePP = list(ave_pp), prop = list(modest_prop)) %>% 
      set_names(c("Title", "class_id", "proportion"))
    
    diagnostics <- diagnostics %>% 
      left_join(diag_df, by = "c(Title") 
  } 
  
  # Join to model
    # diagnostics <- diagnostics %>% 
    #   unite(col = "class_avepp", contains("class_id"), na.rm = TRUE) %>% 
    #   group_by(Classes) %>% 
    #   mutate(class_no = row_number())# %>% 
      # ungroup() %>% 
      # pivot_wider(names_from = class_no, names_prefix = "avepp_c", values_from = class_avepp)
    
   # Extract class proportions

    for (model in 1:length(output)) {
      
      Title <- output[[model]]$input$title
      modest_prop <- output[[model]][["class_counts"]][["modelEstimated"]][["proportion"]]
      
      diag_df2 <- data.frame(Title, prop = list(modest_prop)) %>% 
        set_names(c("Title", "proportion"))
      
      diagnostics <- diagnostics %>% 
        left_join(diag_df2, by = c("Title")) 
    } 
    
    # Join to model
    diagnostics <- diagnostics %>% 
      unite(col = "class_prop", contains("proportion"), na.rm = TRUE) %>% 
      group_by(Classes) %>% 
      mutate(class_no = row_number()) # %>% 
    # ungroup() %>% 
    # pivot_wider(names_from = class_no, names_prefix = "avepp_c", values_from = class_avepp)
    
    
       
  print(diagnostics) 
  }
