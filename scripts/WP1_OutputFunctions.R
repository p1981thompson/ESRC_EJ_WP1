
if (!require("here")) install.packages("here")

# lpa_enum_table: summarise output from class enumeration process ----
### function takes single argument of type mplus.model.list (created by readModels function)
### must also have MplusAutomation package loaded, and scientific notation turned off

lpa_enum_table <- function(output = NA){

  # Extract initial summary table from model output
  mm_summaries <- mixtureSummaryTable(output, keepCols = c("Title", "Classes", "Parameters", "Observations", 
   "LL", "AIC", "BIC", "aBIC", "T11_VLMR_PValue", "T11_LMR_PValue")) %>%

    # Exclude duplicates (e.g., if re-run to extract individual level data)
    distinct(Title, .keep_all = TRUE) %>%

    # Order by number of classes
    arrange(Classes) %>%

    # Create information indices
    mutate(CAIC = -2*LL + Parameters*(log(Observations) + 1),
           AWE = -2*LL + Parameters*(log(Observations) + 1.5),
           SIC = -0.5*BIC,
           BF = round(exp(SIC-lag(SIC)), 3)) 

  # Use SIC to compute approximate correct model probability out of k-class models
    max_sic <- max(mm_summaries$SIC, na.rm = TRUE)

    # Step 1 of cmP computation (top half of equation)
    mm_summaries <- mm_summaries %>%
      mutate(cmp_top = exp(SIC - max_sic))

    cmp_bottom <- sum(mm_summaries$cmp_top, na.rm = TRUE)

  # Add cmP_k to table output, and format for output
  mm_summaries <- mm_summaries %>%
    mutate(cmP_k = round(cmp_top/cmp_bottom, 2)) %>% 
    select(Title, Classes, LL, Parameters, BIC, CAIC, AWE, T11_VLMR_PValue,
           T11_LMR_PValue, BF, cmP_k) %>%
    rename(Specification = Title,
           VLMR_p = T11_VLMR_PValue,
           LMR_p = T11_LMR_PValue) %>%
    mutate(across(where(is.numeric), round, 2)) %>%
    mutate(VLMR_p = pvalue(VLMR_p, accuracy = 0.01),
           LMR_p = pvalue(LMR_p, accuracy = 0.01),
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
                                                           "LL", "AIC", "BIC", "aBIC", "T11_VLMR_PValue", "T11_LMR_PValue")) %>% 
    
    # Exclude duplicates (e.g., if re-run to extract individual level data)
    distinct(Title, .keep_all = TRUE) %>% 
    
    # Order by number of classes
    arrange(Classes) %>% 
    
    # Create information indices
    mutate(SIC = -0.5*BIC,
           BF = round(exp(SIC-lag(SIC)), 3))               
  
  # Use SIC to compute approximate correct model probability out of k-class models
  max_sic <- max(mm_summaries$SIC, na.rm = TRUE) 
  
  # Step 1 of cmP computation (top half of equation)
  mm_summaries <- mm_summaries %>% 
    mutate(cmp_top = exp(SIC - max_sic))
  
  cmp_bottom <- sum(mm_summaries$cmp_top, na.rm = TRUE)
  
  # Add cmP_k to table output, and format for output
  mm_summaries <- mm_summaries %>% 
    mutate(cmP_k = round(cmp_top/cmp_bottom, 2)) %>% 
    select(Title, Classes, LL, Parameters, AIC, BIC, aBIC, T11_VLMR_PValue, 
           T11_LMR_PValue) %>% 
    rename(Specification = Title, 
           VLMR_p = T11_VLMR_PValue,
           LMR_p = T11_LMR_PValue) %>% 
    mutate(across(where(is.numeric), round, 2)) %>% 
    mutate(VLMR_p = pvalue(VLMR_p, accuracy = 0.01), 
           LMR_p = pvalue(LMR_p, accuracy = 0.01))
  
  mm_summaries
}

# add_bLRT: add bLRT from model re-runs to summary table ----
### takes argument of output from the rerun (including bLRTs) plus the original summary table
add_bLRT <- function(rerun_output = NA, orig_summary = NA){
  
  if ("bLRT_p" %in% names(orig_summary)) {
    orig_summary <- select(orig_summary, -bLRT_p)
  }
  
  mm_summaries <- mixtureSummaryTable(rerun_output, keepCols = c("Title", "Classes", "Parameters",
                                                                 "LL", "BLRT_PValue")) %>%
    rename(Specification = Title, 
           bLRT_p = BLRT_PValue) %>% 
    mutate(across(where(is.numeric), round, 2)) %>% 
    mutate(bLRT_p = pvalue(bLRT_p, accuracy = 0.01)) %>% 
    right_join(orig_summary, by = c("Specification", "Classes", "Parameters", "LL"), copy = TRUE) %>% 
    arrange(Classes) %>% 
    relocate(bLRT_p, .after = LMR_p)
  
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

# getSection: extract file info ----
###Borrowed code from here: https://r-forge.r-project.org/scm/viewvc.php/pkg/R/utilityFunctions.R?view=markup&root=mplusautomation
####Needs some refinements.
getSection <- function(sectionHeader, outfiletext, headers="standard", omit=NULL) {
  #encode the top-level major headers here, but allow for custom headers to be passed in
  #omit allows for one or more strings from headers not to be considered
  #just used for factor score statistics at the moment (these include a SAMPLE STATISTICS section)
  if (headers[1L] == "standard") headers <- c("INPUT INSTRUCTIONS", "SUMMARY OF ANALYSIS",
                                              "SUMMARY OF DATA FOR THE FIRST DATA SET", "SUMMARY OF DATA FOR THE FIRST REPLICATION",
                                              "SUMMARY OF MISSING DATA PATTERNS FOR THE FIRST REPLICATION",
                                              "SUMMARY OF MISSING DATA PATTERNS",
                                              "COVARIANCE COVERAGE OF DATA FOR THE FIRST REPLICATION",
                                              "SAMPLE STATISTICS", "SAMPLE STATISTICS FOR THE FIRST REPLICATION",
                                              "CROSSTABS FOR CATEGORICAL VARIABLES", "UNIVARIATE PROPORTIONS AND COUNTS FOR CATEGORICAL VARIABLES",
                                              "RANDOM STARTS RESULTS RANKED FROM THE BEST TO THE WORST LOGLIKELIHOOD VALUES",
                                              "TESTS OF MODEL FIT", "MODEL FIT INFORMATION", "CLASSIFICATION QUALITY",
                                              "FINAL CLASS COUNTS AND PROPORTIONS FOR THE LATENT CLASSES",
                                              "FINAL CLASS COUNTS AND PROPORTIONS FOR THE LATENT CLASS PATTERNS",
                                              "LATENT TRANSITION PROBABILITIES BASED ON THE ESTIMATED MODEL",
                                              "FINAL CLASS COUNTS AND PROPORTIONS FOR EACH LATENT CLASS VARIABLE",
                                              "CLASSIFICATION OF INDIVIDUALS BASED ON THEIR MOST LIKELY LATENT CLASS MEMBERSHIP",
                                              "Average Latent Class Probabilities for Most Likely Latent Class Membership \\(Row\\)",
                                              "MODEL RESULTS", "LOGISTIC REGRESSION ODDS RATIO RESULTS", "RESULTS IN PROBABILITY SCALE",
                                              "IRT PARAMETERIZATION IN TWO-PARAMETER LOGISTIC METRIC",
                                              "LATENT CLASS ODDS RATIO RESULTS", "LOGRANK OUTPUT", "STANDARDIZED MODEL RESULTS", 
                                              "R-SQUARE", "QUALITY OF NUMERICAL RESULTS", "TECHNICAL OUTPUT", "TECHNICAL \\d+ OUTPUT",
                                              "TECHNICAL 5/6 OUTPUT", 
                                              "TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS",
                                              "STANDARDIZED TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS", "CONFIDENCE INTERVALS OF MODEL RESULTS",
                                              "CONFIDENCE INTERVALS FOR THE LOGISTIC REGRESSION ODDS RATIO RESULTS",
                                              "CREDIBILITY INTERVALS OF MODEL RESULTS",
                                              "CONFIDENCE INTERVALS OF STANDARDIZED MODEL RESULTS",
                                              "CREDIBILITY INTERVALS OF STANDARDIZED MODEL RESULTS",
                                              "CONFIDENCE INTERVALS OF TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT, AND DIRECT EFFECTS",
                                              "CONFIDENCE INTERVALS OF STANDARDIZED TOTAL, TOTAL INDIRECT, SPECIFIC INDIRECT,", #omitted "AND DIRECT EFFECTS"
                                              "EQUALITY TESTS OF MEANS ACROSS CLASSES USING POSTERIOR PROBABILITY-BASED",
                                              "THE FOLLOWING DATA SET\\(S\\) DID NOT RESULT IN A COMPLETED REPLICATION:",
                                              "RESIDUAL OUTPUT", "MODEL MODIFICATION INDICES", "MODEL COMMAND WITH FINAL ESTIMATES USED AS STARTING VALUES",
                                              "FACTOR SCORE INFORMATION \\(COMPLETE DATA\\)", "SUMMARY OF FACTOR SCORES", "PLOT INFORMATION", "SAVEDATA INFORMATION",
                                              "SAMPLE STATISTICS FOR ESTIMATED FACTOR SCORES", "DIAGRAM INFORMATION",
                                              "Beginning Time:\\s*\\d+:\\d+:\\d+", "MUTHEN & MUTHEN"
  )
  if (!is.null(omit)) headers <- headers[which(!headers %in% omit)] #drop omit
  beginSection <- grep(sectionHeader, outfiletext, perl=TRUE)[1]
  #if section header cannot be found, then bail out
  if (is.na(beginSection)) return(NULL)
  #form alternation pattern for regular expression (currently adds leading and trailing spaces permission to each header)
  headerRegexpr <- paste("(", paste(gsub("(.*)", "^\\\\s*\\1\\\\s*$", headers, perl=TRUE), sep="", collapse="|"), ")", sep="") 
  headerLines <- grep(headerRegexpr, outfiletext, perl=TRUE)
  subsequentHeaders <- which(headerLines > beginSection)
  if (length(subsequentHeaders) == 0) nextHeader <- length(outfiletext) #just return the whole enchilada
  else nextHeader <- headerLines[subsequentHeaders[1]] - 1
  section.found <- outfiletext[(beginSection+1):nextHeader]
  attr(section.found, "lines") <- beginSection:nextHeader
  return(section.found)
}

# get_optseed: extract optseed from initial models ----
get_optseed <- function(filename){
  ############################################################################
  outfiletext_pt<-readLines(filename)
  test<-getSection("RANDOM STARTS RESULTS RANKED FROM THE BEST TO THE WORST LOGLIKELIHOOD VALUES", outfiletext=outfiletext_pt, headers="standard", omit=NULL)
  ############################################################################
  test<-test[str_detect(test,"[[:alpha:]]")==FALSE]
  for(i in 1:length(test))
  {
    test[i]<-ifelse(test[i]=="",NA,test[i])
  }
  test<-na.omit(test)
  test<-str_split_fixed(test, "\\s+", 4)[,-1]
  test<-data.frame(test)
  test<-sapply(test,as.numeric)
  test<-data.frame(test)
  names(test)<-c("LoMax","seed","starts")
  #optiseed<-test$seed[which.min(test$LoMax)]   # incorrect?
  optiseed<- test$seed[1]  # edited to best value
  optiseed
}

# mm_extract_data: refit subset of models to save out participant classifications ----
### this currently requires the original model list to still be in the environment
### NOTE: increase lrtstarts if best values not replicated https://www.statmodel.com/examples/webnotes/webnote14.pdf
mm_extract_data <- function(orig_mods = NA,        # list of original models (in environment)
                            orig_output = NA, 
                            candidate_mods = NA,   # candidate number of classes
                            filepath = NA,         # folder to save models to
                            analysis_id = "sv",
                            rerun = TRUE,          # whether to re-run models to extract data, set to false if just loading 
                            optseed = TRUE,        # vector of corresponding optseed values (manually identified from output)
                            one_fit = TRUE) {      # whether a one-class model was fitted (if not, adjusts selection from list)
  
  # Make sure output ordered by class number
  orig_output <- list.sort(orig_output, summaries$NLatentClasses)
  
  # Adjust index of model if one-class model not fitted 
  if (one_fit == FALSE){
    candidate_mods = candidate_mods-1
  }
  
  # Extract file name
  model_name <- str_remove(deparse(substitute(orig_mods)), "m_")
  
  # Update scripts 
  if (rerun == TRUE){
    
    
    for (model in candidate_mods){
      
      
      # Adjust name to reflect class n (if necessary)
      n_classes <- ifelse(one_fit == TRUE, model, model+1) 
      
      # If using optseed
      if (optseed == TRUE) {
        
        # Extract filename from original output, extract optseed
        mod_file <- names(orig_output)[[model]]
        optseed_path <- paste0(here::here(), "/scripts/", filepath, "/", mod_file)
        optseed_val <- get_optseed(optseed_path)
        
        body <- update(orig_mods[[model]],
                       ANALYSIS = as.formula(sprintf("~ 'estimator = mlr; type = mixture; 
                                                   starts = 0;
                                                   optseed = %d;
                                                   processors = 4(starts);
                                                   lrtstarts = 0 0 100 20;'", optseed_val)),
                      # VARIABLE = ~ . + "AUXILIARY = cidB3153;",
                       OUTPUT = ~ "TECH7 TECH11 TECH14 stdyx CINTERVAL ENTROPY;",
                       SAVEDATA = as.formula(sprintf(" ~ 'FILE IS sv_%s_%dclass.dat; SAVE = cprobabilities;'", model_name, n_classes)))
      } else {
        body <- update(orig_mods[[model]],
                       OUTPUT = ~ "TECH7 TECH11 TECH14 stdyx CINTERVAL ENTROPY;",
                       ANALYSIS = ~ "estimator = mlr; type = mixture; starts = 1000 250; 
              processors = 4(STARTS);",
                       SAVEDATA = as.formula(sprintf(" ~ 'FILE IS sv_%s_%dclass.dat; SAVE = cprobabilities;'", model_name, n_classes)))
      }
      
      
      mplusModeler(body, sprintf("%s/%s_%s_rerun_%dclass.dat", filepath, analysis_id, model_name, n_classes), run = TRUE)
    }
  }
  filefilter = paste0(model_name, "_rerun")
  fmm_classification <- readModels(target = filepath, filefilter = filefilter)
  fmm_classification
  
}




# plotMixtures_simpView: plot standardised class means for task variables ----
### currently makes use of standard plotMixtures function, but could be extract manually for more flexibility (link saved in slack)
plotMixtures_simpView <- function(output = NA){

  # Title
  plot_title <- deparse(substitute(output))
  
  # Plot
  fig <- plotMixtures(output, variables = c("Combacc","Naraacc", "Naracomp", "Woldcomp"),
               coefficients = "stdyx.standardized") +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1)) +
    scale_x_discrete(limits = c("Combacc","Naraacc", "Naracomp", "Woldcomp")) + 
    ggtitle(plot_title)
  print(fig)  
  
  # Save out for clearer inspection
  plot_title <- substr(plot_title, 5, 6)
  filename <- paste0("../output/figures/SimpView_", plot_title, "candidateProfs.tiff")
  ggsave(filename, plot = fig, dpi = 600, height = 9, width = 10)
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
    group_by(.[[1]], .[[2]], .[[3]], .[[4]]) %>%
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
  # if four
  if (length(output) == 4){
    plot_base <- all_classes %>%
      mutate_at(1:length(output), as.factor) %>% 
      group_by(.[[1]], .[[2]], .[[3]], .[[4]]) %>%
      summarise(n = n(), mean_comp = mean(NARACOMP)) %>% 
      ggplot(aes(y = n, axis1 = .[[1]], axis2 = .[[2]], axis3 = .[[3]], axis4 = .[[4]], fill = mean_comp))  
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



