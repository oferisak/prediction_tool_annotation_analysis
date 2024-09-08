generate_var_set_combinations <- function(...) {
  lists <- list(...)
  names_lists <- lapply(lists, names)
  combination_df <- expand.grid(names_lists)
  colnames(combination_df) <- paste0("filter", seq_along(lists))
  combination_df<-combination_df%>%mutate(across(all_of(colnames(combination_df)),as.character))
  return(combination_df)
}

produce_roc<-function(scores_table,tool){
  proc<-pROC::roc(response=scores_table$clinvar_class,predictor=scores_table%>%pull(tool),ci=TRUE,of='auc')
  return(proc)
}

combine_filter_sets<-function(var_sets, var_set_names) {
  if (!all(var_set_names %in% names(var_sets))) {
    stop("One or more specified elements do not exist in the list")
  }
  # Extract the vectors corresponding to the specified elements
  selected_vectors <- var_sets[var_set_names]
  # Compute the intersection where all are TRUE
  result <- Reduce(`&`, selected_vectors)
  return(result)
}

adaptive_density <- function(score, density_object) {
  eval.points <- score
  density_value <- predict(density_object, x = eval.points)
  return(density_value)
}


# claude solution for posterior prob calculation
calculate_plp_probability <- function(scores, labels, dataset_ratio, pop_ratio,bw_method) {
  # Ensure labels are factor
  labels <- factor(labels, levels = c("P/LP", "B/LB"))
  # Calculate likelihoods
  #density_plp <- density(scores[labels == "P/LP"],bw = bw_method)
  density_plp <- tryCatch({
    density(scores[labels == "P/LP"], bw = bw_method)
  }, error = function(e) {
    message("Error with initial bw_method. Trying a different bw_method.")
    # You can specify an alternative bw_method here
    alternative_bw_method <- "ucv"  # Example alternative method
    density(scores[labels == "P/LP"], bw = alternative_bw_method)
  })
  #max_plp_dens<-max(density_plp$y)
  #density_blb <- density(scores[labels == "B/LB"],bw = bw_method)
  density_blb <- tryCatch({
    density(scores[labels == "B/LB"], bw = bw_method)
  }, error = function(e) {
    message("Error with initial bw_method. Trying a different bw_method.")
    # You can specify an alternative bw_method here
    alternative_bw_method <- "ucv"  # Example alternative method
    density(scores[labels == "B/LB"], bw = alternative_bw_method)
  })
  # if the score is overly dispersed, return NA
  if (sd(density_plp$y)>50){
    message(glue('Score is overly dispersed, will skip posterior probability calculation'))
    return(c())}
  #max_blb_dens<-density_blb$y[which.max(density_plp$y)]
  #print(glue('max plp {max_plp_dens}, plp sd is {sd(density_plp$y)}, IQR is {IQR(density_plp$y)}'))
  
  # Function to get density at a specific point
  get_density <- function(x, d) {
    idx <- which.min(abs(d$x - x))
    return(d$y[idx])
  }
  
  # Calculate posterior probabilities
  posterior_probs <- sapply(scores, function(score) {
    likelihood_plp <- get_density(score, density_plp)
    likelihood_blb <- get_density(score, density_blb)
    
    # likelihood_PLP <- approx(density_plp$x, density_plp$y, xout = score)$y
    # likelihood_BLB <- approx(density_blb$x, density_blb$y, xout = score)$y
    
    # likelihood_plp <- adaptive_density(score, density_plp)
    # likelihood_blb <- adaptive_density(score, density_blb)
    
    # Using dataset_ratio as prior
    # prior_plp <- dataset_ratio
    # prior_blb <- 1 - dataset_ratio
    
    # posterior_plp <- (likelihood_plp * prior_plp) / 
    #   (likelihood_plp * prior_plp + likelihood_blb * prior_blb)
    
    # Adjust for population ratio
    # adjusted_posterior_plp <- (posterior_plp * pop_ratio) / 
    #   (posterior_plp * pop_ratio + (1 - posterior_plp) * (1 - pop_ratio))
    
    adjusted_posterior_plp <- (likelihood_plp * pop_ratio) / 
         (likelihood_plp * pop_ratio + likelihood_blb * (1-pop_ratio))
    
    return(adjusted_posterior_plp)
  })
  
  return(posterior_probs)
}


# the dataset ratio should not be based on the tool's ratio but rather the P/B ratio in the complete var set regardless of each score's missingness
# converted tool - whether the tools scores are converted (lower scores are more likely to be pathogenic)
calculate_dataset_posterior_prob<-function(tool_var_set,
                                           tool,
                                           dataset_ratio,
                                           original_dataset_ratio,
                                           pop_ratio=0.0441,
                                           converted_tool=FALSE,
                                           bw_method){
  data <- tool_var_set
  
  # Calculate dataset-specific ratio
  #dataset_ratio <- mean(data$clinvar_class == "P/LP")
  
  # adjust pop_ratio according to the dataset ratio 
  adjusted_pop_ratio<- pop_ratio * (dataset_ratio/original_dataset_ratio)
  #print(adjusted_pop_ratio)
  
  probabilities <- calculate_plp_probability(
    data%>%pull(tool), 
    data%>%pull(clinvar_class), 
    dataset_ratio, 
    adjusted_pop_ratio,
    bw_method
  )
  # if score overly dispresed, plp will return NA
  if (length(probabilities)==0){return(data.frame(tool=tool))}
  criteria<-cut(probabilities,breaks=c(0.0999,0.2108,0.6073,0.9811,1),labels=c('P','M','S','VS'))
  posterior_prob_df<-data.frame(tool,
                                score=data%>%pull(tool),
                                post_prob = probabilities,
                                criteria,
                                clinvar_class=tool_var_set$clinvar_class)
  if(converted_tool){
    posterior_prob_df_by_criteria<-posterior_prob_df%>%
      group_by(criteria)%>%
      slice_max(score,with_ties = F)
  }else{
    posterior_prob_df_by_criteria<-posterior_prob_df%>%
      group_by(criteria)%>%
      slice_min(score,with_ties = F)
  }
  posterior_prob_df_by_criteria<-posterior_prob_df_by_criteria%>%
    select(-clinvar_class)%>%
    left_join(posterior_prob_df%>%count(criteria,clinvar_class)%>%
                mutate(clinvar_class=make.names(clinvar_class))%>%
                group_by(criteria)%>%mutate(total=sum(n))%>%ungroup()%>%
                pivot_wider(names_from = clinvar_class,values_from = n))%>%
    ungroup()%>%
    filter(!is.na(criteria))
  if (nrow(posterior_prob_df_by_criteria)==0){
    message('Tool did not reach minimal posterior prob for supporting criteria, skipping..')
    return(data.frame(tool=tool))
  }
  return(posterior_prob_df_by_criteria)
}


# calculate_dataset_posterior_prob(tool_var_set,
#                                  tool,dataset_ratio = dataset_ratio,
#                                  original_dataset_ratio = original_dataset_ratio,
#                                  converted_tool=FALSE,bw_method = 'nrd0')

# bandwidth selection explanation: https://aakinshin.net/posts/kde-bw/ the recommended bw is SJ, but it returns sample is too sparse
calculate_roc_metrics<-function(proc_data,var_sets_list,
                                tools_to_test,
                                converted_scores,
                                complete=TRUE,save_rocs=TRUE,
                                calculate_posterior_probs=FALSE,
                                bw_method='nrd0'){
  roc_metrics_table<-NULL
  rocs_table<-NULL
  posterior_probs<-NULL
  procs<-list()
  original_dataset_ratio<-mean(proc_data$clinvar_class=='P/LP',na.rm = T)
  for (i in 1:nrow(var_sets_list)){
    row<-var_sets_list[i,,drop=F]
    message(glue('{i}/{nrow(var_sets_list)}:analyzing {paste0(row%>%slice(1),collapse=":")}'))
    # calculate rocs
    var_set<-proc_data%>%filter(combine_filter_sets(var_sets,row%>%as.character()))
    if (complete){
      original_nrow<-nrow(var_set)
      var_set<-var_set%>%filter(!if_any(all_of(tools_to_test), is.na))
      complete_nrow<-nrow(var_set)
      message(glue('Out of {original_nrow} variants, {complete_nrow} ({round(complete_nrow/original_nrow,3)}) had values in all of the tools to test'))
    }
    var_set<-var_set%>%mutate(revel_pred.dbnsfp4.5a=ifelse(revel_score.dbnsfp4.5a>0.5,'D','T'))
    var_set_n_plp<-sum(var_set$clinvar_class=='P/LP',na.rm = T) # its more correct to use the tool_var_set plp because its not fair to calculate sensitivity if the tool cant see tthe variant
    dataset_ratio<-mean(var_set$clinvar_class=='P/LP',na.rm = T)
    message(glue('With {var_set_n_plp} pathogenic variants, the P/B ratio in the {paste0(row%>%slice(1),collapse=":")} dataset is {dataset_ratio}'))
    cat('\n')
    for (tool in tools_to_test) {
      cat(sprintf("\rAnalyzing %s",tool))
      tool_var_set<-var_set%>%filter(!is.na(!!sym(tool)))
      if(nrow(tool_var_set%>%filter(!is.na(!!sym(tool)))%>%count(clinvar_class))!=2){
        message(glue('There werent enough cases/controls for {tool}'))
        next()
      }
      # calculate posterior probs
      if (calculate_posterior_probs){
        converted_tool<-ifelse(tool %in% converted_scores,TRUE,FALSE)
        tool_posterior_probs<-suppressMessages(calculate_dataset_posterior_prob(tool_var_set,
                                                               tool,
                                                               dataset_ratio,
                                                               original_dataset_ratio,
                                                               pop_ratio=0.0441,
                                                               converted_tool,
                                                               bw_method = bw_method))
                                                               
        if ('P.LP' %in% colnames(tool_posterior_probs)){
          tool_var_set_n_plp<-sum(tool_var_set$clinvar_class=='P/LP',na.rm = T) # its more correct to use the tool_var_set plp because its not fair to calculate sensitivity if the tool cant see tthe variant
          tool_posterior_probs<-tool_posterior_probs%>%
            mutate(sensitivity=P.LP/tool_var_set_n_plp)
        }
        posterior_probs<-posterior_probs%>%bind_rows(
          data.frame(
            row,
            complete=complete,
            total=nrow(tool_var_set),
            tool_posterior_probs
          )
        )
      }
      
      # calculate roc metrics
      proc <- suppressMessages(produce_roc(tool_var_set,tool))
      if (save_rocs){
        procs[[paste0(c(row%>%as.character(),tool,ifelse(complete,'complete','full')),collapse='|')]]<-proc
      }
      # calculate confusion matrix
      tool_cat<-stringr::str_replace(tool,'_score','_pred')
      proc_coords<-proc%>%
        coords(ret='all',transpose=F)
      recall_at_precision90<-proc_coords%>%
        select(threshold,precision,recall)%>%
        slice_min(order_by=abs(0.9-precision),with_ties = F)%>%
        select(threshold_ppv90=threshold,precision90=precision,recall_precision90=recall)
      
      percision_at_recall90<-proc_coords%>%
        select(threshold,precision,recall)%>%
        slice_min(order_by=abs(0.9-recall),with_ties = F)%>%
        select(threshold_recall90=threshold,recall90=recall,precision_recall90=precision)
      percision_at_recall99<-proc_coords%>%
        select(threshold,precision,recall)%>%
        slice_min(order_by=abs(0.99-recall),with_ties = F)%>%
        select(threshold_recall99=threshold,recall99=recall,precision_recall99=precision)
      tool_metrics<-data.frame(row,
                               complete=complete,
                               total=nrow(tool_var_set),
                               tool=tool,
                               data.frame(as.numeric(proc$ci))%>%
                                 mutate(names=c('AUC2.5','AUC50.','AUC97.5'))%>%
                                 pivot_wider(names_from = names,values_from = as.numeric.proc.ci.),
                               recall_at_precision90,
                               percision_at_recall90,
                               percision_at_recall99)
      if (tool_cat%in%colnames(var_set)){
        var_set <- var_set %>% mutate(
          #alphamissense_pred.dbnsfp4.5a=ifelse(alphamissense_pred.dbnsfp4.5a=='A',NA,alphamissense_pred.dbnsfp4.5a),
          preds = !!sym(tool_cat),
          preds = case_when(
            preds=='D'~ 'P/LP',
            preds=='T' ~ 'B/LB',
            .default= NA
          ))%>%
          
          mutate(preds = factor(preds),
                 clinvar_class = factor(clinvar_class))
        if (var_set%>%count(preds)%>%nrow()>=2){
          cm<-yardstick::conf_mat(data=var_set,truth='clinvar_class',estimate='preds')
          cm<-data.frame(TP=cm$table['P/LP','P/LP'],
                         TN=cm$table['B/LB','B/LB'],
                         FP=cm$table['P/LP','B/LB'],
                         FN=cm$table['B/LB','P/LP'])
          metrics_to_calc<-yardstick::metric_set(yardstick::sens,yardstick::spec,yardstick::ppv,yardstick::npv,yardstick::accuracy,yardstick::bal_accuracy)
          added_metrics<-metrics_to_calc(data=var_set,truth='clinvar_class',estimate='preds',event_level = 'second')%>%
            select(-.estimator)%>%pivot_wider(names_from = .metric,values_from = .estimate)
          tool_metrics<-tool_metrics%>%
            bind_cols(data.frame(cm,added_metrics))
        }
      }
      roc_metrics_table<-roc_metrics_table%>%bind_rows(tool_metrics)
    }
    cat("\n")
    to_ret<-list()
    to_ret[['procs']]<-procs
    to_ret[['roc_metrics_table']]<-roc_metrics_table
    to_ret[['posterior_prob']]<-posterior_probs
  }
  return(to_ret)
}

plot_roc<-function(raw_roc_table,overall_performance_table,top_tools,year_to_analyze){
  roc_table<-raw_roc_table%>%
    #filter(tool %in% tools_to_test)%>%
    filter(year==year_to_analyze)%>%
    left_join(overall_performance_table%>%
                rename(year=filter1))%>%
    mutate(tool=stringr::str_replace(tool,'_score.+','')%>%toupper(),
           tool_with_auc=ifelse(tool%in%top_tools$tool,glue('{tool} ({round(AUC50.,3)})'),'Other'),
           tool_with_auc=forcats::fct_reorder(tool_with_auc,desc(AUC50.)))
  colors <- setNames(colorRampPalette(c("darkred", "darkcyan", "darkorange",'darkmagenta'))(length(unique(roc_table$tool_with_auc))), setdiff(unique(roc_table$tool_with_auc),'Other'))
  colors["Other"] <- "lightgray" 
  print(colors)
  # Plot the ROC for all the tools and emphasize the top ones
  roc_plot<-
    roc_table%>%filter(tool_with_auc=='Other')%>%
    ggplot(aes(x=1-specificity,y=sensitivity,color=tool_with_auc,group=tool))+
    geom_line(linewidth=1)+
    geom_line(data=roc_table%>%filter(tool_with_auc!='Other'),linewidth=1)+
    geom_abline(slope = 1,intercept = 0,linetype=2,alpha=0.4)+
    scale_color_manual(values=colors)+
    theme_minimal()+
    theme(legend.position = 'top')+
    labs(color=NULL)
  print(roc_plot)
  return(roc_plot)
}
