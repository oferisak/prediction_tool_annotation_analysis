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

calculate_roc_metrics<-function(proc_data,var_sets_list,tools_to_test,complete=TRUE){
  roc_metrics_table<-NULL
  rocs_table<-NULL
  procs<-list()
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
    for (tool in tools_list) {
      #print(tool)
      tool_var_set<-var_set%>%filter(!is.na(!!sym(tool)))
      if(nrow(tool_var_set%>%filter(!is.na(!!sym(tool)))%>%count(clinvar_class))!=2){
        message(glue('There werent enough cases/controls for {tool}'))
        next()
      }
      proc <- suppressMessages(produce_roc(tool_var_set,tool))
      procs[[paste0(c(row%>%as.character(),tool,ifelse(complete,'complete','full')),collapse='|')]]<-proc
      # calculate confusion matrix
      tool_cat<-stringr::str_replace(tool,'_score','_pred')
      proc_coords<-proc%>%
        coords(ret='all',transpose=F)
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
    to_ret<-list()
    to_ret[['procs']]<-procs
    to_ret[['roc_metrics_table']]<-roc_metrics_table
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
