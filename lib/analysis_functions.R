generate_var_set_combinations <- function(...) {
  lists <- list(...)
  names_lists <- lapply(lists, names)
  combination_df <- expand.grid(names_lists)
  colnames(combination_df) <- paste0("list", seq_along(lists))
  return(combination_df)
}

produce_roc<-function(scores_table,tool){
  proc<-pROC::roc(response=scores_table$clinvar_class,predictor=scores_table%>%pull(tool),ci=TRUE,of='auc')
  return(proc)
}

calculate_roc_metrics<-function(proc_data,selected_main_var_sets,selected_sub_var_sets,selected_missing_var_sets){
  all_var_set_combinations<-generate_var_set_combinations(selected_main_var_sets,selected_sub_var_sets,selected_missing_var_sets)
  colnames(all_var_set_combinations)<-c('main_var_set','sub_var_set','missing_var_set')
  roc_metrics_table<-NULL
  rocs_table<-NULL
  procs<-list()
  #balancing_options<-c('all_class','balanced_class')
  missing_table<-NULL
  for (i in 1:nrow(all_var_set_combinations)){
    main_var_set<-as.character(all_var_set_combinations%>%slice(i)%>%pull(main_var_set))
    sub_var_set<-as.character(all_var_set_combinations%>%slice(i)%>%pull(sub_var_set))
    missing_option<-as.character(all_var_set_combinations%>%slice(i)%>%pull(missing_var_set))
    if (!(main_var_set%in%names(procs))){procs[[main_var_set]]<-list()}
    if (!(sub_var_set%in%names(procs[[main_var_set]]))){procs[[main_var_set]][[sub_var_set]]<-list()}
    if (!(missing_option%in%names(procs[[main_var_set]][[sub_var_set]]))){procs[[main_var_set]][[sub_var_set]][[missing_option]]<-list()}
    message(glue('{i}/{nrow(all_var_set_combinations)}:analyzing {main_var_set}:{sub_var_set}:{missing_option}'))
    # calculate rocs
    var_set<-proc_data%>%filter(main_var_sets[[main_var_set]]&sub_var_sets[[sub_var_set]]&missing_var_sets[[missing_option]])
    for (tool in tools_list) {
      #print(tool)
      tool_var_set<-var_set%>%filter(!is.na(!!sym(tool)))
      if(nrow(tool_var_set%>%filter(!is.na(!!sym(tool)))%>%count(clinvar_class))!=2){
        message(glue('There werent enough cases/controls for {tool}'))
        next()
      }
      procs[[main_var_set]][[sub_var_set]][[missing_option]][[tool]]<-list()
      proc <- suppressMessages(produce_roc(tool_var_set,tool))
      
      procs[[main_var_set]][[sub_var_set]][[missing_option]][[tool]]<-proc
      # calculate confusion matrix
      tool_cat<-stringr::str_replace(tool,'_score','_pred')
      tool_metrics<-data.frame(main_var_set,
                               sub_var_set,
                               missing_option,
                               total=nrow(tool_var_set),
                               tool=tool,
                               data.frame(as.numeric(proc$ci))%>%
                                 mutate(names=c('AUC2.5','AUC50.','AUC97.5'))%>%
                                 pivot_wider(names_from = names,values_from = as.numeric.proc.ci.))
      if (tool_cat%in%colnames(var_set)){
        
        var_set <- var_set %>% mutate(
          revel_pred.dbnsfp4.5a=ifelse(revel_score.dbnsfp4.5a>0.5,'P/LP','B/LB'),
          alphamissense_pred.dbnsfp4.5a=ifelse(alphamissense_pred.dbnsfp4.5a=='A',NA,alphamissense_pred.dbnsfp4.5a),
          preds = !!sym(tool_cat),
          preds = case_when(
            (preds=='P' & tool_cat=='mutationtaster_pred.dbnsfp4.5a')~'B',
            TRUE~preds
          ),
          preds = case_when(
            preds %in% c('A','P', 'D','H','M','P/LP') ~ 'P/LP',
            preds %in% c('B', 'T','N','L','B/LB') ~ 'B/LB',
            TRUE ~ NA
          ))%>%
          
          mutate(preds = factor(preds),
                 clinvar_class = factor(clinvar_class))
        if (var_set%>%count(preds)%>%nrow()<2){next}
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
      roc_metrics_table<-roc_metrics_table%>%bind_rows(tool_metrics)
    }
    to_ret<-list()
    to_ret[['procs']]<-procs
    to_ret[['roc_metrics_table']]<-roc_metrics_table
  }
  return(to_ret)
}
