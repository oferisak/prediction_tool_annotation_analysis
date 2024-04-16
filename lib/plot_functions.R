plot_hist_by_class<-function(proc_data,tools_to_plot,title=NULL){
  hist_plot<-
    proc_data%>%select(tools_to_plot,clinvar_class)%>%
    filter(!is.na(clinvar_class))%>%
    pivot_longer(-clinvar_class)%>%
    mutate(name=stringr::str_replace(name,'_score.+','')%>%toupper())%>%
    ggplot(aes(x=value,fill=clinvar_class))+
    geom_histogram(alpha=0.5,position='dodge')+
    facet_grid(.~name,scale='free')+
    labs(x=NULL,y='Number of variants',fill=NULL,title=title)+
    scale_fill_manual(values=c('darkcyan','darkred'))+
    theme_minimal()+
    theme(legend.position='top')
  box_plot<-
    proc_data%>%select(tools_to_plot,clinvar_class)%>%
    pivot_longer(-clinvar_class)%>%
    filter(!is.na(clinvar_class))%>% 
    mutate(name=stringr::str_replace(name,'_score.+','')%>%toupper())%>%
    ggplot(aes(y=value,fill=clinvar_class))+
    #ggdist::stat_histinterval(alpha=0.5)+
    geom_boxplot(position='dodge',aes(x=clinvar_class),alpha=0.4,outlier.shape=NA)+
    facet_wrap(name~.,scales='free')+
    scale_fill_manual(values=c('darkcyan','darkred'))+
    labs(fill=NULL,x=NULL,y='Score')+
    guides(fill='none')+
    coord_flip()+
    theme_minimal()+
    ylim(0,1)+
    theme(legend.position = 'top',strip.text = element_blank())
  
  to_ret<-hist_plot/box_plot+ plot_layout(heights = c(3, 1))
  return(to_ret)
}

# plot calibration - 
# to_plot - needs to be the filtered proc_data table with only the variant set you want to test on
# facet - a categorical column that should be present in the to_plot table. the calibration plot will have a color for each category
# balanced - if true, will balance P / B variants for each tool and facet combination
plot_calibration<-function(to_plot,tools_to_plot,facet,balanced=F,text_size=10){
  
  balanced_sample_size_per_tool<-
    to_plot %>%
    select(tools_to_plot, clinvar_class, facet) %>%
    pivot_longer(-c(clinvar_class,facet))%>%
    filter(!(is.na(value)|is.na(!!sym(facet)))) %>%
    group_by(name,!!sym(facet), clinvar_class)%>%
    count()%>%ungroup()%>%group_by(name,!!sym(facet))%>%summarize(balanced_n=min(n))
  print(balanced_sample_size_per_tool)
  to_plot<-
    to_plot%>%
    select(clinvar_class,tools_to_plot,facet)%>%
    filter(!is.na(clinvar_class))%>%
    pivot_longer(-c(clinvar_class,facet))%>%
    filter(!is.na(value))%>%filter(!is.na(!!sym(facet)))
  balanced_data<-NULL
  for (tool in tools_to_plot){
    print(tool)
    for (val in unique(to_plot%>%pull(facet))){
      for (cc in c('P/LP','B/LB')){
        balanced_data<-balanced_data%>%
          bind_rows(
            to_plot%>%filter(name==tool,clinvar_class==cc)%>%
              filter(!!sym(facet)==val)%>%
              slice_sample(n=balanced_sample_size_per_tool%>%filter(name==tool,!!sym(facet)==val)%>%pull(balanced_n)))
      }
    }
  }
  if (balanced){
    to_plot<-balanced_data
  }
  # add all to the plot
  to_plot<-to_plot%>%bind_rows(
    to_plot%>%mutate(!!sym(facet):='All')
  )
  to_ret<-
    to_plot%>%
    mutate(name=stringr::str_replace(name,'_score.+','')%>%toupper())%>%
    #mutate(value_cat=cut(value,breaks=c(seq(0,1,0.05)),labels = seq(0.025,0.975,0.05)))%>%
    mutate(value_cat=cut(value,breaks=c(seq(0,1,0.1)),labels = seq(0.05,0.95,0.1)),
           name=stringr::str_replace(name,'_score.+','')%>%toupper())%>%
    group_by(name,!!sym(facet),value_cat)%>%
    count(clinvar_class)%>%
    mutate(value_cat=as.numeric(as.character(value_cat)),
           total=sum(n),
           rate=n/sum(n))%>%
    filter(clinvar_class=='P/LP')%>%
    mutate(broom::tidy(binom.test(n,total)))%>%
    ggplot(aes_string(x='value_cat',y='estimate',ymin='conf.low',ymax='conf.high',color=facet))+
    geom_line()+
    geom_point()+
    facet_wrap(name~.)+
    geom_errorbar(width=0.01)+
    ggsci::scale_color_startrek()+
    ylim(0,1)+xlim(0,1)+
    labs(color=NULL,x=NULL,y=NULL)+
    geom_abline(intercept = 0,slope = 1,linetype=2,alpha=0.5)+
    theme_minimal()+
    theme(legend.position='top',text=element_text(size=text_size))
  return(to_ret)
}
