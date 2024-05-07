setwd('/media/SSD/Bioinformatics/Projects/prediction_tool_annotation_analysis/')
library(ProjectTemplate)
load.project()

# Set Var Sets ####
# should be ran only after running the preprocessing script or loading the saved preprocessed data table
## Definitions ####
data_folder_name<-'./data/preprocessed_data/2024-01-17'
proc_data_file<-grep('processed_data.prediction_tool_annotation',list.files(data_folder_name,full.names = T),value=T)
proc_data<-readr::read_delim(proc_data_file)
### Allele frequency threholds to split on ####
af_threshs<-10^c(-2:-6)
conservation_column<-'phylop100way_vertebrate.dbnsfp4.5a'
num_o_conservation_levels<-3
## The main var sets ####
## these var sets will be further subdivided into the different sub_var_sets
main_var_sets<-list('confident_all'=grepl('(criteria_provided,_multiple_submitters,_no_conflicts)|(criteria_provided,_single_submitter)|(reviewed_by_expert_panel)',proc_data$clinvar_review.dbnsfp4.5a))

# missing values filter ####
tools_list<-grep('_score',colnames(proc_data),value = T)
tools_with_pred<-grep('_pred',colnames(proc_data),value = T)

tools_list<-c(tools_list[which(stringr::str_replace_all(tools_list,'_score','_pred')%in%tools_with_pred)],'revel_score.dbnsfp4.5a')
              


# Excluded tools from missingness calc (because they have to few variants)
high_missingness_tools<-c('eve_score.dbnsfp4.5a',
                          'mutpred_score.dbnsfp4.5a',
                          'esm1b_score.dbnsfp4.5a',
                          'list.s2_score.dbnsfp4.5a')
missing_var_sets<-list()
missing_var_sets[['full']]<-1:nrow(proc_data)
missing_var_sets[['no_missing']]<-apply(proc_data[setdiff(tools_list,high_missingness_tools)], 1, function(row) all(!is.na(row)))


## The sub var sets ####
sub_var_sets<-list('all'=rep(TRUE,nrow(proc_data)),
                   'AD'=proc_data$inheritance.clinical.genomic.database.2023.05.04..ghi=='Autosomal Dominant',
                   'AR'=proc_data$inheritance.clinical.genomic.database.2023.05.04..ghi=='Autosomal Recessive')

sub_var_sets[['multiple_submitters']]<-grepl('(criteria_provided,_multiple_submitters,_no_conflicts)|(reviewed_by_expert_panel)',proc_data$clinvar_review.dbnsfp4.5a)
# Per year var set
years<-2013:2023
date_col<-'date_created'
for (year in years){
  sub_var_sets[[as.character(year)]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))==year
  sub_var_sets[[glue('at_or_after_{year}')]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))>=year
  sub_var_sets[[glue('before_{year}')]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))<year
}
sub_var_sets[['2013']]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))<=2013



### Allele frequency ####
sub_var_sets[['af_unknown']]<-is.na(proc_data$max_af)
for (af_thresh in af_threshs){
  sub_var_sets[[glue('af_above_{af_thresh}')]]<-(!is.na(proc_data$max_af))&(proc_data$max_af>af_thresh)
}

### GOF / LOF #### 
sub_var_sets[['gof']]=(!is.na(proc_data$gof_lof))&(proc_data$gof_lof=='GOF')
sub_var_sets[['lof']]=(!is.na(proc_data$gof_lof))&(proc_data$gof_lof=='LOF')

### Allele frequency cat ####
proc_data<-proc_data%>%mutate(max_af_cat=Hmisc::cut2(max_af,g = 5))
proc_data%>%count(max_af_cat)
sub_var_sets[[glue('af_cat_NA')]]<-is.na(proc_data$max_af_cat)
for (af_level in levels(proc_data$max_af_cat)){
  sub_var_sets[[glue('af_cat_{af_level}')]]<-proc_data$max_af_cat==af_level
}

### Conservation ####
proc_data<-proc_data%>%mutate(cons_cat=Hmisc::cut2(!!sym(conservation_column),g = num_o_conservation_levels))
if (num_o_conservation_levels==3){
  levels(proc_data$cons_cat)<-c('low','intermediate','high')
}
for (cons_level in levels(proc_data$cons_cat)){
  sub_var_sets[[glue('{cons_level}_conservation')]]<-proc_data$cons_cat==cons_level
}

## PANELAPP annotations ####
panelapp_data<-readr::read_delim('./data/accessory_data/clinvar_panelapp.tsv',delim='\t')
colnames(panelapp_data)<-tolower(make.names(colnames(panelapp_data)))
### Disease groups ####
panelapp_disease_group<-panelapp_data%>%select(contains(c('.info','disease.group')))%>%
  separate_rows(disease.group.panel.app.genomic.genes.2023.11.01..ge,sep=',')%>%
  filter(disease.group.panel.app.genomic.genes.2023.11.01..ge!='')%>%distinct()
for (disease_group in unique(panelapp_disease_group$disease.group.panel.app.genomic.genes.2023.11.01..ge)){
  message(disease_group)
  disease_group_vars<-proc_data%>%
    left_join(panelapp_disease_group%>%filter(disease.group.panel.app.genomic.genes.2023.11.01..ge==disease_group))
  sub_var_sets[[glue('panelapp_disease_group:{tolower(make.names(disease_group))}')]]<-!is.na(disease_group_vars$disease.group.panel.app.genomic.genes.2023.11.01..ge)
}
### Panels ####
# panelapp_panels<-panelapp_data%>%select(contains(c('.info','panel.name')))%>%
#   separate_rows(panel.name.panel.app.genomic.genes.2023.11.01..ge,sep=',')%>%
#   filter(panel.name.panel.app.genomic.genes.2023.11.01..ge!='')%>%distinct()
# 
# for (panel_name in unique(panelapp_panels$panel.name.panel.app.genomic.genes.2023.11.01..ge)){
#   message(panel_name)
#   panel_vars<-proc_data%>%
#     left_join(panelapp_panels%>%filter(panel.name.panel.app.genomic.genes.2023.11.01..ge==panel_name))
#   sub_var_sets[[glue('panelapp_panel:{tolower(make.names(panel_name))}')]]<-!is.na(panel_vars$panel.name.panel.app.genomic.genes.2023.11.01..ge)
# }
# 
# # if there are less than 2000 variants in the panel remove it
# for (panelapp_group in grep('panelapp',names(sub_var_sets),value=T)){
#   if (sum(sub_var_sets[[panelapp_group]])<2000){sub_var_sets[[panelapp_group]]<-NULL}
# }

save(main_var_sets,sub_var_sets,missing_var_sets,file=glue('{data_folder_name}/var_sets.RData'))
