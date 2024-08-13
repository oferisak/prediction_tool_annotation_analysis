setwd('/media/SSD/Bioinformatics/Projects/prediction_tool_annotation_analysis/')
library(ProjectTemplate)
load.project()

# Set Var Sets ####
# should be ran only after running the preprocessing script or loading the saved preprocessed data table
## Definitions ####
data_folder_name<-'./data/preprocessed_data/2024-06-10/'
proc_data_file<-grep('processed_data.prediction_tool_annotation',list.files(data_folder_name,full.names = T),value=T)
proc_data<-readr::read_delim(proc_data_file)
### Allele frequency threholds to split on ####
af_threshs<-10^c(-2:-6)
af_col<-'allele.frequency.gnomad.joint.variant.frequencies.4.0.v2..broad'
conservation_column<-'phylop100way_vertebrate.dbnsfp4.5a'
num_o_conservation_levels<-3
tools_list<-grep('_score',colnames(proc_data),value = T)
tools_with_pred<-grep('_pred',colnames(proc_data),value = T)

tools_list<-c(tools_list[which(stringr::str_replace_all(tools_list,'_score','_pred')%in%tools_with_pred)],'revel_score.dbnsfp4.5a')
## The main var sets ####
var_sets<-list('all'=rep(TRUE,nrow(proc_data)),
               'high_confidence'=!grepl('criteria_provided,_single_submitter',proc_data$clnrevstat.variant.info))

# missing values filter ####
var_sets[['AD']]<-grepl('Autosomal dominant',proc_data$moi_gencc)
var_sets[['AR']]<-grepl('Autosomal recessive',proc_data$moi_gencc)

# Per year var set
years<-2013:2024
date_col<-'date_created'
for (year in years){
  var_sets[[as.character(year)]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))==year
  var_sets[[glue('at_or_after_{year}')]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))>=year
  var_sets[[glue('before_{year}')]]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))<year
}
var_sets[['2013']]<-!is.na(proc_data%>%pull(date_col)) & year(proc_data%>%pull(date_col))<=2013

### Allele frequency ####
var_sets[['af_unknown']]<-is.na(proc_data%>%pull(af_col))
for (af_thresh in af_threshs){
  var_sets[[glue('af_above_{af_thresh}')]]<-(!is.na(proc_data%>%pull(af_col)))&(proc_data%>%pull(af_col)>af_thresh)
}
### Allele frequency cat ####
proc_data<-proc_data%>%mutate(af_cat=cut(proc_data%>%pull(af_col),breaks = c(1,0.001,1e-4,1e-5,1e-6,0)))
proc_data%>%count(af_cat)

var_sets[[glue('af_cat_unknown')]]<-is.na(proc_data$af_cat)
for (af_level in levels(proc_data$af_cat)){
  var_sets[[glue('af_cat_{af_level}')]]<-(!is.na(proc_data$af_cat)) & proc_data$af_cat==af_level
}

### Conservation ####
proc_data<-proc_data%>%mutate(cons_cat=Hmisc::cut2(!!sym(conservation_column),g = num_o_conservation_levels))
if (num_o_conservation_levels==3){
  levels(proc_data$cons_cat)<-c('low','intermediate','high')
}
for (cons_level in levels(proc_data$cons_cat)){
  var_sets[[glue('{cons_level}_conservation')]]<-proc_data$cons_cat==cons_level
}

### Haploinsufficiency/Triplosensitivity
var_sets[['haploinsufficiency_genes']]<-proc_data%>%pull(haploinsufficiency.description.clingen.gene.dosage.sensitivity.2024.04.01..ncbi)=='Sufficient evidence for dosage pathogenicity' & !is.na(proc_data%>%pull(haploinsufficiency.description.clingen.gene.dosage.sensitivity.2024.04.01..ncbi))
var_sets[['triplosensitivity_regions']]<-grepl('Sufficient evidence for dosage pathogenicity',proc_data$triplosensitivity.description.clingen.region.dosage.sensitivity.2024.03.01..ncbi)

### Missense constraint (min O/E ratio)
proc_data<-proc_data%>%mutate(missense_min_oe_ratio_cat=Hmisc::cut2(missense_min_oe_ratio,g = 10,levels.mean = T))
for (mm_oe_level in levels(proc_data$missense_min_oe_ratio_cat)){
  var_sets[[glue('mismatch_min_oe_ratio:{mm_oe_level}')]]<-(!is.na(proc_data$missense_min_oe_ratio_cat)) & proc_data$missense_min_oe_ratio_cat==mm_oe_level
}
proc_data<-proc_data%>%mutate(missense_min_oe_ratio_manual=cut(missense_min_oe_ratio,breaks=seq(0,1.5,0.2),labels = seq(0.1,1.4,0.2)))
for (mm_oe_level in levels(proc_data$missense_min_oe_ratio_manual)){
  var_sets[[glue('mismatch_min_oe_ratio_manual:{mm_oe_level}')]]<-(!is.na(proc_data$missense_min_oe_ratio_manual)) & proc_data$missense_min_oe_ratio_manual==mm_oe_level
}

# Missense constraint Regeneron
for (regeneron_level in levels(proc_data$regeneron_constraint)){
  var_sets[[glue('regeneron_mc:{regeneron_level}')]]<-(!is.na(proc_data$regeneron_constraint)) & proc_data$regeneron_constraint==regeneron_level
}
var_sets[['regeneron_constrained']]<-(!is.na(proc_data$shet_constrained)) & proc_data$shet_constrained
var_sets[['regeneron_not_constrained']]<-(!is.na(proc_data$shet_constrained)) & !proc_data$shet_constrained

# pLI and pRec

var_sets[['pli_high']]<-proc_data$pli_gnomad_cat=='pli_high'
var_sets[['prec_high']]<-proc_data$prec_gnomad_cat=='prec_high'

save(var_sets,file=glue('{data_folder_name}/var_sets.RData'))
