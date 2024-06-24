setwd('/media/SSD/Bioinformatics/Projects/prediction_tool_annotation_analysis/')
library(ProjectTemplate)
load.project()

# Setup ####
## Set analysis name #####
output_folder_name<-glue('./data/preprocessed_data/{Sys.Date()}')
if (!dir.exists(output_folder_name)){dir.create(output_folder_name)}
# Load the raw annotated clinvar data exported from VarSeq
#raw_data<-readr::read_delim('./data/accessory_data/clinvar_annotated_varseq_export.20231213.tsv')
raw_data<-readr::read_delim('./data/prediction_tools_analysis_hg38_20240610.tsv')
original_colnames<-colnames(raw_data)
colnames(raw_data)<-tolower(make.names(original_colnames))
# add cadd phred to the scores
raw_data$cadd_phred_score.dbnsfp4.5a<-raw_data$cadd_phred.dbnsfp4.5a

# Definitions ####
## dbnsfp score columns ####
score_names<-stringr::str_match(grep('_score',colnames(raw_data),value=T),'(.+)_score')[,2]
cons_columns<-grep('rankscore',grep('phast|gerp|phylop',colnames(raw_data),value = T),invert = T,value=T)
all_score_columns<-grep('_score|_pred|_rankscore',colnames(raw_data),value = T)
all_score_columns<-c(all_score_columns,cons_columns)
numeric_columns_to_fix<-unique(c(grep('_score|_rankscore',colnames(raw_data),value = T),cons_columns))
pred_columns_to_fix<-c(grep('_pred',colnames(raw_data),value = T))
clinvar_classification_col<-'clnsig.variant.info'

# Preprocessing ####
# remove all variants without dbnsfp annotation
raw_data<-raw_data%>%filter(!is.na(hgvsc_vep.dbnsfp4.5a))
# make sure the clinvar ids match between the clinvar annotation and the original vcf variant (so that when i grab the hgvsc annotation it corresponds to the correct variant)
# all the variants that do not match are those without clinvar annotation (because it is an earlier version)
raw_data<-raw_data%>%filter(identifier.variant.info==variant.id.clinvar.2024.04.04..ncbi)

# grab cds from clinvar hgvs
proc_data<-raw_data%>%
  mutate(clinvar_annotation_cds=stringr::str_match(hgvs.c..name.clinvar.2024.04.04..ncbi,'NM_[\\d\\.]+:(c\\.[^\\;]+)')[,2])

# check which scores need to be converted
converted_scores<-stringr::str_match(numeric_columns_to_fix,'(.+)_converted')[,2]
converted_scores<-glue('{converted_scores[!is.na(converted_scores)]}_score.dbnsfp4.5a')

# fix scores by taking the most damaging prediction for each score
fixed_numeric_scores<-process_numeric_scores_by_most_damaging(proc_data,numeric_columns_to_fix,delim = ';',converted_scores)
fixed_pred_scores<-process_pred_scores_by_most_damaging(proc_data,pred_columns_to_fix,delim = ';')
fixed_score_columns<-bind_cols(fixed_numeric_scores,fixed_pred_scores)

# replace the columns in the original table
proc_data<-proc_data%>%select(-all_score_columns)%>%bind_cols(fixed_score_columns)

# number of missing values for each score
missing_values_before_preprocessing<-t(raw_data%>%select(all_score_columns)%>%summarize(across(all_score_columns,~sum(is.na(.))/n())))
missing_values_before_preprocessing<-missing_values_before_preprocessing%>%as.data.frame()%>%mutate(tool=rownames(missing_values_before_preprocessing))%>%rename(before=V1)
missing_values_after_preprocessing<-t(proc_data%>%summarize(across(all_score_columns,~sum(is.na(.))/n())))
missing_values_after_preprocessing<-missing_values_after_preprocessing%>%as.data.frame()%>%mutate(tool=rownames(missing_values_after_preprocessing))%>%rename(after=V1)
comp_missing<-missing_values_before_preprocessing%>%left_join(missing_values_after_preprocessing,by='tool')%>%
  mutate(diff=after-before)

# look at one conversion for example
tool_to_comp<-'eve_score.dbnsfp4.5a'
tool_comp<-raw_data%>%select(hgvs.c..name.clinvar.2024.04.04..ncbi,before=tool_to_comp)%>%bind_cols(proc_data%>%select(vep_canonical.dbnsfp4.5a,hgvsc_vep.dbnsfp4.5a,clinvar_annotation_cds,after=tool_to_comp))


# fix clinvar and add accessory columns
proc_data<-proc_data%>%
  mutate(clinvar_class=case_when(
           grepl('Pathogenic|Likely_pathogenic',!!sym(clinvar_classification_col))~'P/LP',
           grepl('benign',!!sym(clinvar_classification_col),ignore.case=TRUE)~'B/LB',
           grepl('conflict|uncertain',!!sym(clinvar_classification_col),ignore.case=TRUE)~'VUS',
           grepl('not_provided',!!sym(clinvar_classification_col),ignore.case=TRUE)~'UNKNOWN_SIGNIFICANCE',
           is.na(!!sym(clinvar_classification_col))~'UNKNOWN_SIGNIFICANCE',
           TRUE~'other'
         ))%>%
  mutate(chr=stringr::str_replace(chr.pos.variant.info,':.+',''),
         chr_pos_refalt_hg19=glue('{chr}:{hg19_pos.1.based..dbnsfp4.5a}:{ref.alt.variant.info}'))
# set min-date
# load the preprocessed VCV file
# downloading https://ftp.ncbi.nlm.nih.gov/pub/clinvar/xml/ClinVarVCVRelease_00-latest.xml.gz
# there is a VariationArchive element with a variable called DateCreated which seems to correspond to the First-created date
# generated it using: 
#   zcat ClinVarVCVRelease_2024-05.xml.gz | grep "<VariationArchive" > ClinVarVCVRelease_2024-05.VariationArchive.txt
# and then converted into a csv using:
#   grep 'VariationArchive' ClinVarVCVRelease_2024-05.VariationArchive.txt | awk -F'"' '{print $4"\t"$18"\t"$20"\t"$22"\t"$6}' > ClinVarVCVRelease_2024-05.VariationArchive.tsv
# the date_created date corresponds to the First in Clinvar date in the variant page
vcv_file<-readr::read_delim('./data/accessory_data/ClinVarVCVRelease_2024-05.VariationArchive.tsv',col_names = c('variant_id','date_last_updated','date_created','most_recent_submission','variant_name'))
proc_data<-proc_data%>%left_join(vcv_file,by=c('identifier.variant.info'='variant_id'))

# Mismatch constraint ####
## gnomAD ####
proc_data$missense_min_oe_ratio<-apply(proc_data,1,function(row){
  mmoe<-as.numeric(unlist(strsplit(as.character(row['missense.variants...observed.expected.gnomad...gene.constraint.4.0..broad']), ",")))
  return(min(mmoe,na.rm = T))
})
proc_data$missense_min_oe_ratio[is.infinite(proc_data$missense_min_oe_ratio)]<-NA

## Regeneron ####
regeneron_mc<-readr::read_delim('./data/accessory_data/regeneron_million_exomes.media-1.txt')
regeneron_mc<-regeneron_mc%>%mutate(regeneron_constraint=Hmisc::cut2(mean,g=10))
proc_data<-proc_data%>%left_join(regeneron_mc%>%select(GeneName,regeneron_constraint,shet_constrained),by=c('gene.names.clinvar.2024.04.04..ncbi'='GeneName'))

# Mode of inheritance ####
moi_gencc<-apply(proc_data,1,function(row){
  mois<-unique(unlist(strsplit(as.character(row['mode.of.inheritance.gencc.db.2024.02.27..gencc']), ",")))
  if (length(mois)==1){
    return(mois)
  }else{return('undertermined')}
})
proc_data$moi_gencc<-moi_gencc

# pLI and pRec - ####
# pLI - probability of being LOF intolerant
# pRec - probability of being intolerant to homozygous but not heterozygous variants
proc_data$prec_gnomad<-apply(proc_data,1,function(row){
  prec_gnomad<-as.numeric(unlist(strsplit(as.character(row['prec.gnomad...gene.constraint.4.0..broad']), ",")))
  return(max(prec_gnomad,na.rm = T))
})
proc_data$prec_gnomad[is.infinite(proc_data$prec_gnomad)]<-NA

proc_data$pli_gnomad<-apply(proc_data,1,function(row){
  pli_gnomad<-as.numeric(unlist(strsplit(as.character(row['pli.gnomad...gene.constraint.4.0..broad']), ",")))
  return(max(pli_gnomad,na.rm = T))
})
proc_data$pli_gnomad[is.infinite(proc_data$pli_gnomad)]<-NA
proc_data$pli_gnomad_cat<-ifelse(proc_data$pli_gnomad>=0.9,'pli_high','pli_low')
proc_data%>%count(pli_gnomad_cat)
proc_data$prec_gnomad_cat<-ifelse(proc_data$prec_gnomad>=0.9,'prec_high','prec_low')
proc_data%>%count(prec_gnomad_cat)

# add gene ontology info
min_plp_blb_counts<-2000

# add processes
gene2go_process<-readr::read_delim('/media/SSD/Bioinformatics/Projects/ncbi_utils_2024/output/gene2go_process.csv')
proc_data<-proc_data%>%left_join(gene2go_process,by=c('gene.names.clinvar.2024.04.04..ncbi'='gene_symbol'))
# now change every missing value to false (if its missing it means the gene is not in the GO value)
proc_data<-proc_data%>%mutate(across(all_of(setdiff(colnames(gene2go_process),'gene_symbol')), ~replace_na(., FALSE)))

counts<-NULL
for (colname in colnames(gene2go_process%>%select(-gene_symbol))){
  print(colname)
  if (nrow(proc_data%>%filter(!!sym(colname))%>%count(clinvar_class)%>%as.data.frame())==0){
    counts<-counts%>%bind_rows(data.frame(colname,B.LB=0,P.LP=0))
  }
  counts<-counts%>%bind_rows(
    data.frame(colname,proc_data%>%filter(!!sym(colname))%>%count(clinvar_class)%>%as.data.frame()%>%pivot_wider(names_from = clinvar_class,values_from = n))
  )
}
low_count_processes<-counts%>%filter(B.LB<min_plp_blb_counts | P.LP<min_plp_blb_counts)%>%pull(colname)
proc_data<-proc_data%>%select(-low_count_processes)

# add functions
gene2go_function<-readr::read_delim('/media/SSD/Bioinformatics/Projects/ncbi_utils_2024/output/gene2go_function.csv')
proc_data<-proc_data%>%left_join(gene2go_function,by=c('gene.names.clinvar.2024.04.04..ncbi'='gene_symbol'))
# now change every missing value to false (if its missing it means the gene is not in the GO value)
proc_data<-proc_data%>%mutate(across(all_of(setdiff(colnames(gene2go_function),'gene_symbol')), ~replace_na(., FALSE)))
counts<-NULL
for (colname in colnames(gene2go_function%>%select(-gene_symbol))){
  print(colname)
  if (nrow(proc_data%>%filter(!!sym(colname))%>%count(clinvar_class)%>%as.data.frame())==0){
    counts<-counts%>%bind_rows(data.frame(colname,B.LB=0,P.LP=0))
  }else{
    counts<-counts%>%bind_rows(
      data.frame(colname,proc_data%>%filter(!!sym(colname))%>%count(clinvar_class)%>%as.data.frame()%>%pivot_wider(names_from = clinvar_class,values_from = n))
    )
    }
}
low_count_functions<-counts%>%filter(B.LB<min_plp_blb_counts | P.LP<min_plp_blb_counts)%>%pull(colname)
proc_data<-proc_data%>%select(-low_count_functions)

# if you want to remove existing process and functions
#for(go_proccess in setdiff(colnames(gene2go_process),'gene_symbol')){if(go_proccess%in%colnames(proc_data)){proc_data<-proc_data%>%select(-go_proccess)}}
#for(go_function in setdiff(colnames(gene2go_function),'gene_symbol')){if(go_function%in%colnames(proc_data)){proc_data<-proc_data%>%select(-go_function)}}

# add primate-ai3d scores (hg38)
# according to https://primad.basespace.illumina.com/help - 0.8 should be considered the cutoff for pathogenicity
primateai3d_table<-readr::read_delim('/media/SSD/Bioinformatics/Databases/primateai3d/PrimateAI-3D_scores (1).csv.gz',delim=',')
primateai3d_table<-primateai3d_table%>%rename(primateai3d_score.dbnsfp4.5a=score_PAI3D,
                                              reference.variant.info=non_flipped_ref,
                                              alternates.variant.info=non_flipped_alt)%>%
  mutate(primateai3d_pred.dbnsfp4.5a=ifelse(primateai3d_score.dbnsfp4.5a>0.8,'D','T'),
         chr.pos.variant.info=glue('{stringr::str_replace(chr,"chr","")}:{pos}'))

head(primateai3d_table)
proc_data<-proc_data%>%left_join(primateai3d_table%>%select(chr.pos.variant.info,
                                                            reference.variant.info,
                                                            alternates.variant.info,
                                                            primateai3d_score.dbnsfp4.5a,
                                                            primateai3d_pred.dbnsfp4.5a))


# generate panelapp gene-panel table
panelapp_raw<-proc_data%>%select(contains('panel.app'))
panelapp_raw<-panelapp_raw%>%mutate(panel.name.panel.app.genomic.genes.2023.11.01..ge=stringr::str_replace(panel.name.panel.app.genomic.genes.2023.11.01..ge,'von Willebrand','Von Willebrand'),
                                    panel.name.panel.app.genomic.genes.2023.11.01..ge=stringr::str_replace_all(panel.name.panel.app.genomic.genes.2023.11.01..ge,",(?![A-Z])", ""))
panelapp_long<-panelapp_raw%>%separate_longer_delim(colnames(panelapp_raw)[1:4],delim = ',')
# remove all Red associations
panelapp_long<-panelapp_long%>%
  rename(gene_symbol=gene.name.panel.app.genomic.genes.2023.11.01..ge,
         panel_name=panel.name.panel.app.genomic.genes.2023.11.01..ge)%>%
  filter(!confidence.level.panel.app.genomic.genes.2023.11.01..ge%in%c('Red','No List'))%>%distinct()
# remove panels that include less than 5 genes
panel_counts<-panelapp_long%>%select(panel_name,gene_symbol)%>%distinct()%>%
  count(panel_name)
excluded_panels<-panel_counts%>%filter(n<10|n>1000)%>%pull(panel_name)
excluded_panels<-c(excluded_panels,grep('COVID|Additional',panel_counts$panel_name,value = T))

panelapp_wide<-panelapp_long%>%
  filter(!panel_name%in%excluded_panels)%>%
  select(gene_symbol,panel_name)%>%
  mutate(value=TRUE)%>%
  distinct()%>%
  pivot_wider(names_from = panel_name,values_from = value,values_fill = list(value=FALSE))
colnames(panelapp_wide)<-paste0('panelapp:',make.names(colnames(panelapp_wide)))
colnames(panelapp_wide)[1]<-'gene_symbol'
proc_data<-proc_data%>%left_join(panelapp_wide,by=c('gene.names.clinvar.2024.04.04..ncbi'='gene_symbol'))

count_clinvar_class <- function(data, panel_names) {
  result <- NULL
  
  for (panel in panel_names) {
    counts <- data %>%
      filter(get(panel) == TRUE) %>%
      group_by(clinvar_class) %>%
      summarise(count = n()) %>%
      pivot_wider(names_from = clinvar_class, values_from = count, values_fill = list(count = 0))%>%
      mutate(panel_name=panel)
    
    result <- bind_rows(result, counts)
  }
  
  result[is.na(result)] <- 0
  return(result)
}
# keep only panels that have at least 750 variants for B and for P
panel_clinvar_counts<-count_clinvar_class(proc_data,panel_names = setdiff(colnames(panelapp_wide),'gene_symbol'))
panels_to_remove<-panel_clinvar_counts%>%filter(`B/LB`<750 | `P/LP`<750 )%>%pull(panel_name)
proc_data<-proc_data%>%select(-panels_to_remove)

# save the processed clinvar data
processed_file_name<-glue('{output_folder_name}/processed_data.prediction_tool_annotation.{Sys.Date()}.tsv')
write.table(proc_data,file=processed_file_name,sep='\t',row.names = F)
