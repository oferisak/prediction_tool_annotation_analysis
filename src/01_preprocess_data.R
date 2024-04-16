setwd('/media/SSD/Bioinformatics/Projects/prediction_tool_annotation_analysis/')
library(ProjectTemplate)
load.project()

# Setup ####
## Set analysis name #####
output_folder_name<-glue('./data/preprocessed_data/{Sys.Date()}')
if (!dir.exists(output_folder_name)){dir.create(output_folder_name)}
# Load the raw annotated clinvar data exported from VarSeq
raw_data<-readr::read_delim('./data/accessory_data/clinvar_annotated_varseq_export.20231213.tsv')
original_colnames<-colnames(raw_data)
colnames(raw_data)<-tolower(make.names(original_colnames))
## Prep clinvar submission date data ####
clinvar_submission_dates<-readr::read_delim('./data/accessory_data/clinvar_assertions.csv',col_names = c('id','date_created','date_updated'))
clinvar_min_submissions<-clinvar_submission_dates%>%group_by(id)%>%slice_min(n=1,order_by=date_created,with_ties=F)%>%ungroup()
submission_summary<-readr::read_delim('/media/SSD/Bioinformatics/Databases/clinvar/submission_summary.txt.gz',skip = 15)
submission_summary<-submission_summary%>%mutate(DateLastEvaluated=lubridate::mdy(DateLastEvaluated))
submission_summary<-submission_summary%>%select(id=`#VariationID`,DateLastEvaluated)%>%distinct()%>%group_by(id)%>%slice_min(n=1,order_by=DateLastEvaluated,with_ties=F)
## grab HG19 positions ####
## prepared using VarSeq - for the same file, just added the hg19 positions and exported only those
hg19_positions<-readr::read_delim('./data/accessory_data/hg19_positions.tsv')
colnames(hg19_positions)<-tolower(make.names(colnames(hg19_positions)))
## GOF LOF variants #### 
## downloaded from https://itanlab.shinyapps.io/goflof/
## prepared using VarSeq - for the same file, just added the annotations and exported only those
gof_lof_variants<-readr::read_delim('./data/accessory_data/goflof_ClinVar_v062021.csv')
gof_lof_variants<-gof_lof_variants%>%mutate(chr_pos_refalt=glue('{Chromosome}:{PositionVCF}:{ReferenceAlleleVCF}/{AlternateAlleleVCF}'))

# Definitions ####
## dbnsfp score columns ####
score_columns<-grep('hgvsc_vep|_score|_pred',colnames(raw_data),value = T)
## conservation score columns ####
cons_columns<-grep('phast|gerp|phylop',colnames(raw_data),value = T)
## what columns should be split ####
cols_to_fix<-c(score_columns,cons_columns)
## the Allele frequency columns ####
af_cols<-grep('^alt.allele.freq|_af',colnames(raw_data),value=T)

# Preprocessing ####
# number of missing values for each score
missing_values_before_preprocessing<-t(raw_data%>%select(cols_to_fix)%>%summarize(across(cols_to_fix,~sum(is.na(.))/n())))
fixed_score_columns<-process_scores(raw_data,cols_to_fix)
fixed_score_columns[fixed_score_columns=="."]<-NA
missing_values_after_preprocessing<-t(fixed_score_columns%>%select(cols_to_fix)%>%summarize(across(cols_to_fix,~sum(is.na(.))/n())))
score_cols<-grep('_score',score_columns,value = T)
pred_cols<-grep('_pred',score_columns,value = T)
fixed_score_columns<-fixed_score_columns%>%
  mutate(across(all_of(score_cols), as.numeric))%>%
  mutate(across(all_of(cons_columns),as.numeric))%>%
  mutate(across(all_of(pred_cols), function(x) {ifelse(grepl(",", x), NA, x)}))# some of the values are B,B or A,A - will convert them to NA
# replace the columns in the original table
proc_data<-raw_data%>%select(-cols_to_fix)%>%bind_cols(fixed_score_columns)
# add hg19 positions
proc_data<-proc_data%>%left_join(hg19_positions)
# fix clinvar and add accessory columns
proc_data<-proc_data%>%
  mutate(across(all_of(af_cols),as.numeric),
         clinvar_class=case_when(
           grepl('Pathogenic|Likely_pathogenic',classification.clinvar.2023.11.02..ncbi)~'P/LP',
           grepl('benign',classification.clinvar.2023.11.02..ncbi,ignore.case=TRUE)~'B/LB',
           grepl('conflict|uncertain',classification.clinvar.2023.11.02..ncbi,ignore.case=TRUE)~'VUS',
           grepl('not_provided',classification.clinvar.2023.11.02..ncbi,ignore.case=TRUE)~'UNKNOWN_SIGNIFICANCE',
           is.na(classification.clinvar.2023.11.02..ncbi)~'UNKNOWN_SIGNIFICANCE',
           TRUE~'other'
         ))%>%
  mutate(max_af = do.call(pmax, c(select(., all_of(af_cols)), na.rm = TRUE)),
         chr=stringr::str_replace(chr.pos.variant.info,':.+',''),
         chr_pos_refalt_hg19=glue('{chr}:{hg19_pos.1.based..dbnsfp4.5a}:{ref.alt.variant.info}'))
# add dates (there are many enteries for the same variant with different submission dates. will take the earliest)
proc_data<-proc_data%>%left_join(clinvar_min_submissions,by = c('identifier.variant.info'='id'))
proc_data<-proc_data%>%left_join(submission_summary,by=c('identifier.variant.info'='id'))
proc_data$min_date = do.call(pmin, c(proc_data[c("date_created", "DateLastEvaluated")], na.rm = TRUE))
head(proc_data%>%select(identifier.variant.info,min_date))
# add gof/lof info
gof_lof_variants%>%count(chr_pos_refalt%in%proc_data$chr_pos_refalt_hg19)
proc_data<-proc_data%>%left_join(gof_lof_variants%>%select(chr_pos_refalt,gof_lof=LABEL)%>%distinct(),
                                 by=c('chr_pos_refalt_hg19'='chr_pos_refalt'))
proc_data%>%count(gof_lof,clinvar_class)
# save the processed clinvar data
processed_file_name<-glue('{output_folder_name}/processed_data.prediction_tool_annotation.{Sys.Date()}.tsv')
write.table(proc_data,file=processed_file_name,sep='\t',row.names = F)

