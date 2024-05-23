process_numeric_scores_by_most_damaging<-function(df,numeric_score_cols,delim=';',converted_scores){
  # first replace all ',' with ';' - not sure its needed, the added value is very minor and probably false negatives are entered in the mix
  # df[numeric_score_cols] <- lapply(df[numeric_score_cols], function(x) gsub(",", ";", x))
  
  result <- t(apply(df, 1, function(row) {
    sapply(numeric_score_cols, function(col) {
      score_values <- as.numeric(unlist(strsplit(as.character(row[col]), ";")))
      if (all(is.na(score_values))){return(NA)}
      if (col%in%converted_scores){
        return(min(score_values,na.rm = T))
      }else{
        return(max(score_values,na.rm = T))
      }
    })
  }))
  
  # Converting the result to a data frame and setting column names
  result_df <- as.data.frame(result, stringsAsFactors = FALSE)
  return(result_df)
}

process_pred_scores_by_most_damaging<-function(df,pred_score_cols,delim=';'){
  df[pred_score_cols] <- lapply(df[pred_score_cols], function(x) gsub(",", ";", x))
  
  result <- t(apply(df, 1, function(row) {
    sapply(pred_score_cols, function(col) {
      score_values <- unlist(strsplit(as.character(row[col]), ";"))
      score_values[score_values=='.']<-NA
      score_values[grepl(',',score_values)]<-NA
      if (all(is.na(score_values))){return(NA)}
      to_ret<-
      case_when(
        'D'%in%score_values~'D',
        'P'%in%score_values & grepl('polyphen2|alphamissense',col,ignore.case = T)~'D',
        'A'%in%score_values~'D',
        (('H'%in%score_values)|('M'%in%score_values))&grepl('mutationassessor',col,ignore.case = T)~'D',
        .default = 'T'
      )
      return(to_ret)
    })
  }))
  
  # Converting the result to a data frame and setting column names
  result_df <- as.data.frame(result, stringsAsFactors = FALSE)
  return(result_df)
}

extract_unique_values <- function(df) {
  # List to store results
  results <- list()
  
  # Loop through each column of the dataframe
  for (col_name in names(df)) {
    # Combine all strings in the column, split them by ';', and find unique elements
    unique_values <- unique(unlist(strsplit(unlist(df[[col_name]]), ";")))
    unique_values[unique_values=='.']<-NA
    unique_values[grepl(',',unique_values)]<-NA
    # Add the result to the list with the column name as the key
    results[[col_name]] <- unique(unique_values)
  }
  
  # Return the list of results
  return(results)
}

#extract_unique_values(raw_data%>%select(pred_columns_to_fix))

process_numeric_scores_by_cds<-function(df,numeric_score_cols,delim=';',converted_scores){
  result <- t(apply(df, 1, function(row) {
    clinvar_cds<-row[['clinvar_annotation_cds']]
    dbnsfp_cds_values <- unlist(strsplit(row[['hgvsc_vep.dbnsfp4.5a']], ";"))

    matching_cds_values <- which(dbnsfp_cds_values == clinvar_cds)
    
    sapply(numeric_score_cols, function(col) {
      score_values <- as.numeric(unlist(strsplit(as.character(row[col]), ";")))
      if (all(is.na(score_values))){return(NA)}
      if (length(score_values)==length(dbnsfp_cds_values)){ #if the length of the scores matches the number of cds options in dbnsfp
        if (length(matching_cds_values)==1){# if there is matching cds - use it (even if it NA)
          return(score_values[matching_cds_values])
        }
        if (length(matching_cds_values)==0){# if there is matching cds - use it (even if it NA)
          return(NA)
        }
        # otherwise - take the max / min out of the matching 
        if (col%in%converted_scores){
          return(min(score_values[matching_cds_values],na.rm = T))
        }else{
          return(max(score_values[matching_cds_values],na.rm = T))
        }
      }
      # if there is no clinvar cds defined, or there is a matching cds and there is only one unique value, use it:
      if ((is.na(clinvar_cds)|length(matching_cds_values)>0) & length(unique(score_values[!is.na(score_values)]))==1){
        return(unique(score_values[!is.na(score_values)]))
      }else{# otherwise, it means that either there is no matching cds or there is more than one value
        return(NA)
      }
      

    })
  }))
  
  # Converting the result to a data frame and setting column names
  result_df <- as.data.frame(result, stringsAsFactors = FALSE)
  return(result_df)
}

find_multivalue_list_columns <- function(df) {
  multivalue_columns <- sapply(df, function(column) {
    if (is.list(column)) {
      any(sapply(column, function(item) length(item) > 1))
    } else {
      FALSE
    }
  })
  names(multivalue_columns)[multivalue_columns]
}

process_pred_scores_by_cds<-function(df,pred_score_cols,delim=';'){
  result <- t(apply(df, 1, function(row) {
    clinvar_cds<-row[['clinvar_annotation_cds']]
    dbnsfp_cds_values <- unlist(strsplit(row[['hgvsc_vep.dbnsfp4.5a']], ";"))
    
    matching_cds_values <- which(dbnsfp_cds_values == clinvar_cds)
    
    sapply(pred_score_cols, function(col) {
      score_values <- unlist(strsplit(as.character(row[col]), ";"))
      score_values[score_values=='.']<-NA
      if (all(is.na(score_values))){return(NA)}
      if (length(score_values)==length(dbnsfp_cds_values)){ #if the length of the scores matches the number of cds options in dbnsfp
        if (length(matching_cds_values)==1){# if there is matching cds - use it (even if it NA)
          return(score_values[matching_cds_values])
        }
        if (length(matching_cds_values)==0){# if there is no matching cds - return NA
          return(NA)
        }
        if (all(is.na(score_values[matching_cds_values]))){return(NA)}# if all the matching cds scores are NA return NA
        # otherwise - return the most common pred of the matching values
        most_common_pred <- names(which.max(table(score_values[matching_cds_values])))
        return(most_common_pred)
      }
      # if (there is no clinvar cds defined, or there is a matching cds) and there is only one unique value, use it:
      if ((is.na(clinvar_cds)|length(matching_cds_values)>0) & length(unique(score_values[!is.na(score_values)]))==1){
        return(unique(score_values[!is.na(score_values)]))
      }else{# otherwise, it means that either there is no matching cds or there is more than one value
        return(NA)
      }
      
      
    })
  }))
  
  # Converting the result to a data frame and setting column names
  result_df <- as.data.frame(result, stringsAsFactors = FALSE)
  return(result_df)
}




count_unique_scores <- function(df, score_columns) {
  for (col in score_columns[grepl('_score',score_columns)]) {
    message(col)
    count_col_name <- paste(col, "count", sep = "_")
    df[[count_col_name]] <- sapply(df[[col]], function(cell) {
      # Split the cell by ';', remove '.', and count unique values
      length(unique(na.omit(unlist(strsplit(cell, ";"))[unlist(strsplit(cell, ";")) != '.'])))
    })
  }
  return(df)
}
#df<-count_unique_scores(df,score_columns)
# count_cols<-grep('_count',colnames(df),value = T)
# df%>%summarize(across(count_cols,~sum(.==1)/n()))
