# Takes in a varseq table with the dbnsfp scores, finds the cannonical transcript location and grabs all the scores corresponding to that 
# transcript. if there is only one score in the cell, use it. 
# if there is a different number of scores that the number of transcripts - 
# -- if its the same value, use it
# -- if its not the same value, change to NA
process_scores <- function(df, score_columns,delim=';') {
  # Applying the function to each row of the dataframe
  result <- t(apply(df, 1, function(row) {
    vep_values <- unlist(strsplit(as.character(row['vep_canonical.dbnsfp4.5a']), ";"))
    yes_indices <- which(vep_values == "YES")[1]
    if (is.na(yes_indices)){yes_indices<-1}# if there are no cannonical transcripts in the vep cannonical
    
    sapply(score_columns, function(col) {
      score_values <- unlist(strsplit(as.character(row[col]), ";"))
      
      if (length(score_values) == 1) {
        return(score_values)
      } else if (length(score_values) != length(vep_values)) {
        if (length(unique(score_values[score_values!='.'])) == 1) {
          return(score_values[score_values!='.'][1])
        }
        return(NA)
      } else {
        if (length(unique(score_values[score_values!='.'])) == 1){
          return(score_values[score_values!='.'][1])
        }
        return(score_values[yes_indices])
        # selected_value <- score_values[yes_indices]
        # return(ifelse(length(selected_value) > 0, selected_value, NA))
      }
    })
  }))
  
  # Converting the result to a data frame and setting column names
  result_df <- as.data.frame(result, stringsAsFactors = FALSE)
  colnames(result_df) <- score_columns
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
