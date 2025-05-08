#setwd("I:/Bin/F/LSCC.zzy/write/code")

library(dplyr)
### 
#df_file <- "I:/Bin/F/LSCC.zzy/bio_analysis/0.raw_data/P/LUAD2_expr.txt"
#expr_process(df_file)

expr_process <- function(df_file){
  
  out_dir <- "1.Signatrue"
  dir.create(out_dir,recursive = T)
  
  df <- read.table(df_file,header = TRUE,sep = "\t")
  colnames(df)[1] <- "geneSymbol"
  
  # Set the expression of repetitive proteins as their mean
  df_avg <- df %>%
    group_by(geneSymbol) %>%
    summarise(across(everything(), ~ {
      if (all(is.na(.x))) {
        NA
      } else {
        .x[is.na(.x)] <- 0
        mean(.x)
      }
    }))
  
  # Delete proteins that are missing in at least 50% of the samples
  na_rate <- rowMeans(is.na(df_avg))
  df_del <- df_avg[na_rate < 0.5,]
  gene_symbol <- df_del$geneSymbol
  
  # Fill in NA with the average expression value of protein
  df_fil <- t(apply(df_del[,-1], 1, function(x) {
    x[is.na(x)] <- mean(x, na.rm = TRUE)
    return(x)
  })) %>%
    as.data.frame() 
  row.names(df_fil) <- gene_symbol
  
  # Save data
  write.table(df_fil, paste0(out_dir,"/expression.txt"),row.names = T, sep = "\t",quote = F)
  
}




