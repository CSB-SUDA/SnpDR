
library(impute)
library(pROC)
library(ggplot2)
library(cowplot)

sensPredict <- function(node_file,cell="LUSC"){
  out_dir <- "result/4.Prediction/A.sensitive"
  dir.create(out_dir,recursive = T)
  
  # Import proteins
  node_df <- read.table(node_file,header = T,sep = "\t",col.names = "node")
  
  ctrp_exp <- read.table("data/Signatrue/Pharmacogenomic/CTRP/data_exp.txt",header = T,sep = "\t",check.names = F)
  ctrp_act <- read.table("data/Signatrue/Pharmacogenomic/CTRP/data_act.txt",header = T,sep = "\t",check.names = F)
  cell_line <- read.table("data/Signatrue/Pharmacogenomic/CTRP/cell_line_annotation.txt",header = T,sep = "\t",quote="\"",check.names = F)
  cell_line_sel <- cell_line[grepl(cell,cell_line$OncoTree3),]
  
  # Filter drug activity data
  row.names(ctrp_act) <- ctrp_act$ID
  drug_mat <- ctrp_act[,cell_line_sel$Name]
  # Filter drugs with activity scores exceeding 50% in cell line loss
  drug_mat <- drug_mat[rowSums(is.na(drug_mat))<(0.5*ncol(drug_mat)),]
  # Filter cell lines with activity scores exceeding 80% drug deficiency
  drug_mat <- drug_mat[,colSums(is.na(drug_mat))<(0.8*nrow(drug_mat))]
  # Estimating missing values using KNN method
  drug_mat_imp <- impute.knn(as.matrix(drug_mat))$data
  
  # Filter expression data
  row.names(ctrp_exp) <- ctrp_exp$ID
  expr_mat <- ctrp_exp[,colnames(drug_mat_imp)]
  comm_node <- intersect(node_df$node,row.names(expr_mat))
  expr_mat <- as.matrix(expr_mat[comm_node,])
  
  results <- c()
  
  # Traverse each drug
  for (drug in rownames(drug_mat_imp)) {
    y <- drug_mat_imp[drug, ]
    
    # Initialize the analysis results of each gene
    gene_result <- data.frame(Gene=character(),
                              AUC=numeric(),
                              ANOVA_P=numeric(),
                              Label=character(),
                              stringsAsFactors = FALSE)
    
    for (gene in rownames(expr_mat)){
      x <- expr_mat[gene,]
      
      # Constructing a single factor linear regression model
      fit <- lm(y ~ x)
      aov_p <- anova(fit)$`Pr(>F)`[1]
      
      # Calculate AUC based on median drug activity grouping
      group <- ifelse(y > median(y), "High", "Low")
      roc_obj <- tryCatch({
        roc(group, x)
      }, error=function(e) NA)
      
      auc_val <- if (!is.na(roc_obj)[1]) as.numeric(auc(roc_obj)) else NA
      status <- "NS"
      if (!is.na(auc_val) && aov_p < 0.05) {
        if (auc_val >= 0.75) status <- "Sensitive"
        if (auc_val <= 0.25) status <- "Resistant"
      }
      
      gene_result <- rbind(gene_result,
                           data.frame(Drug=drug,
                                      Gene=gene,
                                      AUC=auc_val,
                                      ANOVA_P=aov_p,
                                      Label=status))
    }
    
    results <- rbind(results,gene_result)
  }
  
  write.csv(results,paste0(out_dir,"/result_auc_pvalue.csv"),row.names = F,quote = T)
  results_filter <- results[results$Label=="Sensitive",]
  write.csv(results_filter,paste0(out_dir,"/result_sensDPI.csv"),row.names = F,quote = T)
  
  p1 <- ggplot(results,aes(x=AUC,y=-log10(ANOVA_P)))+
    geom_point(color="gray",alpha=0.5)+
    geom_vline(xintercept = 0.75,color="gray",linetype="dashed",linewidth=1)+
    geom_hline(yintercept = -log10(0.05),color="gray",linetype="dashed",linewidth=1)+
    geom_point(data=results_filter,aes(x=AUC,y=-log10(ANOVA_P)),color="#a82424")+
    labs(x="AUC",y="ANOVA -log10(pvalue)")+
    theme_classic(base_size = 14)+
    theme(axis.text = element_text(color="black",size = 12))
  
  p2 <- ggplot(results_filter,aes(x=Gene))+
    geom_bar(stat = "count",width = 0.6,color="#a82424",fill="#e6adad")+
    #scale_y_continuous(limits = c(0,10),breaks = seq(0,10,2))+
    coord_flip()+
    labs(x="Protein",y="The number of drug")+
    theme_classic(base_size = 14)+
    theme(axis.text = element_text(color="black",size = 12))
  p <- plot_grid(p1, p2,ncol = 2,align = "h")#labels = LETTERS[1:4]
  ggsave(filename = paste0(out_dir,"/Plot.pdf"),p,width = 8,height = 4)
  
}

#sensPredict(node_file="node1.txt",cell="LUSC")
