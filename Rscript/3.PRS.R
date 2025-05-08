
library(oncoPredict)
library(reticulate)
library(dplyr)

getDRN <- function(expr_file){
  
  testData <- read.table(expr_file,header=T,sep = "\t",row.names = 1)
  GDSC1_Expr <- read.table("data/Signatrue/Pharmacogenomic/GDSC/data_exp.txt",header = T,sep = "\t",check.names = F)
  GDSC1_Res <- read.table("data/Signatrue/Pharmacogenomic/GDSC/data_act.txt",header = T,sep = "\t",check.names = F)
  
  out_dir <- "result/3.PRS_analysis/A.DRN"
  dir.create(out_dir,recursive = T)
  setwd(out_dir)
  
  # Predict patients' clinical responses based on module expressions.
  testExpr <- as.matrix(testData)
  calcPhenotype(trainingExprData = GDSC1_Expr,
                trainingPtype = GDSC1_Res,
                testExprData = testExpr,
                batchCorrect = 'eb',
                powerTransformPhenotype = TRUE,
                removeLowVaryingGenes = 0,
                minNumSamples = 10,
                printOutput = TRUE,
                removeLowVaringGenesFrom = 'rawData',
                cc = TRUE)
  
  moduleData <- data.frame(t(testData))
  moduleLabel <- apply(moduleData, 2, function(x) {ifelse(x > median(x), "high", "low")})
  
  drugData <- read.csv("./calcPhenotype_Output/DrugPredictions.csv",row.names = 1,check.names = F)
  dat <- merge(moduleLabel,drugData,by="row.names")
  rownames(dat) <- dat$Row.names
  dat <- dat[,-1]
  
  protein <- colnames(moduleData)
  drug <- colnames(drugData)
  
  proteins <- c()
  drugs <- c()
  pvalue <- c()
  for (pro in protein) {
    for (dg in drug) {
      proteins <- c(proteins,pro)
      drugs <- c(drugs,dg)
      
      temp = dat[,c(pro,dg)]
      dat.high = temp[temp[,1] == "high",2]
      dat.low = temp[temp[,1] != "high",2]
      if(all(is.na(dat.high))|all(is.na(dat.low))){
        pvalue <- c(pvalue,NA)
      }else{
        wilcox <- wilcox.test(dat.high,dat.low)
        pvalue <- c(pvalue,wilcox$p.value)
      }
    }
  }
  
  DPI_all <- data.frame("protein"=proteins,"drug"=drugs,"pvalue"=pvalue)
  DPI_all$label <- ifelse(DPI_all$pvalue < 0.05,"sig","non-sig")
  DRN <- subset(DPI_all,DPI_all$pvalue < 0.05)
  write.table(DRN,"DRN.txt",sep = "\t",quote = F,row.names = F)
  DRN.info <- data.frame(node = c(protein,drug),type = c(rep("protein",length(protein)),rep("drug",length(drug))))
  write.table(DRN.info,"DRN_info.txt",sep = "\t",quote = F,row.names = F)
  
  setwd("../..")
  
}

#--------------------------------------------------------------------------------------------------------
getAffinity <- function(smiles_file,seq_file,DRN_file,py_env,pretrained_model='MPNN_CNN_BindingDB_IC50'){
  
  out_dir <- "result/3.PRS_analysis/B.DRN_weight"
  dir.create(out_dir,recursive = T)
  
   #Drug smiles
  smiles <- read.table(smiles_file,header = T,sep = "\t",comment.char = "")
  colnames(smiles) <- c("drug","SMILES")
  row.names(smiles) <- smiles$drug
  
  #Protein sequences
  seqs <- read.table(seq_file,header = T,sep = "\t",comment.char = "")
  colnames(seqs) <- c("protein","Sequence")
  row.names(seqs) <- seqs$protein
  
  #Drug-Protein Interactions
  DPI <- read.table(DRN_file,header = T,sep = "\t")
  colnames(DPI) <- c("protein","drug")
  DPI <- DPI[DPI$protein %in% seqs$protein & DPI$drug %in% smiles$drug,]

  setwd(out_dir)
  
  # Integrate SMILES and Sequence to DPI
  DPI$Sequence <- apply(DPI, 1, function(p) seqs[p[1],"Sequence"])
  DPI$SMILES <- apply(DPI, 1, function(p) smiles[p[2],"SMILES"])
  DPI_info <- DPI[,c("drug","protein","Sequence","SMILES")]
  DPI_info <- DPI_info[DPI_info$SMILES !="",]
  colnames(DPI_info) <- c("drug_name","target_name","target_seq","drug_smiles")
  write.table(DPI_info,"./DPI_info.txt",sep = "\t",quote = F,row.names = F)
  
  #-----------
  # Activate the virtual environment for DeepPurpose
  use_condaenv(py_env)
  py_run_string("from DeepPurpose import DTI as models")
  
  # Import DPI information
  df = read.table("DPI_info.txt",header = T, sep="\t",comment.char = "")
  target_seq = as.list(df$target_seq)
  drug_smiles =  as.list(df$drug_smiles)
  target_name =  as.list(df$target_name)
  drug_name =  as.list(df$drug_name)
  # Import pre-trained model
  model = py$models$model_pretrained(model = pretrained_model)
  # Predict binding affinity
  my_predict = py$models$virtual_screening(drug_smiles, target_seq, model, drug_name, target_name,verbose=FALSE)
  
  setwd("../..")
}

getPS <- function(BA_file,PPIN_file,py_env="C:/Users/Bin/.conda/envs/DeepPurpose"){
  
  out_dir <- "result/3.PRS_analysis/C.PS"
  dir.create(out_dir,recursive = T)
  
  use_condaenv(py_env)
  py_run_string("import os")
  py_run_string("import numpy as np")
  py_run_string("import networkx as nx")
  py_run_string("import pandas as pd")
  py_run_string("import sys")
  
  enm_path <- "inst/python"
  py_run_string(paste0("sys.path.append(r'",enm_path,"')"))
  py_run_string("from enm.Enm import *")
  
  enm = py$Enm('PPIN')
  enm$read_network(PPIN_file, sep='\t')
  enm$gnm_analysis(normalized=FALSE)
  enm$cluster_matrix(enm$prs_mat)
  write.csv(enm$df,paste0(out_dir,"/pcc_df.csv"),row.names = T)
  write.table(enm$prs_mat_df,paste0(out_dir,"/prs_mat_df.txt"),quote = F)
  
  BA <- read.table(BA_file,skip = 3,sep = "|",comment.char ="+")
  BA <- BA[,2:5]
  colnames(BA) <- c("Rank", "Drug.Name", "Target.Name", "Binding.Score")
  BA$Drug.Name <- trimws(BA$Drug.Name)
  BA$Target.Name <- trimws(BA$Target.Name)
  
  Sens <- read.csv(paste0(out_dir,"/pcc_df.csv"),row.names = 1)
  colnames(Sens)[1] <- "Target.Name"
  
  ps <- merge(BA,Sens,by="Target.Name")
  ps <- ps[,c('Drug.Name', 'Target.Name', 'Binding.Score', 'sens')]
  ps$ps <- ps$Binding.Score*ps$sens
  ps$ps_norm <- (ps$ps-min(ps$ps))/(max(ps$ps)-min(ps$ps))
  ps <- ps[order(ps$ps_norm,decreasing = T),]
  write.csv(ps,paste0(out_dir,"/prs_dti_score.csv"),row.names = F)
  
  drug_PS <-  ps %>%
    group_by(Drug.Name) %>%
    summarize(Target.Name=paste(sort(Target.Name),collapse = ", "),
              PS = sum(ps_norm))
  drug_PS <- drug_PS %>% arrange(-PS)
  write.csv(drug_PS,paste0(out_dir,"/drug_score.csv"),row.names = F)
}

