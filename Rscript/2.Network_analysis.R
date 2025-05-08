
#expr_file  <- "I:/Bin/F/LSCC.zzy/bio_analysis/0.fromzzy/expr/1_na_replace_mean.xlsx"
#group_file <- "I:/Bin/F/LSCC.zzy/bio_analysis/0.raw_data/P/LUAD1_group.txt"
#ppi_file <- "physical_interactome.txt"
library(dplyr)
library(ggplot2)
library(ggrepel)
library(igraph)
library(openxlsx)
library(org.Hs.eg.db)
library(clusterProfiler)
library(AnnotationDbi)
library(R.utils)
library(tidyverse)
R.utils::setOption("clusterProfiler.download.method",'auto')
library(msigdbr)
library(survminer)
library(survival)
library(WGCNA)
library(glmnet)
library(logistf)
library(randomForest)
#library(caret)
library(pROC)

#---------------------------------------------------------------------------------------
getModule <- function(expr_file,group_file,node_cutoff=10){
  
  out_dir <- "result/2.Network_analysis/A.module_division"
  dir.create(out_dir,recursive = T)
  
  expr <- read.table(expr_file,header=T,sep = "\t",row.names = 1)
  group <- read.table(group_file,header=T,sep = "\t")
  colnames(group) <- c("Sample","Group")
  tSample = group[group$Group == "Tumor","Sample"]
  nSample = group[group$Group == "NAT","Sample"]
  
  pvalue = rep(NA,nrow(expr))
  tMean = rep(NA,nrow(expr))
  nMean = rep(NA,nrow(expr))
  for (i in 1:nrow(expr)) {
    tExpr = as.numeric(expr[i,tSample])
    nExpr = as.numeric(expr[i,nSample])
    
    Test = wilcox.test(tExpr,nExpr)
    pvalue[i] = Test$p.value
    tMean[i] = mean(tExpr)
    nMean[i] = mean(nExpr)
    
    if(i%%1000 == 0)cat(paste0(i,"/",nrow(expr),"\r"))
  }
  result = data.frame(protein=row.names(expr),tMean=tMean,nMean=nMean,p.value = pvalue)
  result$log2FC <- ifelse(result$tMean > 10000, log2(result$tMean/result$nMean), tMean-nMean)
  result$label <- rep("non-significant",nrow(result))
  result$label[result$p.value < 0.01 & result$log2FC > 1] = "up"
  result$label[result$p.value < 0.01 & result$log2FC < -1] = "down"
  #write.csv(result,paste0(out_dir,"/all_proteins.csv"),row.names = F,quote = F)
  
  DEPs <- subset(result,result$label != "non-significant")
  write.csv(DEs,paste0(out_dir,"/DEPs.csv"),row.names = F)
  
  # Map DEPs to physical protein-protein interaction network
  ppi_net <- read.table("data/Signatrue/Physical/physical_interactome.txt",header = T,sep = "\t")
  colnames(ppi_net) <- c("node1","node2")
  network <- ppi_net[ppi_net$node1%in%DEPs$protein & ppi_net$node2%in%DEPs$protein,]
  
  DEPs_expr <- t(expr[DEPs$protein,])
  DEPs_cor <- cor(DEPs_expr,method = "pearson")
  network$cor <- apply(network,1,function(x) DEPs_expr[x[1],x[2]])
  
  for (i in 1:nrow(network)) {
    network$cor[i] <- DEPs_cor[network[i,1],network[i,2]]
  }
  network_final <- network[abs(network$cor)>0.3,]
  write.table(network_final,paste0(out_dir,"/network.txt"),row.names = F,sep = "\t",quote = F)
  
  # Create the network
  edges <- network_final[,1:2]
  g_net <- graph_from_data_frame(edges, directed=FALSE)
  
  # Perform community detection using GN and LPA
  set.seed(123)
  GN <- cluster_edge_betweenness(g_net,weights=NULL)
  set.seed(123)
  LP <- cluster_label_prop(g_net,weights=NULL)
  
  # Associate each node with its corresponding module
  GN_label <- data.frame(node = get.vertex.attribute(g_net)[[1]],
                         module = GN$membership)
  LP_label <- data.frame(node = get.vertex.attribute(g_net)[[1]],
                         module = LP$membership)
  
  # Extract objects in different clusters
  GN_label_list <- split(GN_label$node,GN_label$module)
  LP_label_list <- split(LP_label$node,LP_label$module)
  
  # Remove outliers that belong to no cluster
  GN_label_list2<- GN_label_list[(lengths(GN_label_list) >1)]
  LP_label_list2<- LP_label_list[(lengths(LP_label_list) >1)]
  
  GN_LP <- as.data.frame(matrix(nrow=length(LP_label_list2),
                                ncol=length(GN_label_list2)))
  GN_LP_union <- GN_LP
  GN_LP_phyper <- GN_LP
  
  # Hypergeometric test between each cluster pairs
  for(i in 1:length(GN_label_list2)){
    for(j in 1:length(LP_label_list2)){
      GN_LP[j,i]=length(intersect(GN_label_list2[[i]],
                                  LP_label_list2[[j]]))
      GN_LP_union[j,i]=length(union(GN_label_list2[[i]],
                                    LP_label_list2[[j]]))
    }
  }
  for(i in 1:ncol(GN_LP_phyper)){
    for(j in 1:nrow(GN_LP_phyper)){
      GN_LP_phyper[j,i]=1-phyper(GN_LP[j,i],
                                 length(GN_label_list2[[i]]),GN_LP_union[j,i],
                                 length(LP_label_list2[[j]]))
    }
  }
  colnames(GN_LP_phyper) <- paste("GN",1:length(GN_label_list2),sep="_")
  rownames(GN_LP_phyper) <- paste("LP",1:length(LP_label_list2),sep="_")
  
  # Extract significantly correlated clusters
  sig <- data.frame(which(GN_LP_phyper < 0.05, arr.ind = TRUE))
  sig_list <- apply(sig,1,function(x) intersect(GN_label_list2[[x[2]]],LP_label_list2[[x[1]]]))
  
  # Integrate the nodes included in robustness modules
  modules <- data.frame(node = unlist(sig_list),
                        module = rep(1:length(sig_list),times=as.vector(sapply(sig_list, length))))
  row.names(modules) <- modules$node
  modules$module <- paste0("M",modules$module)
  
  # Get module label of each interaction.
  edges <- edges[,1:2]
  edges$module <- apply(edges, 1, function(x) ifelse((!x[1]%in%modules$node)|(!x[2]%in%modules$node),"M0",
                                                     ifelse(modules[x[1],2]==modules[x[2],2],modules[x[1],2],"M0")))
  # Extract interactions of robustness modules
  patternM <- edges[!edges$module=="M0",]
  
  # Write the modules and edges data frames to separate output files
  write.table(modules,paste0(out_dir,"/node_Module.txt"),sep="\t",quote=F,row.names=F,col.names=c("node","module"))
  write.table(patternM,paste0(out_dir,"/edge_Module.txt"),sep="\t",quote=F,row.names=F,col.names=c("node1","node2","module"))
  write.table(edges,paste0(out_dir,"/edges.txt"),sep="\t",quote=F,row.names=F,col.names=c("node1","node2","module"))
  message("Module identification succeeded!")
  
  # Select the modules
  node_count <- data.frame(table(modules$module))
  colnames(node_count) <- c("module","count")
  count1 <- node_count[node_count$count >= node_cutoff,]
  node_secl <- modules[modules$module %in% count1$module,]
  write.table(node_secl,paste0(out_dir,"/node_Module_select.txt"),sep = "\t",quote = FALSE,row.names = FALSE)
  message("Modules with no less than 10 nodes are selected!")
}

#---------------------------------------------------------------------------------------
modl_Annot <- function(node_file){
  out_dir <- "result/2.Network_analysis/A.module_annotation"
  dir.create(out_dir,recursive = T)
  
  node_modl <- read.table(node_file,header = T)
  modl_list <- unique(node_modl$module)
  
  for (i in 1:length(modl_list)){
    #i=1
    cat(paste0(i,"/",length(modl_list),"\r"))
    modl_file <- paste0(out_dir,"/",modl_list[i])
    dir.create(path = modl_file,recursive = T)
    perM <-  modl_list[i]
    
    #GO
    perM_gene <- node_modl[which(node_modl$module == perM),"node"]
    go <- enrichGO(perM_gene, OrgDb = org.Hs.eg.db, ont='ALL', 
                   pAdjustMethod = "none", pvalueCutoff = 0.05, 
                   qvalueCutoff = 1, keyType = 'SYMBOL')
    write.xlsx(go,paste0(modl_file,"/GO.xlsx"))
    
    #KEGG
    entrezid <- AnnotationDbi::select(org.Hs.eg.db, keys=perM_gene, columns = 'ENTREZID', keytype = 'SYMBOL')
    no_map <- sort(as.character(entrezid[is.na(entrezid$ENTREZID),'SYMBOL']))
    alias <- AnnotationDbi::select(org.Hs.eg.db, keys=no_map, columns = c('SYMBOL', 'ENTREZID'), keytype = 'ALIAS')
    no_map_id <-  alias[,-2]
    colnames(no_map_id) <- colnames(entrezid)
    geneID <- rbind(entrezid[!is.na(entrezid$ENTREZID),], no_map_id)
    
    perM_id <- geneID$ENTREZID
    ekk <- enrichKEGG(gene = perM_id, organism  = 'hsa', 
                      pAdjustMethod = "none", pvalueCutoff = 0.05, 
                      qvalueCutoff = 1,keyType = 'kegg') 
    kegg <- setReadable(ekk,OrgDb = "org.Hs.eg.db",keyType = "ENTREZID")
    write.xlsx(kegg,paste0(modl_file,"/KEGG.xlsx"))
  }
}

#---------------------------------------------------------------------------------------
Hallmark <- function(node_file,FC_file){
  out_dir <- "result/2.Network_analysis/B.hallmark_gsea"
  dir.create(out_dir,recursive = T)
  
  node_modl <- read.table(node_file,header = T)
  colnames(node_modl) <- c("node","module")
  
  node_FC <- read.table(node_file,header = T)
  colnames(node_modl) <- c("node","log2FC")
  
  node_df <- merge(node_FC,node_modl,by="node")
  node_df <- node_df[order(node_df$log2FC, decreasing = TRUE),]
  gene_ranks <- node_df$log2FC
  names(gene_ranks) <- node_df$node
  
  # Hallmark genesets 
  Hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
  hallmark_gene <- Hallmark %>% 
    dplyr::select(gs_name, gene_symbol)
  
  # For all nodes
  gsea_results <- GSEA(
    geneList = gene_ranks, 
    TERM2GENE = hallmark_gene, 
    minGSSize = 5,
    seed = T,
    pvalueCutoff = 0.05,
    pAdjustMethod ="none"
  )
  write.xlsx(gsea_results,paste0(out_dir,"/hallmark.xlsx"))
  
  # For nodes in each module
  modl_list <- unique(node_df$module)
  for (i in 1:length(modl_list)){
    #i=1
    cat(paste0(i,"/",length(modl_list),"\r"))
    modl_file <- paste0(out_dir,"/",modl_list[i])
    dir.create(path = modl_file,recursive = T)
    perM <-  modl_list[i]
    
    perM_gene <- gene_ranks[names(gene_ranks)%in%node_modl[node_modl$module==perM,1]]
    
    hallmarker <- GSEA(
      geneList = perM_gene, 
      TERM2GENE = hallmark_gene, 
      minGSSize = 5,
      pvalueCutoff = 0.05,
      pAdjustMethod ="none"
    )
    write.xlsx(hallmarker,paste0(modl_file,"/Hallmark.xlsx"))
    rm(list = c("perM_gene","hallmarker"))
  }
}

#---------------------------------------------------------------------------------------

Prognosis <- function(node_file,expr_file,surv_file,hit_time=1000){
  out_dir <- "result/2.Network_analysis/C.prognosis_association"
  dir.create(out_dir,recursive = T)
  
  exprData <- read.table(expr_file,header = T,row.names = 1,sep = "\t")
  survData <- read.table(surv_file,header = T,sep = "\t")
  colnames(survData) <- c("Sample","Time","Status")
    
  node_info <- read.table(node_file,header = T)
  colnames(node_info) <- c("node","module","regulate")
  
  modl_list <- unique(node_info$module)
  surv_info <- c()
  for (i in 1:length(modl_list)) {
    cat(paste0(i,"/",length(modl_list),"\r"))
    
    perM <- modl_list[i]
    modl_file <- paste0(out_dir,"/",perM) 
    dir.create(modl_file,recursive = T)
    
    perM_info <- node_info[node_info$module %in% perM,]
    perM_pro <- perM_info[,1]
    
    perM_expr <- data.frame(t(exprData[perM_pro,]))
    perM_expr$Sample <- row.names(perM_expr)
    perM_surv <- merge(survData,perM_expr,by="Sample")
    
    res.cut <- surv_cutpoint(perM_surv, 
                             time = "Time",
                             event = "Status",
                             variables = colnames(perM_surv)[-c(1:3)],
                             minprop = 0.1,
                             progressbar =F
    )
    res_label <- surv_categorize(res.cut)
    
    res_surv <- c()
    pros <- c()
    specifc_time <- hit_time 
    for (index in 3:ncol(res_label)) {
      pro <- colnames(res_label)[index]
      tempOS <-  res_label[,c(1,2,index)]
      colnames(tempOS)[3] <- "Group"
      tempOS <- na.omit(tempOS)
      if(any(grepl(paste0(as.character(c(1:9)),collapse = "|"),tempOS$Group))==T){
        tempOS$Group <- ifelse(tempOS$Group > median(as.numeric(tempOS$Group)),"high","low")
      }
      tempOS$Group <- factor(tempOS$Group,levels = c("high","low"))
      if (length(unique(tempOS$Group)) < 2 | length(unique(tempOS$Status)) < 2) {
        surv_pval <- NA
        surv_prob <- c(NA, NA)
      } else {
        fit <- survfit(Surv(Time, Status) ~ Group, data = tempOS)
        diff <- tryCatch({
          survdiff(Surv(Time, Status) ~ Group, data = tempOS)
        }, error = function(e) return(NULL))
        
        if (is.null(diff) || is.null(diff$chisq)) {
          surv_pval <- NA
        } else {
          surv_pval <- 1 - pchisq(diff$chisq, length(diff$n) - 1)
        }
        
        surv_prob <- tryCatch({
          summary(fit, times = specifc_time)$surv
        }, error = function(e) {
          c(NA, NA)
        })
        pros <- c(pros,pro)
        res_surv <- rbind(res_surv,c(surv_pval,specifc_time,surv_prob))
        if(surv_pval < 0.05 & any(is.na(surv_prob))){
          p <- ggsurvplot(data=tempOS, 
                          fit =fit ,
                          pval = TRUE,
                          pval.method = TRUE,
                          pval.size = 4,
                          pval.coord = c(3, 0.08),
                          pval.method.coord = c(3, 0.18),
                          palette = "simpsons",
                          legend.title = pro
          )
          ggsave(paste0(modl_file,"/",pro,".pdf"),p$plot,width = 4,height = 4)
        }
      }
    }
    result <- cbind(pros,data.frame(res_surv)) %>%
      setNames(c("node","surv_pvalue","specifc_time","surv_prob_high","surv_prob_low"))%>%
      arrange(surv_pvalue >= 0.05, node) %>%
      mutate(surv_type = if_else(surv_pvalue > 0.05, "non-sig",
                                 if_else(surv_prob_high < surv_prob_low, "unfavor", "favor")))
    write.xlsx(result,paste0(modl_file,"/survival.xlsx"),rowNames = FALSE)
    surv_info <- rbind(surv_info,result)
  }
  node_surv <- merge(node_info,surv_info,by="node")
  node_surv$effect <- "non-effect"
  node_surv$effect[(node_surv$regulate=="up"&node_surv$surv_type=="unfavor")|
                     (node_surv$regulate=="down"&node_surv$surv_type=="favor")] <- "negative"
  node_surv$effect[(node_surv$regulate=="up"&node_surv$surv_type=="favor")|
                     (node_surv$regulate=="down"&node_surv$surv_type=="unfavor")] <- "positive"
  write.xlsx(node_surv,paste0(out_dir,"/node_effect.xlsx"),rowNames = FALSE)
  
  count <- data.frame(table(node_surv$module,node_surv$effect))
  colnames(count) <- c("module","label","count")
  count$module <- paste0("M",count$module)
  count$label <- factor(count$label,levels = c("positive","non-effect","negative"))
  
  count_pct <- count %>% 
    group_by(module) %>% 
    mutate(Percentage=round(count/sum(count),2),
           Label =  paste(Percentage*100,"%",sep = "")) %>% 
    arrange(label)
  write.xlsx(count_pct,paste0(out_dir,"/module_effect.xlsx"),rowNames = FALSE)
}
#Prognosis(node_file="node_Module_info.txt",expr_file="expression.txt",surv_file="LUAD_survival.txt")

#---------------------------------------------------------------------------------------
Preserve <- function(Pdata,Tdata,Pmodule){
  out_dir <- "result/2.Network_analysis/D.module_preservation"
  dir.create(out_dir,recursive = T)
  
  datExprP1 <- read.table(Pdata,header = T,row.names = 1,sep = "\t")
  datExprT1 <- read.table(Tdata,header = T,row.names = 1,sep = "\t")
  
  commonProbes1 <- intersect(rownames(datExprP1), rownames(datExprT1)) 
  datExprP1 <- datExprP1[commonProbes1,] 
  datExprT1 <- datExprT1[commonProbes1,] 
  
  modulesP <- read.table(Pmodule,header = T,sep = "\t",row.names = 1)
  modulesP1 <- modulesP[,-1]
  
  multiExpr <- list(P1=list(data=t(datExprP1)),T1=list(data=t(datExprT1))) 
  multiColor <- list(P1 = modulesP1)
  mp <- modulePreservation(multiExpr, multiColor, 
                          referenceNetworks = 1, 
                          verbose = 3) 
  stats <- mp$preservation$Z$ref.P1$inColumnsAlsoPresentIn.T1
  stats_module <- stats[order(-stats[,2]),c(1:2)]
  
  rownames(modulesP) <- NULL 
  modulesP <- unique(modulesP)
  modulesP <- modulesP[modulesP$color!="grey",]
  rownames(modulesP) <- modulesP$color
  stats_module <- stats_module[row.names(stats_module)%in%modulesP$color,]
  rownames(stats_module) <- modulesP[rownames(stats_module),"module"]
  stats_module$module <- rownames(stats_module)
  write.xlsx(stats_module,paste0(out_dir,"/Zsummary.xlsx"), colNames = T, rowNames = T)
  
  g <- ggplot(data=stats_module,aes(x=moduleSize,y=Zsummary.pres,col=module))+
    geom_point(alpha=0.8, size=5) +
    theme_bw(base_size=15)+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.line = element_line(colour = "black"))+
    xlab("Module Size") + ylab("Zsummary") +
    theme(plot.title = element_text(size=15,hjust = 0.5))+
    theme(legend.position='none')+
    geom_hline(yintercept = c(2,10),lty=4,lwd=1,col=c("blue","red"))+
    geom_text_repel(aes(label=module),color="black",alpha = 0.8)
  ggsave(g,filename = paste0(out_dir,"/Preservation Zsummary.pdf"),height = 5,width = 5)
}
#Preserve(Pdata="1_na_replace_mean.txt",Tdata="1_na_replace_meanT.txt",Pmodule="protein_intersect.txt")


#---------------------------------------------------------------------------------------
Similarity <- function(node_file1, node_file2){
  out_dir <- "result/2.Network_analysis/E.module_similarity"
  dir.create(out_dir,recursive = T)
  
  node1 <- read.table(node_file1,header = T,sep = "\t")
  colnames(node1) <- c("node","module")
  node2 <- read.table(node_file2,header = T,sep = "\t")
  colnames(node2) <- c("node","module")
  
  modl1 <- unique(node1$module)
  modl2 <- unique(node2$module)
  
  M1 <- c()
  M2 <- c()
  Jaccard <- c()
  for(m1 in modl1){
    for (m2 in modl2) {
      N1 <- node1[node1$module==m1,1]
      N2 <- node2[node2$module==m2,1]
      jac <- length(intersect(N1,N2))/length(union(N1,N2))
      Jaccard <- c(Jaccard,jac)
      M1 <- c(M1,m1)
      M2 <- c(M2,m2)
    }
  }
  res_Jac <- data.frame(M1,M2,Jaccard)
  text <- res_Jac[Jaccard!=0,]
  p <- ggplot(res_Jac,aes(x=M2,y=M1,fill=Jaccard))+
    geom_tile(color="gray",size=0.1) +
    scale_fill_gradient(low = "white", high = "red",limit= c(0,1))+
    geom_text(data=text,aes(label=round(Jaccard,3)),col ="black",size = 3)+
    theme_minimal(base_size = 12)+
    theme(axis.text = element_text(colour = "black",size = 11),
          axis.text.x = element_text(angle = 90,hjust =1,vjust = 0.5),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank())
  ggsave(p,filename=paste0(out_dir,"/similarity.pdf"),width = 6,height = 6)
  
}
#Similarity(node_file1 = "F:/LSCC.zzy/bio_analysis/LUAD1/0.module/node_Module.txt",node_file2 = "F:/LSCC.zzy/bio_analysis/LSCC/0.module/node_Module.txt")

#---------------------------------------------------------------------------------------
Dynamics <- function(edge_file,py_env="C:/Users/Bin/.conda/envs/DeepPurpose"){
  out_dir <- "result/2.Network_analysis/F.module_dynamics"
  dir.create(out_dir,recursive = T)
  
  library(reticulate)
  library(dplyr)
  
  use_condaenv(py_env)
  py_run_string("import os")
  py_run_string("import numpy as np")
  py_run_string("import networkx as nx")
  py_run_string("import pandas as pd")
  py_run_string("import sys")
  
  enm_path <- "inst/python"
  py_run_string(paste0("sys.path.append(r'",enm_path,"')"))
  py_run_string("from enm.Enm import *")

  edge_modl <- read.table(edge_file,header = T,sep="\t")
  colnames(edge_modl) <- c("gene1","gene2","module")
  modl_list <- unique(edge_modl$module)
  
  sens_df <- c()
  for (modl in modl_list) {
    modl_file <- paste0(out_dir,"/",modl) 
    dir.create(modl_file,recursive = T)
    
    per_edge <- edge_modl[edge_modl$module == modl,1:2]
    write.table(per_edge,paste0(modl_file,"/edge.txt"),sep = "\t",quote = F,row.names = F)
    
    enm = py$Enm('PPIN')
    PPIN_file <- paste0(modl_file,"/edge.txt")#paste0(modl,"_edge.txt")
    enm$read_network(PPIN_file, sep='\t')
    enm$gnm_analysis(normalized=FALSE)
    enm$cluster_matrix(enm$prs_mat)
    write.csv(enm$df,paste0(modl_file,"/pcc_df.csv"),row.names = T)
    write.table(enm$prs_mat,paste0(modl_file,"/prs_df.txt"),quote = F)
    write.table(enm$prs_mat_df,paste0(modl_file,"/prs_mat_df.txt"),quote = F)
    
    temp_df <- enm$df
    temp_df$module <- modl
    sens_df <- rbind(sens_df,temp_df)
  }
  p <- ggplot(sens_df,aes(x=reorder(module, sens, FUN = median),y=sens))+
    geom_boxplot()+#outlier.shape = NA
    #coord_cartesian(ylim =  c(0, 1))+
    labs(x="Module",y="Sensitivity")+
    theme_bw(base_size = 12)+
    theme(axis.text = element_text(colour = "black",size = 11))
  ggsave(filename = paste0(out_dir,"/sensitivity_module.pdf"),width = 6,height = 4)
}
Dynamic(edge_file="I:/Bin/F/LSCC.zzy/bio_analysis/LSCC/0.module/edge_Module.txt")

#---------------------------------------------------------------------------------------
Independence <- function(node_file,expr_file1,expr_file2){
  out_dir <- "result/2.Network_analysis/G.module_similarity"
  dir.create(out_dir,recursive = T)
  
  node_data <- read.table(node_file,header = T,sep = "\t")
  colnames(node_data) <- "node"
  nodes <- node_data$node
  
  expr_data <- read.table(expr_file1,header = T,row.names = 1,sep = "\t")
  group <- ifelse(grepl("N$",colnames(expr_data))==T,"Normal","Tumor")
  group_data <- data.frame(Sample=colnames(expr_data),Group=group)
  expr_matrix <- t(expr_data[nodes,])
  
  # Divide the training set and testing set
  tumor_data <- expr_matrix[group_data$Group=="Tumor",]
  normal_data <- expr_matrix[group_data$Group=="Normal",]
  split_ratio <- 0.5
  set.seed(123) 
  tumor_split <- sample(nrow(tumor_data), size = round(nrow(tumor_data) * split_ratio))
  normal_split <- sample(nrow(normal_data), size = round(nrow(normal_data) * split_ratio))
  # Training set
  tumor_train <- tumor_data[tumor_split, ]
  normal_train <- normal_data[normal_split, ]
  # Testing set
  tumor_test <- tumor_data[-tumor_split, ]
  normal_test <- normal_data[-normal_split, ]
  # Merge
  train_data <- rbind(tumor_train, normal_train)
  test_data <- rbind(tumor_test, normal_test)
  train_group <- group_data[group_data$Sample %in% row.names(train_data),]
  test_group <- group_data[group_data$Sample %in% row.names(test_data),]
  train_data <- scale(train_data[train_group$Sample,])# 对表达矩阵进行标准化
  test_data <- scale(test_data[test_group$Sample,])# 对表达矩阵进行标准化
  
  # Feature selection on the training set
  response <- factor(ifelse(train_group$Group == "Tumor", 1, 0))
  
  # Regression
  set.seed(123)
  fit <- glmnet(train_data, response, alpha = 0, family = "binomial")
  cv_fit <- cv.glmnet(train_data, response, alpha = 8, family = "binomial")
  oneSElambda <- cv_fit$lambda.1se
  final_model <- glmnet(train_data, response, alpha = 0, family = "binomial", lambda = oneSElambda)
  
  # Extract genes with non-zero coefficients
  nonzero_coef <- coef(final_model)
  nonzero_genes <- rownames(nonzero_coef)[which(nonzero_coef != 0)][-1]
  write.table(nonzero_genes, paste0(out_dir, "/features.txt"), row.names = FALSE, col.names = FALSE, quote = FALSE)
  
  # Calculate Odds Ratio on the training set
  train_data_nonzero <- train_data[, nonzero_genes]
  OR_results <- data.frame(Protein = character(), OddsRatio = c(), CI_Lower = c(), CI_Upper = c(),p=c(), stringsAsFactors = FALSE)
  
  for (protein in colnames(train_data_nonzero)) {
    tryCatch({
      fit <- logistf(response ~ train_data_nonzero[, protein], plcontrol = logistpl.control(maxit = 500))
      coef <- exp(fit$coefficients[2])
      ci_lower <- exp(fit$ci.lower[2])
      ci_upper <- exp(fit$ci.upper[2])
      p <- fit$prob[2]
      OR_results <- rbind(OR_results, data.frame(Protein = protein, OddsRatio = coef, CI_Lower = ci_lower, CI_Upper = ci_upper,p=p, stringsAsFactors = FALSE))
    }, error = function(e) {})
  }
  
  rownames(OR_results) <- OR_results$Protein
  write.table(OR_results, paste0(out_dir, "/OR_results.txt"), row.names = FALSE, col.names = TRUE, quote = FALSE)
  
  # Calculate the Risk Score on the testing set
  log_OR <- OR_results$OddsRatio
  names(log_OR) <- OR_results$Protein
  test_data_nonzero <- test_data[, OR_results$Protein]
  
  z_scores <- data.frame(test_data_nonzero)

  risk_scores_manual <- numeric(nrow(z_scores))
  for (i in 1:nrow(z_scores)) {
    risk_score_i <- 0
    for (j in 1:length(log_OR)) {
      risk_score_i <- risk_score_i + z_scores[i, j] * log_OR[j]
    }
    risk_scores_manual[i] <- risk_score_i
  }
  risk_scores_manual <- as.data.frame(risk_scores_manual)
  colnames(risk_scores_manual) <- "score"
  risk_scores_manual$Sample <- row.names(z_scores)
  risk_scores_manual <- merge(risk_scores_manual,test_group,by="Sample")
  risk_scores_manual$group <- ifelse(risk_scores_manual$Group == "Tumor", "Tumor", "Normal")
  write.table(risk_scores_manual, paste0(out_dir, "/risk_scores_test.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  p1 <- ggplot(risk_scores_manual, aes(x = group, y = score, fill = group)) +
    geom_violin(alpha=0.2,width=0.9,position=position_dodge(width=0.8),size=0.75,trim = FALSE,color="white")+
    geom_boxplot(alpha = 0.7,width = 0.08) +
    scale_fill_manual(values = c("#E69F00","#009E73")) +
    labs(title = "Test cohort", x = "Group", y = "Risk Score") +
    theme_test(base_size = 12)+
    theme(axis.text = element_text(color="black",size=11),
          legend.position = "none")+
    stat_compare_means()
  #p1
  ggsave(filename=paste0(out_dir, "/risk_score_test.pdf"),p1,width = 3.2,height = 3.3)
  
  # ROC analysis
  roc_curve <- roc(risk_scores_manual$group, risk_scores_manual$score)
  auc_value <- auc(roc_curve)
  p2 <- ggroc(roc_curve,legacy.axes = TRUE,size = 1,color="red")+
    geom_abline(color = "grey")+
    annotate("text",x=0.6,y=0.5,label=paste0("AUC:",round(auc_value,3)),color="red")+
    theme_test(base_size = 12)+
    theme(axis.text = element_text(color = "black",size=11))
  #p2
  ggsave(filename=paste0(out_dir, "/ROC_test.pdf"),p2,width = 3,height = 3)

  ###### Independent cohort
  expr_data_indep <- read.table(expr_file2,header = T,row.names = 1,sep = "\t")
  group_indep <- ifelse(grepl("N$",colnames(expr_data_indep))==T,"Normal","Tumor")
  group_data_indep <- data.frame(Sample=colnames(expr_data_indep),Group=group_indep)
  
  comm_node <- intersect(nodes,row.names(expr_data_indep))
  expr_matrix_indep <- t(expr_data_indep[comm_node,])
  
  OR_results_comm <- OR_results[comm_node,]
  log_OR_indep <- OR_results_comm$OddsRatio
  names(log_OR_indep) <- OR_results_comm$Protein
  
  z_scores_indep <- data.frame(expr_matrix_indep)
  
  # Calculate the Risk Score
  risk_scores_indep <- numeric(nrow(z_scores_indep))
  for (i in 1:nrow(z_scores_indep)) {
    risk_score_i <- 0
    for (j in 1:length(log_OR_indep)) {
      risk_score_i <- risk_score_i + z_scores_indep[i, j] * log_OR_indep[j]
    }
    risk_scores_indep[i] <- risk_score_i
  }
  risk_scores_indep <- as.data.frame(risk_scores_indep)
  risk_scores_indep$Sample <- row.names(z_scores_indep)
  risk_scores_indep <- merge(risk_scores_indep,group_data_indep,by="Sample")
  risk_scores_indep$group <- ifelse(risk_scores_indep$Group == "Tumor", "Tumor", "Normal")
  write.table(risk_scores_indep, paste0(out_dir, "/risk_scores_indep.txt"), row.names = TRUE, col.names = TRUE, quote = FALSE)
  
  p3 <- ggplot(risk_scores_indep, aes(x = group, y = risk_scores_indep, fill = group)) +
    geom_violin(alpha=0.2,width=0.9,position=position_dodge(width=0.8),size=0.75,trim = FALSE,color="white")+
    geom_boxplot(outlier.shape = NA, alpha = 0.7,width = 0.08) +
    scale_fill_manual(values = c("#E69F00","#009E73")) +
    labs(title = "Independent cohort", x = "Group", y = "Risk Score") +
    theme_test(base_size = 12)+
    theme(axis.text = element_text(color="black",size=11),
          legend.position = "none")+
    stat_compare_means()
  ggsave(filename=paste0(out_dir, "/risk_score_indep.pdf"),p3,width = 3.2,height = 3.3)
  
  # ROC analysis
  roc_curve_indep <- roc(risk_scores_indep$group, risk_scores_indep$risk_scores_indep)
  auc_value_indep <- auc(roc_curve_indep)
  p4 <- ggroc(roc_curve_indep,legacy.axes = TRUE,size = 1,color="#005cff")+
    geom_abline(color = "grey")+
    annotate("text",x=0.6,y=0.5,label=paste0("AUC:",round(auc_value_indep,3)),color="#005cff")+
    theme_test(base_size = 12)+
    theme(axis.text = element_text(color = "black",size=11))
  ggsave(filename=paste0(out_dir, "/ROC_indep.pdf"),p4,width = 3,height = 3)
}

#Independence(node_file="node.txt",expr_file1="1_na_replace_mean.txt",expr_file2="1_na_replace_meanT.txt")


