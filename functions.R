##-----------------------------edgeR---------------------------##
volcanoplotR <- function(dge.obj, logfc, p.adj, colup = "#d84b47", coldown = "#82bf56", colns = "#6e6d6e") {
  ##This function adds a new column (T or F) according to the FDR and LFC of each gene in edgeR list of DEG
  ##dge.obj: List with DEG
  ##logfC: logFC threshold used for the differential expression test
  ##p.adj: p.adj or FDR threshold to obtain significant genes
  ##Updated 15-abr-2021 by Rodolfo Chavez
  #Rename the lfc and padj columns from DESeq2
  #Color update 2-ago-2021 by Rodolfo Chavez
  if("log2FoldChange" %in% colnames(dge.obj) & "padj" %in% colnames(dge.obj)) {
    dge.obj <- dge.obj %>% rename(logFC = log2FoldChange, FDR = padj)
  }
  #Add an extra column according to the logfc and p.adj thresholds
  dge.obj <- dge.obj %>%
    mutate(condition = ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR < p.adj), "Over-expressed",
                              ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR < p.adj), "Sub-expressed",
                                     ifelse((dge.obj$logFC > logfc) & (dge.obj$FDR > p.adj), "NS",
                                            ifelse((dge.obj$logFC < -logfc) & (dge.obj$FDR > p.adj), "NS",
                                                   ifelse((dge.obj$logFC < logfc) & (dge.obj$FDR > p.adj), "NS", "NS"))))))
  ##Eliminate NA values
  if(TRUE %in% unique(is.na(dge.obj$FDR))) {
    dge.obj <- dge.obj %>% drop_na()
  }
  #Plot the results
  volcano_plot <- ggplot(dge.obj)+
    geom_point(aes(x = logFC, y = -log10(FDR), color = condition))+
    scale_color_manual(name = "Condition",
                       labels = paste(c("NS", "Over-expressed", "Sub-expressed"), c(sum(dge.obj$condition == "NS"), sum(dge.obj$condition == "Over-expressed"), sum(dge.obj$condition == "Sub-expressed"))),
                       values = c(colns, colup, coldown))+
    geom_vline(aes(xintercept = logfc), linetype = "dashed")+
    geom_vline(aes(xintercept = -logfc), linetype = "dashed")+
    geom_hline(aes(yintercept = -log10(p.adj)), linetype = "dashed")+
    theme(plot.title = element_text(face = "bold", size = 18), 
          axis.title = element_text(size = 18),
          legend.title = element_text(face = "bold", size = 15),
          legend.text = element_text(size = 15), 
          legend.position = "bottom")
  return(volcano_plot)
}


