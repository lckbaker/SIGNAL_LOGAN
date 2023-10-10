ComputeEnrichment_KEGG <- function(pathway, hits, non.hits, file.name, siRNA.Score, iteration)
{
  hits <- intersect(hits,pathway$EntrezID)
  non.hits <- intersect(non.hits,pathway$EntrezID)
  pathway <- pathway[pathway$EntrezID %in% union(hits,non.hits),]
  unique.pathways <- unique(as.character(pathway$PathwayName))
  

  #Figure out what below does
  p.val <- pathway.genes.number <- hit.genes <- hit.gene.names <- rep(NA,length(unique.pathways))
  
  #This section calculates the p-value for each of the different pathways
  for(i in 1:length(unique.pathways))
  {
    #Subsets to just genes in the pathway
    pathway.genes <- pathway$EntrezID[which(as.character(pathway$PathwayName) == unique.pathways[i])]

    #Worst case, implement fora package to speed this up
    #Creates a 2x2 matrix that contain on/off pathways hits/not-hits
    contingency <- matrix(NA,nrow = 2, ncol = 2)
    contingency[1,1] <- length(intersect(pathway.genes,hits)) # pathway.genes.hits
    contingency[1,2] <- length(intersect(pathway.genes,non.hits)) # pathway.genes.non.hits
    contingency[2,1] <- length(setdiff(hits,pathway.genes)) # non.pathway.hits 
    contingency[2,2] <- length(setdiff(non.hits,pathway.genes)) # non.pathway.non.hits
    
    #This is apparently the same as a hypergeom test
    p.val[i] <- fisher.test(contingency,alternative = "greater")$p.value
    pathway.genes.number[i] <- contingency[1,1] + contingency[1,2]
    hit.genes[i] <- contingency[1,1]
    hit.gene.names[i] <- paste(unique(pathway$GeneSymbol[match(intersect(pathway.genes,hits),pathway$EntrezID)]),collapse = ", ")
    
    #Unassigned?
    #length(unique.pathways)
    p.val.FDR <- p.adjust(p.val,method = "BH")
    p.val.FWER <- p.adjust(p.val,method = "bonferroni")
   # message('unique')
    
    results <- data.frame(Pathway = unique.pathways, 
                          pVal = round(p.val, digits = 3), 
                          pValFDR = round(p.val.FDR, digits = 3), 
                          pValBonferroni = round(p.val.FWER, digits = 3), 
                          Genes = pathway.genes.number, 
                          HitGenes = hit.genes, 
    HitGeneNames = hit.gene.names)
  }

  
  #This section writes out the results into a CSV for reading later
  results <- results[with(results, order(results$pValBonferroni,results$pValFDR,results$pVal)),]
  print('prewrite')
  write.csv(results,file = paste0(file.name,".Enrichment_", iteration, ".csv"), row.names = FALSE)
  print('postwrite')
  print(head(results)) 
  sigPathways <- results$Pathway[which(results$pVal < 0.05)]
   sigPathwaysGenes <- unique(pathway$EntrezID[pathway$PathwayName %in% sigPathways])
   if(length(sigPathwaysGenes) > 0){
     nonSigPathwaysGenes <- setdiff(pathway$EntrezID, sigPathwaysGenes)
   } else {
     nonSigPathwaysGenes <- pathway$EntrezID 
   }
   
   print('Sig genes')
   print(sigPathways)
   print(sigPathwaysGenes)


  tempPathwayGenes <- matrix("Missing",nrow(siRNA.Score))
  message("sigPathways")
  message(sigPathways, "\n")
  message(file.name)
  tempPathwayGenes[siRNA.Score$EntrezID %in% sigPathwaysGenes] <- "Yes"
  tempPathwayGenes[siRNA.Score$EntrezID %in% nonSigPathwaysGenes] <- "No"
  if(length(grep("KEGG", file.name))>0){
    siRNA.Score$KEGG <- tempPathwayGenes
  } else if(length(grep("REACTOME", file.name))>0){
    siRNA.Score$REACTOME <- tempPathwayGenes
  } else if(length(grep("GO", file.name))>0){
    #siRNA.Score$GO[intersect(which(siRNA.Score$GO == "Missing"),which(tempPathwayGenes == "No"))] <- "No"
    #siRNA.Score$GO[siRNA.Score$EntrezID %in% sigPathwaysGenes] <- "Yes"
    
    siRNA.Score$KEGG <- tempPathwayGenes

  }
  return(siRNA.Score)
}

if(2 == 1) {
  #Not run on source
  temp <- as.list(org.Hs.egGO2ALLEGS)
  
  GO_Subset = as.list(GOBPOFFSPRING)
  temp2 <- temp[names(temp) %in% names(GO_Subset)]
  
  hits <- sample(unique(unlist(temp)), 500, replace = FALSE)
  nonhits <- sample(unique(unlist(temp)), 1500, replace = FALSE)
  nonhits <- nonhits[!nonhits %in% hits]
  
  out <- ComputeEnrichment_KEGG(pathway = temp, hits = hits, non.hits = nonhits)
  
  
  gene_univ = unique(unlist(temp$`GO:0008150`))
  
  #Creates list of genes in GeneSymbol form 
  #General trend on how to create a version of this thats faster 
    #Translate hit genes in output to entrez at end
  hold <- lapply(temp2, function(x) humanGenes[humanGenes$EntrezID %in% x, 'GeneSymbol'])
  gene_univ = unique(unlist(hold$`GO:0008150`))
  
  hits_GS <- humanGenes[humanGenes$EntrezID %in% hits, 'GeneSymbol']
  
  out_alt <- fgsea::fora(pathways = hold, genes = c(hits_GS) , 
                         universe = gene_univ, minSize = 10, maxSize = 500 ) %>%
    #dplyr::select(pathway, pval, overlap, size) %>%
    dplyr::rename(Pathway = pathway, pVal = pval, HitGenes = overlap, Genes = size) %>%
    mutate(
      pValFDR = p.adjust(pVal,method = "BH"),
      pValBonferroni = p.adjust(pVal,method = "bonferroni")
    ) %>%
    filter(HitGenes > 1)
  system.time({
    out_alt2 <- fgsea::fora(pathways = temp2, genes = hits, universe = gene_univ, minSize = 10, maxSize = 500 ) %>%
      #dplyr::select(pathway, pval, overlap, size) %>%
      dplyr::rename(Pathway = pathway, pVal = pval, HitGenes = overlap, Genes = size) %>%
      mutate(
        pValFDR = p.adjust(pVal,method = "BH"),
        pValBonferroni = p.adjust(pVal,method = "bonferroni")
      ) %>%
      filter(HitGenes > 1) %>%
      rowwise() %>%
      mutate(
        HitGeneNames = paste0(humanGenes[humanGenes$EntrezID %in% overlapGenes,"GeneSymbol"], collapse = ", ")
      )
  })
  
  system.time({
    out2 <- clusterProfiler::enrichGO(gene = hits,
                                      OrgDb = 'org.Hs.eg.db',
                                      ont = 'BP', 
                                      readable = TRUE,
                                      minGSSize = 10, maxGSSize = 500,
                                      pvalueCutoff = 1, qvalueCutoff = 1)
    out2 <- as.data.frame(out2)
  })
  
  out3 <- out2[1:5,] %>%
    rowwise() %>%
    mutate(
      TEST = paste0(Entrez2Sym(strsplit(geneID, "/")[[1]]), collapse = "," )
    )
  
  #TODO
  #need to fix universe to just be the whole dataset, will be necessary 
  #speed up id -> name or move it
  
  #hit.gene.names[i] <- paste(unique(pathway$GeneSymbol[match(intersect(pathway.genes,hits),pathway$EntrezID)])
  #,collapse = ", ")
  out2 <- clusterProfiler::enrichGO(gene = hits,
                                    OrgDb = 'org.Hs.eg.db',
                                    ont = 'BP', 
                                    readable = TRUE,
                                    pvalueCutoff = 1, qvalueCutoff = 1)
  out2 <- as.data.frame(out2)
  
  out3 <- clusterProfiler::enrichGO(gene = hits,
                                    OrgDb = 'org.Hs.eg.db',
                                    ont = 'BP', 
                                    readable = FALSE,
                                    pvalueCutoff = 1, qvalueCutoff = 1)
  out3 <- as.data.frame(out3)
  
  
  system.time({out <- ComputeEnrichment(pathway = temp, 
                                        hits = hits, non.hits = nonhits) })
  
  
  
  out2 <- clusterProfiler::enrichGO(gene = hits,
                                    OrgDb = 'org.Hs.eg.db',
                                    ont = 'BP', 
                                    readable = TRUE,
                                    minGSSize = 10, maxGSSize = 500,
                                    pvalueCutoff = 2, qvalueCutoff = 2)
  out2 <- as.data.frame(out2)
  
  
  rownames(out2)[!(rownames(out2) %in% out_alt$pathway)] %>% length()
  
  out_alt$pathway[!(out_alt$pathway) %in% rownames(out2)] %>% length()
  
  
  out_combine = out2[,c('ID', 'Description', 'pvalue')] %>%
    dplyr::rename(pathway = ID, pEnrGO = pvalue) %>%
    dplyr::full_join(out_alt[,c('pathway', 'pval')], by = c("pathway"))
}
