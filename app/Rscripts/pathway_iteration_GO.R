ComputeEnrichment_GO <- function(pathway, pathway_Sym, univ,  
                                 hits, non.hits, file.name, siRNA.Score, iteration)
{
  # using list format
  #Note, pathways seem to have some redudant terms due to evidence code, so value will be smaller
  hits <- intersect(hits, univ) 
  non.hits <- intersect(non.hits, univ)
  #pathway <- lapply(pathway, function(x) {x[x %in% c(hits, non.hits) ]} )
  
  #Subsets out empty pathways
  unique.pathways <- as.character(names(pathway))
  
  #Currently consistent with other method!
  results_Entrez = fgsea::fora(pathways = pathway, genes = hits, universe = univ, minSize = 10, maxSize = 500)
  
  #return(results_Entrez)
  
  #Repeats the above with symbols
  hits_Sym = Entrez2Sym(hits)
  univ_Sym = Entrez2Sym(univ)
  
  #Its just faster to run it twice somehow then convert entrez -> symbol lol
  results_GeneSyms = fgsea::fora(pathways = pathway_Sym, genes = hits_Sym, universe = univ_Sym,  
                                 minSize = 10, maxSize = 500) %>%
    dplyr::select(pathway, overlapGenes) %>%
    dplyr::rename(HitGeneNames = overlapGenes)
    
  results <- results_Entrez %>%
    full_join(results_GeneSyms, by = 'pathway') %>%
    dplyr::select(-padj) %>%
    dplyr::rename(Pathway = pathway, pVal = pval, HitGenes = overlap, Genes = size) %>%
    mutate(
      pValFDR = p.adjust(pVal,method = "BH"),
      pValBonferroni = p.adjust(pVal,method = "bonferroni")
    ) %>%
    rowwise() %>%
    mutate(
      HitGeneNames = paste(HitGeneNames, collapse = ', ')
    ) %>%
    ungroup() %>%
    dplyr::select(Pathway, pVal, pValFDR, pValBonferroni, Genes, HitGenes, HitGeneNames)
  
  # results <- results %>%
  #   rowwise() %>%
  #   mutate(
  #     HitGeneNames = paste0(humanGenes[humanGenes$EntrezID %in% overlapGenes,"GeneSymbol"], collapse = ", ")
  #   ) %>%
  #   dplyr::select(-overlapGenes) %>%
  #   ungroup()
    
    #TODO
    #figure out adjustment problem
    #see if can move hitgene names out of this section (or convert everything to symbols earlier?)
    #update so it can work for other species
  
  #This section writes out the results into a CSV for reading later
  results <- results[with(results, order(results$pValBonferroni,results$pValFDR,results$pVal)),]

  write.csv(results,file = paste0(file.name,".Enrichment_", iteration, ".csv"), row.names = FALSE)
  print('Determind sigs')
  
  sigPathways <- results$Pathway[which(results$pVal < 0.05)]
   
  sigPathwaysGenes <- unique(unlist(pathway[sigPathways]))
  
  if(length(sigPathwaysGenes) > 0){
    nonSigPathwaysGenes <- setdiff(unique(unlist(pathway)), sigPathwaysGenes)
  } else {
    nonSigPathwaysGenes <- unique(unlist(pathway))
  }
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

if(2 ==1) {
  
  temp <- as.list(org.Hs.egGO2ALLEGS)
  
  GO_Subset = as.list(GO.db::GOBPOFFSPRING)
  temp2 <- temp[names(temp) %in% names(GO_Subset)]
  
  temp3 <- select(org.Hs.eg.db, names(GO_Subset), "SYMBOL", "GOALL")
  
  hits <- sample(unique(unlist(temp)), 500, replace = FALSE)
  nonhits <- sample(unique(unlist(temp)), 1500, replace = FALSE)
  nonhits <- nonhits[!nonhits %in% hits]
  
  #Creates list of genes in GeneSymbol form 
  #General trend on how to create a version of this thats faster 
  #Translate hit genes in output to entrez at end
  hold <- lapply(temp2, function(x) humanGenes[humanGenes$EntrezID %in% x, 'GeneSymbol'])
  
  out_alt <- fgsea::fora(pathways = temp2, genes = hits, 
              universe = unique(unlist(temp$`GO:0008150`)), 
              minSize = 10, maxSize = 500 ) 
  
  system.time({out <- ComputeEnrichment_GO(pathway = temp2, pathway_Sym = hold, 
                                        univ = unique(unlist(temp$`GO:0008150`)),
                                        hits = hits, non.hits = nonhits) })
}


