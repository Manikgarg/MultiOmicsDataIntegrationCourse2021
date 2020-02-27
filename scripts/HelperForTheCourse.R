calculateMode = function(x, na.rm = FALSE) {
  if(na.rm){
    x = x[!is.na(x)]
  }
  
  ux <- unique(x)
  return(ux[which.max(tabulate(match(x, ux)))])
}

# Modified from similar function in MOFA
visualizeWeightsHeatmap <- function(weights, view, threshold = 0){
  W = weights[[view]]
  # apply thresholding of loadings
  #W <- W[!apply(W,1,function(r) all(abs(r)<threshold)),]
  #W <- W[,!apply(W,2,function(r) all(abs(r)<threshold))]
  
  # Define breaks and center colorscale to white at zero
  minW <- min(W)
  maxW <- max(W)
  paletteLength <- 100
  colors <- c("black", "blue", "white", "orange","red")
  color <- colorRampPalette(colors)(paletteLength)
  breaks <- c(seq(minW, 0, length.out=ceiling(paletteLength/2) + 1), 
              seq(maxW/paletteLength,maxW, length.out=floor(paletteLength/2)))
  
  # Plot heatmap
  show_colnames <- ifelse(nrow(W)<80, TRUE, FALSE)
  cellwidth <- ifelse(show_colnames == TRUE, 6, NA)
  pheatmap(t(W), color=color, breaks=breaks, show_colnames = show_colnames, 
           fontsize_col = 6, 
           cellwidth = cellwidth)
}

# Modified from similar function in MOFA
plotTopWeightsUsingSeparateWeightsAndFactors <- function(weights, view, factor, nfeatures = 10, abs = TRUE, scale = TRUE, sign = "both") {
  
  # Sanity checks
  # if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  # stopifnot(view %in% viewNames(object))
  # # if(!is.null(manual_features)) { stopifnot(is(manual_features, "list"))
  # # stopifnot(all(Reduce(intersect,manual_features) %in% featureNames(object)[[view]]))  }
  
  # Collect expectations  
  W <- data.frame("value" = weights[[view]][, factor], "feature" = rownames(weights[[view]]))
  #W = data.frame("value" = rbind(W1[, 1], "feature" = rep(rownames(W1), 2), 
  #               "factor" = c(rep("LF1", nrow(W1)), rep("LF2", nrow(W1))))
  
  # Scale values by loading with highest (absolute) value
  if(scale) W$value <- W$value/max(abs(W$value))
  
  # store sign
  W <- W[W$value!=0,]
  W$sign <- ifelse(W$value>0, "+", "-")
  
  # select subset of only positive or negative loadings
  if (sign=="positive") { W <- W[W$value>0,] } else if (sign=="negative") { W <- W[W$value<0,] }
  
  # Absolute value
  if (abs) W$value <- abs(W$value)
  
  
  # Extract relevant features
  W <- W[with(W, order(-abs(value))), ]
  if (nfeatures>0) features <- head(W$feature,nfeatures) # Extract top hits
  # if (!is.null(manual_features)) 
  #    features <- W$feature[W$feature %in% manual_features] # Extract manual hits
  W <- W[W$feature %in% features,]
  
  # Sort according to loadings
  W <- W[with(W, order(-value, decreasing = TRUE)), ]
  W$feature <- factor(W$feature, levels=W$feature)
  W$start <- 0
  
  p <- ggplot(W, aes_string(x="feature", y="value")) +
    geom_point(size=2) +
    geom_segment(aes_string(yend="start", xend="feature"), size=0.75) +
    scale_colour_gradient(low="grey", high="black") +
    # scale_colour_manual(values=c("#F8766D","#00BFC4")) +
    # guides(colour = guide_legend(title.position="top", title.hjust = 0.5)) +
    coord_flip() +
    theme(
      axis.title.x = element_text(size=rel(1.5), color='black'),
      axis.title.y = element_blank(),
      axis.text.y = element_text(size=rel(1.2), hjust=1, color='black'),
      axis.text.x = element_text(size=rel(1.5), color='black'),
      axis.ticks.y = element_blank(),
      axis.ticks.x = element_line(),
      legend.position='top',
      # legend.title=element_text(size=rel(1.5), color="black"),
      legend.title=element_blank(),
      legend.text=element_text(size=rel(1.3), color="black"),
      legend.key=element_rect(fill='transparent'),
      panel.background = element_blank(),
      aspect.ratio = .7
    )
  
  if (sign=="negative") {p <- p + scale_x_discrete(position = "top")}
  
  # If absolute values are used, add the corresponding signs to the plot
  if (abs) {p <- p +  ylim(0,max(W$value)+0.1)+ geom_text(label=W$sign,y=max(W$value)+0.1, size=10)}
  
  if(abs & scale){ p <-  p + ylab(paste("Absolute loading on factor", factor))} 
  else if(abs & !scale) {p <- p + ylab(paste("Absolute loading on factor", factor))}
  else if(!abs & scale) {p <- p + ylab(paste("Loading on factor", factor))}
  else {p <- p + ylab(paste("Loading on factor", factor))}
  
  return(p)
  }

# Modified from Dietrich et al. (2018) BloodCancerMultiOmics2017 vignette  
extractMetadata <- function(lpdCLL, patmeta, drugs){
  lpdCLL <- lpdAll[ , lpdAll$Diagnosis=="CLL"   ]
  # data rearrangements
  survT = patmeta[colnames(lpdCLL),]
  survT[which(survT[,"IGHV"]=="U") ,"IGHV"] = 0
  survT[which(survT[,"IGHV"]=="M") ,"IGHV"] = 1
  survT$IGHV = as.numeric(survT$IGHV)
  
  colnames(survT) = gsub("Age4Main", "age", colnames(survT))
  
  survT$ibr45  <- 1-Biobase::exprs(lpdCLL)[ "D_002_4:5", rownames(survT)  ]  
  survT$ide45  <- 1-Biobase::exprs(lpdCLL)[ "D_003_4:5", rownames(survT)  ] 
  survT$prt45  <- 1-Biobase::exprs(lpdCLL)[ "D_166_4:5", rownames(survT)  ] 
  survT$selu45 <- 1-Biobase::exprs(lpdCLL)[ "D_012_4:5", rownames(survT)  ]
  survT$ever45 <- 1-Biobase::exprs(lpdCLL)[ "D_063_4:5", rownames(survT)  ]
  
  survT$nut15 <- 1-Biobase::exprs(lpdCLL)[ "D_010_1:5", rownames(survT)  ] 
  survT$dox15 <- 1-Biobase::exprs(lpdCLL)[ "D_159_1:5", rownames(survT)  ] 
  survT$flu15 <- 1-Biobase::exprs(lpdCLL)[ "D_006_1:5", rownames(survT)  ] 
  
  survT$SF3B1      <- Biobase::exprs(lpdCLL)[ "SF3B1",      rownames(survT)  ]
  survT$NOTCH1     <- Biobase::exprs(lpdCLL)[ "NOTCH1",     rownames(survT)  ]
  survT$BRAF       <- Biobase::exprs(lpdCLL)[ "BRAF",       rownames(survT)  ]
  survT$TP53       <- Biobase::exprs(lpdCLL)[ "TP53",       rownames(survT)  ]
  survT$del17p13   <- Biobase::exprs(lpdCLL)[ "del17p13",   rownames(survT)  ]
  survT$del11q22.3 <- Biobase::exprs(lpdCLL)[ "del11q22.3", rownames(survT)  ]
  survT$trisomy12 <-  Biobase::exprs(lpdCLL)[ "trisomy12", rownames(survT)  ]
  survT$IGHV_cont <- patmeta[ rownames(survT) ,"IGHV Uppsala % SHM"]
  
  # competinting risk endpoint fpr 
  survT$compE <- ifelse(survT$treatedAfter == TRUE, 1, 0)
  survT$compE <- ifelse(survT$treatedAfter == FALSE & survT$died==TRUE,
                        2, survT$compE )
  survT$T7  <- ifelse(survT$compE == 1, survT$T5, survT$T6 )
  
  return(survT)
}

forest <- function(Time, endpoint, factors, survT) {  
  
  res <- lapply(factors, function(g) {
    resultPerFactor <- apply(g, 2, function(x){
      surv <- coxph(Surv(survT[,Time], survT[,endpoint]) ~ x)  
      sumsu <- summary(surv) 
      return(c(p      = sumsu[["coefficients"]][, "Pr(>|z|)"], 
               HR = sumsu[["coefficients"]][, "exp(coef)"], 
               lower  = sumsu[["conf.int"]][, "lower .95"], 
               higher = sumsu[["conf.int"]][, "upper .95"]))
    })
    resultPerFactor <- t(rbind("padj" = p.adjust(resultPerFactor["p", ], method = "bonferroni"), resultPerFactor))
    return(resultPerFactor)
  })
  
  return(res)
}

# Modified from Dietrich et al. (2018) BloodCancerMultiOmics2017 vignette
plotSurvivalAnalysisResults <- function(factorsList, survT, title) {  
  ttt = forest(Time="T5", endpoint="treatedAfter", factors = factorsList, survT)
  os = forest(Time="T6", endpoint="died", factors = factorsList, survT)
  
  allSurvivalResults = data.frame()
  
  for (methodName in names(factorsList)){
    survivalResults = data.frame("padj" = ttt[[methodName]][, "padj"], 
                                 "method" = methodName, 
                                 "survivalType"= "TTT",
                                 "factor" = rownames(ttt[[methodName]]))
    allSurvivalResults = rbind(allSurvivalResults, survivalResults)
    
    survivalResults = data.frame("padj" = os[[methodName]][, "padj"], 
                                 "method" = methodName, 
                                 "survivalType"= "OS",
                                 "factor" = rownames(os[[methodName]]))
    allSurvivalResults = rbind(allSurvivalResults, survivalResults)
  }
  
  allSurvivalResults$significant = ifelse((allSurvivalResults$padj < 0.05), "Yes", "No")
  
  p1 <- ggplot(data = allSurvivalResults, aes(x = method, y = -log10(padj), colour = significant)) +
    geom_point() +
    scale_color_manual(values = c("black", "red"))+
    geom_hline(yintercept = -log10(0.05), linetype = 2, size = 0.2) +
    xlab("Method") + 
    ylab("-log (Bonferroni-corrected survival P-value)")+
    theme(legend.position="none") +
    facet_grid(.~ survivalType)
  plot(p1)
  
  for (methodName in names(factorsList)) {
    
    n <- c( padj = NA, p=NA, HR=NA, lower=NA, higher=NA )
    
    od <- order( c(seq(nrow(ttt[[methodName]])), seq(nrow(os[[methodName]])) ))
    
    s <- data.frame( rbind(ttt[[methodName]], os[[methodName]])[od,  ] )
    s$Name <- rownames(s)
    s$x <- 1:nrow(s)
    s$col <- rep(c("black", "darkgreen"), nrow(ttt[[methodName]]) )
    s$Endpoint <- factor( c(rep("TTT", nrow(ttt[[methodName]]) ),
                            rep("OS", nrow(os[[methodName]]) ) )[od] )
    s$features <- "";  s[ which(s$Endpoint=="OS"),"features"] <- rownames(ttt[[methodName]])
    
    p2 <- ggplot(data=s ,aes(x=x, y=HR, ymin=lower, ymax=higher,
                            colour=Endpoint)) +  geom_pointrange() + 
      ggtitle(methodName)+
      theme(legend.position="right", legend.text = element_text(size = 20) ) +
      scale_x_discrete(limits=s$x, labels=s$features ) +
      expand_limits(y=c(0,10)) +
      scale_y_log10(breaks=c(0.01,0.1,0.5,1,2,5,10),
                    labels=c(0.01,0.1,0.5,1,2,5,10)) +
      theme(
        panel.grid.minor = element_blank(),
        axis.title.x = element_text(size=16),
        axis.text.x = element_text(size=16, colour="black"),
        axis.title.y  = element_blank(),
        axis.text.y = element_text(size=12, colour="black"),
        legend.key = element_rect(fill = "white"),
        legend.background = element_rect(fill = "white"),
        legend.title = element_blank(),
        panel.background = element_rect(fill = "white", color="black"),
        panel.grid.major = element_blank(),
        axis.ticks.y = element_blank() 
      ) +
      coord_flip() + 
      scale_color_manual(values=c("OS"="darkgreen", "TTT"="black"),
                         labels=c("OS", "TTT")) +
      geom_hline(aes(yintercept=1), colour="black", size=1.5,
                 linetype="dashed", alpha=0.3)  +
      annotate("text", x = 1:nrow(s)+0.5, y = s$HR+0.003,
               label = ifelse( s$padj<0.001, paste0("padj<","0.001"), 
                               paste0("padj=", round(s$padj,3) ) ), colour=s$col)
    
    plot(p2)
    
  }
}


# Modified from Cantini et al. (2020) momix pipeline
## Perform biological annotation-based comparison 
## INPUTS:
# factorizations = already computed factirizations
# path.database = path to a GMT annotation file
# pval.thr = p-value threshold (default to 0.05)
## OUPUTS: a list containing output values
# selectivity = Selectivity (fraction of significant annotations per all significant factors)
# nonZeroFacs = Number of unique factors that were significant at least once
# total_pathways = Number of clinical annotations with at least one significant factor component
biologicalComparisonModifiedWithMofaCode <- function(metagenesList, factorsList, data, pathToFseaScript, database, pval.thr=0.01){
  
  # Load annotation database
  #pathways <- gmtPathways(pathToDatabase)
  
  # Containers to report results
  report_number <- numeric(0)
  report_nnzero <- numeric(0)
  report_select <- numeric(0)
  
  # For each factorization method
  for(methodName in names(metagenesList)){
    
    # Extract metagenes found by factorization method
    metagenes<-metagenesList[[methodName]]
    factors <- factorsList[[methodName]]
    # Number of factors
    num.factors <- ncol(factors)
    # # Rename columns
    # #colnames(metagenes)<-1:num.factors
    # #Rename rows to gene names
    # genSymbols <- read.table(pathToGeneSymbols, header = TRUE, sep = "\t", quote = "")
    # metagenes = metagenes[rownames(metagenes) %in% rownames(genSymbols), ]
    # rownames(metagenes) = genSymbols$hgnc_symbol[rownames(metagenes) %in% rownames(genSymbols)]
    # # Remove duplicated gene names that could confuse fgsea
    # duplicated_names <- unique(rownames(metagenes)[duplicated(rownames(metagenes))])
    # metagenes <- metagenes[!(rownames(metagenes) %in% duplicated_names), ]
    
    # Variables
    min_pval <- numeric(0)
    path <- numeric(0)
    n <- 0
    
    source(pathToFseaScript)
    gsea <- runEnrichmentAnalysisUsingSeparateWeightsAndFactors(
      metagenes, factors, data,
      view = "mRNA",
      feature.sets = database,
      alpha = pval.thr
    )
    gsea_padj = data.frame(gsea[["pval.adj"]])
    # Calculate biological annotation enrichment.
    # For each factor,
    for(j in 1:num.factors){
      # Assign gene names
      #nrnk <- setNames(as.matrix(metagenes[,j]), rownames(metagenes))
      # Compute fgsea
      # fgseaRes <- fgsea(pathways, rnk, minSize=15, maxSize=500, nperm=1000)
      # If at least one pathway is significant
      if(sum(gsea_padj[, j] < pval.thr)!=0){
        # Count this factor
        n <- n+1
        # Keep min adjusted p-value
        min_pval <- rbind(min_pval, min(gsea_padj[, j]))
        # Keep names of significant pathways
        path <- c(path, rownames(gsea_padj)[which(gsea_padj[, j]<pval.thr)])
      } else {
        min_pval <- rbind(min_pval, NA)
      }
    }
    
    # Report number of unique significant pathways  
    if(length(path)==0){
      report_number <- rbind(report_number, 0)
    }else{
      report_number <- rbind(report_number, length(unique(path)))
    }
    
    # Report selectivity 
    if(length(unique(path))==0){
      report_select <- rbind(report_select, NA)
    }else{
      al<-length(unique(path))/length(path)
      fl<-length(which(!is.na(min_pval)))/length(path)
      report_select <- rbind(report_select, (al+fl)/2)
    }
    # Report number of factors associated with at least one significant pathway
    report_nnzero<-rbind(report_nnzero, n)    
    
  }
  
  out <- data.frame(selectivity=report_select, nonZeroFacs=report_nnzero, total_pathways=report_number)
  rownames(out) <- names(metagenesList)
  #    print(out)
  return(out)
}

# Modified from Cantini et al. (2020) momix pipeline
calculateSelectivityScoreForClinicalAnnotations <- function(factorsList, clinicalData, clinicalFeaturesOfInterest, pval.thr){
  # Empty containers
  line <- numeric(0)
  line2 <- numeric(0)
  line3 <- numeric(0)
  
  # For each factorization
  for(methodName in names(factorsList)){
    factors <- factorsList[[methodName]]
    # Stor all p-values
    pvalues <- numeric(0)
    # Store number of significant annotations
    clin_erich <- 0 
    
    # Test significance association with clinical annotations
    for(j in clinicalFeaturesOfInterest){
      
      # Perform the analysis if there is more than one possible value in current column
      table_values <- table(clinicalData[,j])
      if(sum(table_values>0)>1){
        pvalues_col <- apply(factors, MARGIN=2, 
                             function(x) wilcox.test(x~clinicalData[,j])$p.value)
        pvalues <- c(pvalues, pvalues_col)
        if(min(pvalues_col)<pval.thr){
          clin_erich <- clin_erich+1
        }
      }
    }
    # Number of clinical annotations with at least one significant p-value
    line3 <- rbind(line3, clin_erich)
    
    # Total number of significant factors in all tested columns
    column <- names(pvalues)[pvalues<pval.thr]
    
    # Number of unique factors that were significant at least once
    line2<-rbind(line2, length(unique(column)))
    
    # Number of times a p-value was found significant
    signif <- length(column)
    f<-length(unique(column))   
    
    # Selectivity 
    if(signif!=0){
      line <- rbind(line,((clin_erich/signif)+(f/signif))/2)
    }else{
      line <- rbind(line,0)
    }
  }
  
  # Store and return results
  out <- data.frame(selectivity=line, nonZeroFacs=line2, total_annotations=line3)
  rownames(out) <- names(factorsList)
  return(out)
}
##Convert ENSEMBL ids to gene names
# # extract all ensembl ids
# allGenes = unique(unlist(lapply(metagenesList, function(x) rownames(x[["mRNA"]]))))
# mart = useDataset("hsapiens_gene_ensembl", useMart("ensembl"))
# # get gene ids for ensembl ids
# genSymbols = getBM(filters="ensembl_gene_id",
#                    attributes=c("ensembl_gene_id", "hgnc_symbol"),
#                    values=allGenes, mart=mart)
# # select first id if more than one is present
# genSymbols = genSymbols[!duplicated(genSymbols[,"ensembl_gene_id"]),]
# # set rownames to ens id
# rownames(genSymbols) = genSymbols[,"ensembl_gene_id"]
# write.table(genSymbols, "./data/GeneSymbols.tsv", sep = "\t", col.names = TRUE, row.names = TRUE, quote = FALSE)