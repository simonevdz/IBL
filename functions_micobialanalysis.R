
## remove ineffective samples and normalize
do.CSS <- function(source, title){
  x11()
  phenotypeData = AnnotatedDataFrame(source$meta)
  OTUdata = AnnotatedDataFrame(source$taxa)
  obj = newMRexperiment(source$counts, phenoData = phenotypeData, featureData = OTUdata)
  p2 = cumNormStat(obj, pFlag = TRUE, main = "Trimmed obj1")
  obj = cumNorm(obj, p = p2)
  print("next")

  settings = zigControl(tol=1e-10, maxit = 20, verbose = TRUE)
  mod = model.matrix(~pData(obj)$condition)
  res = fitZig(obj = obj, mod = mod, useCSSoffset = TRUE, control = settings)
  print(res)
  res1 <- calculateEffectiveSamples(res)
  dir.create("normalisation")
  write.table(res1, paste0("normalisation/effsamples-",title,".txt"), sep= "\t")
  
  
  Efmat <- merge(as.data.frame(source$counts), as.data.frame(res1), by="row.names")
  lower <- mean(Efmat$res1,na.rm=TRUE)
  Efmat.filt <- subset(Efmat,res1>lower)
  Effeatuers <- Efmat.filt$Row.names
  counts_EF <- subset(source$counts,rownames(source$counts) %in% Effeatuers)
  tax_EF <- subset(source$taxa,rownames(source$taxa) %in% Effeatuers)
  write.table(counts_EF,paste0("normalisation/count_matrix_EF-",title,".txt"),sep = "\t",row.names = T, col.names = NA)
  write.table(tax_EF,paste0("normalisation/taxonomy_EF-",title,".txt"),sep = "\t",row.names = T, col.names = NA)
  
  #calculating the normalization factors
  OTUdata = AnnotatedDataFrame(tax_EF)
  obj = newMRexperiment(counts_EF, phenoData = phenotypeData, featureData = OTUdata) 
  obj = cumNorm(obj, p = cumNormStatFast(obj))
  
  #Normalization export
  mat = MRcounts(obj, norm = TRUE, log = T)
  exportMat(mat, file = file.path(paste0("normalisation/CSS_EF-",title,".txt")))
  head(mat)
  
  result <- list()
  result$obj <- obj
  result$mat <- mat
  result$taxa <- tax_EF
  
  return(result)
}

## create a phyloseq object
make.phylo <- function(source,EF){
  ps.css <- phyloseq(otu_table(as.matrix(EF$mat), taxa_are_rows=T), sample_data(source$meta), tax_table(as.matrix(EF$tax)))
  dna <- Biostrings::DNAStringSet(taxa_names(ps.css))
  names(dna) <- taxa_names(ps.css)
  ps.css <- merge_phyloseq(ps.css)
  taxa_names(ps.css) <- rownames(raw.data$counts[match(rownames(EF$mat),rownames(raw.data$counts.seq)),])
  return(ps.css)
}

## get beta diversity measures
get.beta <- function(ps.css){
  beta <- list()
  
  ## PCoA
  ordu_PCOA = ordinate(ps.css, "PCoA", "bray")
  beta$pcoa <- plot_ordination(ps.css, ordu_PCOA, color="stage", shape="status") + geom_point(size=2)
  
  ## CAP
  ordu_CAP <- ordinate(ps.css,"CAP","bray", ~condition)
  beta$cap <- plot_ordination(ps.css, ordu_CAP, color="stage", shape="status") + geom_point(size=2)
  
  # NMDS
  ord.nmds.bray <- ordinate(ps.css, method="NMDS", distance="bray")
  beta$nmds <- plot_ordination(ps.css, ord.nmds.bray, color="stage", title="Bray NMDS DNA samples",shape = "status")+geom_point(size=5)
  
  
  #get bray dissimilarity matrix
  sampledf <- data.frame(sample_data(ps.css))
  dist_bray <- phyloseq::distance(ps.css, method = "bray")
  
  ## permanova of bray dissimilarity 
  if(length(unique(sample_data(ps.css)$status)) > 1 & length(unique(sample_data(ps.css)$stage)) > 1){
    beta$adonis_two <- adonis2(dist_bray~stage+status, data=sampledf) 
    beta$adonis_int <- adonis2(dist_bray~stage*status, data=sampledf) 
  }
  if(length(unique(sample_data(ps.css)$status)) > 1){
    beta$adonis_stat <- adonis2(dist_bray~status, data=sampledf) 
  }
  if(length(unique(sample_data(ps.css)$stage)) > 1){
    beta$adonis_stage <- adonis2(dist_bray~stage, data=sampledf)
  }
  
  ## bar plot of phyla 
  top20 <- names(sort(taxa_sums(ps.css), decreasing=TRUE))[1:20]
  ps.top20 <- transform_sample_counts(ps.css, function(x) x/sum(x))
  ps.top20 <- prune_taxa(top20, ps.top20)
  par(mfrow=c(5,1))
  beta$abundance <- ggplot(psmelt(ps.top20), aes_string(x = "stage", y="Abundance", fill = "Phylum")) + 
    geom_bar(stat = "identity", position = "fill") + 
    theme(axis.text.x = element_text(angle = -90, hjust = 0)) +
    ggtitle("DNA relative abundance (50)")
  
  return(beta)
}

## make fitzig model
make.model <- function(obj, C, contrast) {
  settings = zigControl(maxit = 30, verbose = TRUE)
  mod = model.matrix(~0+C)
  res = fitZig(obj = obj, mod = mod, useCSSoffset = T, control = settings)
  zigFit = res@fit
  finalMod= res@fit[["cov.coefficients"]]

  
  contrast.matrix = makeContrasts(contrasts = contrast, levels = finalMod)
  fit2 = contrasts.fit(zigFit, contrast.matrix)
  fit3 = eBayes(fit2)
  res <- topTable(fit3,coef=1,adjust="BH",n=Inf)
  write.table(res, paste0("contrasts/",contrast,"-DNA.txt"), sep="\t")
  return(res)
} 

## calculate the fold changes for all contrasts
calculate.FC <- function(obj,contrasts,sig,fc){
  FC <- list()
  
  # create family aggregated object with new normalisation factors
  counts <- as.data.frame(obj@assayData[["counts"]])
  tax <- as.data.frame(obj@featureData@data)
  merge <- aggregate(counts,by = list(tax$Family,tax$Phylum), FUN=sum)
  counts <- merge[,-c(1,2)]
  tax <- merge[,1:2]
  rownames(counts) <- paste("agg",1:nrow(merge),sep="")
  rownames(tax) <- paste("agg",1:nrow(merge),sep="")
  colnames(tax) <- c("Family","Phylum")
  pheno <- AnnotatedDataFrame(obj@phenoData@data)
  OTU <- AnnotatedDataFrame(tax)
  obj1 <- newMRexperiment(counts, phenoData = pheno,featureData = OTU)
  obj1 <- cumNorm(obj1, p = cumNormStatFast(obj1))
  
  for(i in contrasts){
    #make model for contrast and filter on p-val and logFC
    model <- make.model(obj,pData(obj)$condition,i) 
    FC$raw[[i]] <- model
    FC$filt[[i]] <- model %>%
      filter(adj.P.Val < sig & (logFC > fc | logFC < -fc))
    
    #make family aggregated model and filter on p-val and logFC
    model1 <- make.model(obj1,pData(obj1)$condition,i)
    model1$Family <- tax[match(rownames(model1),rownames(tax)),"Family"]
    model1$Phylum <- tax[match(rownames(model1),rownames(tax)),"Phylum"]
    model1 <- model1 %>%
      filter(adj.P.Val<sig & (logFC > fc | logFC < -fc))
    if(nrow(model1)>0){
      FC$agg[[i]] <- model1
    }
  }
  return(FC)
}

## create barplot of fold change
create.bar <- function(FC){
  barplots <- list()
  for(i in names(FC)){
    FC[[i]]$Phylum <- factor(FC[[i]]$Phylum,levels=unique(FC[[i]]$Phylum))
    g <- ggplot(data=drop_na(FC[[i]]), aes(x=Family, y=logFC, fill=Phylum)) + 
      geom_bar(stat="identity", position=position_dodge(),width = 0.75) + 
      labs(title=i) + 
      facet_grid(Phylum~., scales = "free_y", space = "free_y",switch = "y") + 
      theme(strip.placement = "outside", strip.background = element_rect(fill = "white"),  axis.title = element_blank()) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + coord_flip()
    barplots[[i]] <- g 
  }
  return(barplots)
}

## create the rarefaction curves to identify proper depth
calculate_rarefaction_curves <- function(psdata, measures, depths) {
  require('plyr') # ldply
  require('reshape2') # melt
  
  estimate_rarified_richness <- function(psdata, measures, depth) {
    if(max(sample_sums(psdata)) < depth) return()
    psdata <- prune_samples(sample_sums(psdata) >= depth, psdata)
    
    rarified_psdata <- rarefy_even_depth(psdata, depth, verbose = FALSE)
    
    alpha_diversity <- estimate_richness(rarified_psdata, measures = measures)
    
    # as.matrix forces the use of melt.array, which includes the Sample names (rownames)
    molten_alpha_diversity <- melt(as.matrix(alpha_diversity), varnames = c('Sample', 'Measure'), value.name = 'Alpha_diversity')
    
    molten_alpha_diversity
  }
  
  names(depths) <- depths # this enables automatic addition of the Depth to the output by ldply
  rarefaction_curve_data <- ldply(depths, estimate_rarified_richness, psdata = psdata, measures = measures, .id = 'Depth', .progress = ifelse(interactive(), 'text', 'none'))
  
  # convert Depth from factor to numeric
  rarefaction_curve_data$Depth <- as.numeric(levels(rarefaction_curve_data$Depth))[rarefaction_curve_data$Depth]
  
  rarefaction_curve_data
}

## plot the rarefaction curves
plot.rarefaction <- function(rare){
  ggplot(
    data = rare,
    mapping = aes(
      x = Depth,
      y = Alpha_diversity_mean,
      ymin = Alpha_diversity_mean - Alpha_diversity_sd,
      ymax = Alpha_diversity_mean + Alpha_diversity_sd,
      group = Sample
    )
  ) + geom_line(
  ) + geom_pointrange(
  ) + facet_wrap(
    facets = ~ Measure,
    scales = 'free_y'
  )
}

## get anova results for different alpha diversity measures
get.alpha <- function(data, measures,title,q){
  
  alpha.plot <- plot_richness(data, x="status", measures= measures, color="stage") + geom_boxplot(alpha=0.5,position = "identity")
  alpha.table <- alpha.plot[['data']]
  write.csv(alpha.table,paste0("alpha_div/alpha_table_",title,".csv"))
  
  results <- list()
  for(i in measures){
    table <- alpha.table %>% filter(variable==i)
    
    if(length(unique(table$stage)) == 1){
      model <- aov(value ~ status, data = table) #one way anonva for single stage
    }
    else{
      model <- aov(value ~ stage + status, data=table) #two way anova addititve
      model1 <- aov(value ~ stage * status, data=table) # two way anova with interaction
      if(summary(model1)[[1]][["Pr(>F)"]][3] < 0.05){
        model <- model1
      }
    }
    results[[i]]$model <- model
    results[[i]]$tukey <- TukeyHSD(model)
    
    if(q == T){   ## quality inspection
      x11()
      par(mfrow=c(2,2))
      plot(model)
      par(mfrow=c(1,1))
    }
    
  }
  results$table <- alpha.table
  return(results)
}

## plot fold changes per ASV
plot_phyla <- function(table,title){
  x11() 
  ggplot(data=drop_na(table), aes(x=reorder(asv,-logFC), y=logFC, fill=Phylum, group=yaxis)) + 
    geom_bar(stat="identity", position=position_dodge(),width = 0.75) + 
    labs(title=title, xlba="") + 
    facet_grid(~Phylum, scales = "free_x", space = "free_x",switch = "x")  + 
    theme(strip.placement = "outside", strip.background = element_rect(fill = "white")) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

## plot fold change per ASV wihtin given phyla
plot_zoom <- function(table, phylum, title) {
  x11() 
  ggplot(data=drop_na(table[table$Phylum == phylum,]), aes(x=reorder(asv,-logFC), y=logFC, fill=Family)) + 
    geom_bar(stat="identity", position=position_dodge(),width = 0.75) + 
    labs(title=title, xlba="") + 
    facet_grid(~Family, scales = "free_x", space = "free_x",switch = "x")  + 
    theme(strip.placement = "outside", strip.background = element_rect(fill = "white"),  axis.title = element_blank()) +
    theme(axis.title.x=element_blank(),axis.text.x=element_blank(), axis.ticks.x=element_blank())
}

## add taxonomy to contrast files
get_tax <- function(table, tax, selection,a,b){
  x <- table[rownames(table) %in% selection,]
  x$Family <- tax[match(rownames(x),rownames(tax)),"Family"]
  x$Phylum <- tax[match(rownames(x),rownames(tax)),"Phylum"]
  x$asv <- rownames(x)
  x$yaxis <- a
  x[x$logFC > 0,]$yaxis <- b
  return(x)
}

## export relative abundance files 
export_abundance <- function(nuc.acid,ps,ext){
  abund <- list()
  
  # samplewise 
  for(i in c("Phylum","Family","Genus")){
    counts <- as.data.frame(ps@otu_table)
    counts[,i] <- nuc.acid[["taxa"]][match(rownames(counts),rownames(nuc.acid[["taxa"]])),i]  %>% replace(is.na(.), "Unknown")
    taxcount <- aggregate(counts[,1:ncol(counts)-1], by = list(counts[,i]),FUN=sum)
    taxcount <- as.data.frame(taxcount)
    rownames(taxcount) <- taxcount$Group.1
    taxcount <- taxcount[,-1]
    taxcount <- data.frame(apply(t(taxcount),1,function(x) (x/sum(x))*100))
    abund[[i]] <-taxcount
    write.table(taxcount,paste0("rel_abundance/relative_abundance_",i,"_",ext,".txt"))
  }
  
  # create count matrix aggregegated by condition
  condition.count <- as.data.frame(t(ps@otu_table))
  condition.count$condition <- nuc.acid[["meta"]][match(rownames(condition.count),rownames(nuc.acid[["meta"]])),"condition"]
  condition.count <- aggregate(condition.count[,1:ncol(condition.count)-1],by=list(condition.count$condition),FUN=sum)
  rownames(condition.count) <- condition.count$Group.1
  condition.count <- t(condition.count[,-1])

  for(i in c("Phylum","Family","Genus")){
    counts <- as.data.frame(condition.count)
    counts[,i] <- nuc.acid[["taxa"]][match(rownames(counts),rownames(nuc.acid[["taxa"]])),i]  %>% replace(is.na(.), "Unknown")
    taxcount <- aggregate(counts[,1:ncol(counts)-1], by = list(counts[,i]),FUN=sum)
    taxcount <- taxcount[,1:7]
    rownames(taxcount) <- taxcount$Group.1
    taxcount <- taxcount[,-1]
    taxcount <- data.frame(apply(t(taxcount),1,function(x) (x/sum(x))*100))
    abund$condition[[i]] <-taxcount
    write.table(taxcount,paste0("rel_abundance/relative_abundance_condition_",i,"_",ext,".txt"))
  }
  
  return(abund)
}

## export differential abundances
export_diff <- function(source, names, FCs, title, frac, overwrite){
  for(i in c("Phylum", "Family", "Genus")){
    df <- as.data.frame(unique(source$taxa[[i]]))
    df <- replace(df,is.na(df), "unknown")
    colnames(df) <- c(i)
    rownames(df) <- df[[i]]
    for(j in names){
      x <- FCs[[j]]
      x <- x[frac[[j]],]
      x[[i]] <- source$taxa[match(x$names,rownames(DNA$taxa)),i]
      x <- as.data.frame(table(x[[i]]))
      df[[j]] <- x[match(df[[i]],x$Var1),"Freq"]
    }
    
    #add overlap
    a <- paste0("overlap_spermo_",title)
    b <- paste0("modern_spermo_",title)
    x <- FCs[[b]][frac[[a]],]
    x[[i]] <- DNA$taxa[match(x$names,rownames(DNA$taxa)),i]
    x <- as.data.frame(table(x[[i]]))
    df$overlap_spermo <- x[match(df[[i]],x$Var1),"Freq"]
    
    a <- paste0("overlap_rhizo_",title) 
    b <- paste0("modern_rhizo_",title)
    x <- FCs[[b]][frac[[a]],]
    x[[i]] <- DNA$taxa[match(x$names,rownames(DNA$taxa)),i]
    x <- as.data.frame(table(x[[i]]))
    df$overlap_rhizo <- x[match(df[[i]],x$Var1),"Freq"]
    
    df <- replace(df,is.na(df),0)
    df <- df[apply(df[,-1], 1, function(x) !all(x==0)),]
    
    if(i == "Family"){
      df$Phylum <- source$taxa[match(df$Family,source$taxa$Family),"Phylum"]
    }
    
    write.xlsx(df,file= paste0("differential_ASVs/venn_",i,"_",title,".xlsx"), overwrite = overwrite)
  }
}

## get t statistics for relative abundance
get.tstat <- function(df,meta,level,p){
  asvs <- list()
  for(i in unique(df[[level]])){
    a1 <- na.omit(unlist(c(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="modern_spermosphere",])])))
    b1 <- na.omit(unlist(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="modern_rhizosphere",])]))
    a2 <- na.omit(unlist(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="wild_spermosphere",])]))
    b2 <- na.omit(unlist(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="wild_rhizosphere",])]))
    a3 <- na.omit(unlist(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="bulk_bulkS",])]))
    b3 <- na.omit(unlist(df[df[[level]]==i,colnames(df)%in%rownames(meta[meta$condition=="bulk_bulkR",])]))
    if(length(a1) > 1 & length(b1)>1){
      x1 <- t.test(a1,b1)
      if(is.na(x1$p.value)){
        # nothing
      }
      else if(x1$p.value < p){
        asvs$modern[[i]] <- x1
      }
    }
    if(length(a2) > 1 & length(b2)>1){
      x2 <- t.test(a2,b2)
      if(is.na(x2$p.value)){
        #nothing
      }
      else if(x2$p.value < p){
        asvs$wild[[i]] <- x2
      }
    }
    if(length(a1) > 1 & length(a2)>1){
      x3 <- t.test(a1,a2)
      if(is.na(x3$p.value)){
        #nothing
      }
      else if(x3$p.value < p){
        asvs$sperm[[i]] <- x3
      }
    }
    if(length(b1) > 1 & length(b2)>1){
      x4 <- t.test(b1,b2)
      if(is.na(x4$p.value)){
        #nothing
      }
      else if(x4$p.value < p){
        asvs$rhiz[[i]] <- x4
      }
    }
    if(length(a3) > 1 & length(b3)>1){
      x5 <- t.test(a3,b3)
      if(is.na(x5$p.value)){
        # nothing
      }
      else if(x5$p.value < p){
        asvs$bulk[[i]] <- x5
      }
    }
  }
  return(asvs)
}
