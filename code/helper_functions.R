# Helper functions for T cell DI project

# Plotting theme
theme_custom <- function(base_size = 12, base_family = "Arial", ...) {
  theme(
    text = element_text(size = base_size + 2, family = base_family, face = "plain"),
    strip.background = element_blank(),
    panel.background = element_rect(fill = "transparent", color = NA),          
    plot.background = element_rect(fill = "transparent", color = NA),            
    panel.grid.major  = element_line(color = "transparent", linetype = "solid"),  
    panel.grid.minor = element_line(color = "transparent", linetype = "solid"),   
    panel.border = element_rect(color = "black", fill = NA, size = 1),            
    legend.background = element_rect(fill = "transparent", color = NA),           
    legend.box.background = element_rect(fill = "transparent", color = NA),       
    legend.key = element_rect(fill = "transparent", color = NA),                  
    axis.line = element_line(color = "transparent", size = 1),
    axis.line.x.top = element_line(color = "transparent", size = 1),
    axis.line.x.bottom = element_line(color = "transparent", size = 1),
    axis.line.y.left = element_line(color = "transparent", size =1),
    axis.line.y.right = element_line(color = "transparent", size =1),
    axis.ticks = element_line(color = "black", size = 1),
    axis.ticks.length = unit(4, "pt"),
    axis.text = element_text(color = "black", size = base_size, family = base_family),
    legend.text = element_text(color = "black", size = base_size, family = base_family),
    ...
  )
}

# Makes duplicated elements in vector unique
uniquify <- function(v, sep = '.'){
  ave(v, v, FUN = function(x) {
    if(length(x) > 1) {
      paste(x, 1:length(x), sep = sep)
    } 
    else x
  }
  )
}

# Normalize data to counts per million (CPM)
cpm.normalize <- function(df, libSize) {
  x <- as.matrix(df)
  if(missing(libSize)) {libSize <- colSums(x)}
  res <- t(t(x) * 1e6 / libSize)
  as.data.frame(res)
}

# Normalize data to transcripts per million (TPM)
tpm.normalize <- function(df, len, libSize) {
  if(sum(len == 0) > 0) warning('Length vector has 0 values. TPM set to 0.\n')
  x <- as.matrix(df) / (len / 1e3)
  x[is.na(x)] <- 0
  x[is.infinite(x)] <- 0
  if(missing(libSize)) {libSize <- colSums(x)}
  res <- t(t(x) * 1e6 / libSize)
  as.data.frame(res)
}

# Apply transformation to grouped data
group.transform <- function(df, group, FUN) {
  group.lbls <- unique(group)
  res <- vector(mode = "list")
  for (lbl in group.lbls) {
    res[[lbl]] <- FUN(df[,group==lbl])
  }
  res <- bind_cols(res)
  row.names(res) <- row.names(df)
  res
}

# Performs principal component analysis
doPCA <- function(x, topVar=500) {
  rv <- matrixStats::rowVars(as.matrix(x))
  select <- order(rv, decreasing = TRUE)[seq_len(min(topVar, length(rv)))]
  pca <- prcomp(t(x[select, ]), center = T, scale. = F)
  
  res <- vector(mode = 'list')
  res[['percentVar']] <- pca$sdev^2/sum(pca$sdev^2)*100
  res[['pcs']] <- as.data.frame(pca$x)
  res[['rotation']] <- as.data.frame(pca$rotation)
  
  res
}

# Function to load GMT files (author: Minghui Wang)
readGMT <- function(file) {
  if(!file.exists(file)) stop('File ',file,' not available\n')
  x <- readLines(file)
  n <- length(x)
  res <- list(genesets = vector(mode = "list", length = n),
              geneset.names = vector(mode = "character", length = n),
              geneset.descriptions = vector(mode = "character", length = n))
  for(i in 1:n) {
    s <- strsplit(x[i],'\t')[[1]]
    res$genesets[[i]] <- s[-c(1:2)]
    res$geneset.names[i] <- s[1]
    res$geneset.descriptions[i] <- s[2]
  }
  names(res$genesets) <- res$geneset.names
  res
}

calcPerVolunteerFCs <- function(counts, # Normalised count data
                                sampleInfo, # sample annotation data frame
                                comparison # list of vectors of comparisons for which the FCs should be calculated reference column should be in position 2
){
  # get all different volunteers
  p <- levels(as.factor(sampleInfo$volunteer))
  print(p)
  fcList <- list(n=length(p))
  for( i in 1:length(p)){
    print(p[i])
    # For each volunteers get the gene expression data
    tmpDat <- counts[,sampleInfo$sample[sampleInfo$volunteer == p[i]]]

    tmpInfo <- sampleInfo[sampleInfo$volunteer == p[i],]

    tmpFCList <- list()
    for( c in 1:length(comparison)){
      # collect the time points in the selected volunteer and calculate the fold difference
      tmpFC <- tmpDat[,as.character(tmpInfo$sample[tmpInfo$time == comparison[[c]][1]])] - tmpDat[,as.character(tmpInfo$sample[tmpInfo$time == comparison[[c]][2]])]
      tmpFCList[[c]] <- tmpFC
      names(tmpFCList)[c] <- paste(p[i],paste(comparison[[c]],collapse = "_"),collapse ="_")
    }
    # Combine all fold changes from one volunteer into a list
    tmpOut <-do.call(cbind,tmpFCList)
    fcList[[i]] <- tmpOut
  }
  # Collapse the list (i.e. each entry is a volunteer) into a matrix and return the data
  
  combined <- do.call(cbind,fcList)
  rownames(combined) <- rownames(counts)
  return(combined)
}

runGage <- function(expList, #list where each element is a vector of values on which to perform gage analysis.
                    gSet,
                    cutOff){
  require(gage)
  gageCollect <- list()
  for(i in 1:length(expList)){
    rankDat <- apply(as.matrix(expList[[i]]),2,rank)
    
    tmpGage <- gage(rankDat, gsets = gSet, ref = NULL, samp = NULL,
                    rank.test=T)
    
    tmpSig <- sigGeneSet(tmpGage,cutOff)
    gageCollect[[i]] <- list("all"=tmpGage,"sig"=tmpSig)
  }
  
  names(gageCollect) <- names(expList)
  return(gageCollect)
}
