
#Overal WGCNA correlation network#####
rm(list = ls())

wkdir <- 'C:/Users/evely/Google Drive/summer intern/Preterm birth metabolomics project'

setwd(wkdir)

pd = readRDS("clinic.rds")
meta = readRDS("quant.metabolitemeasurements.rds")

newdat$Label = meta$Label

ctrldat <- subset(newdat, Label == 'Control')
preedat <- subset(newdat, Label != 'Control')

overaldat <- newdat

#Choose the soft-thresholding power#
softval <- NULL

topquantile <- 0.95
internum <- 1000
largevexsize <- 15
largeedgsize <- 15
absolutecut <- 0.4

tag <- 'Preterm to Control'

library(WGCNA)

choosepower <- function(datres = ctrldat, rsqcutline = 0.7){
  
  library(WGCNA)
  
  datExpr <- datres[-1]
  
  powers <- c(c(1:20))
  
  sft <- pickSoftThreshold(datExpr, powerVector = powers, RsquaredCut = rsqcutline, verbose = 5)
  
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], 
       xlab = "Soft Threshold (power)", 
       ylab = "Scale Free Topology Model Fit, signed R^2", type = "n", main = paste("Scale independence"))
  
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*(sft$fitIndices[,2]), labels = powers, col = 'red')
  
  abline(h = rsqcutline, col = 'red')
  
  
  
  plot(sft$fitIndices[,1], sft$fitIndices[,5], 
       xlab = "Soft Threshold (power)", ylab = "Mean Connectivity", type = "n", 
       main = "Mean connectivity")
  
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels = powers, col = 'red')
  
  return(sft)
  
}

if(is.null(softval)){
  
  ctrlpower <- choosepower(rsqcutline = 0.8)
  
  casepower <- choosepower(datres = preedat, rsqcutline = 0.8)
  
  finalsft <- max(ctrlpower$powerEstimate, casepower$powerEstimate)
  
}else{
  finalsft <- softval
}

##finalsft = 8

makemain <- function(datres = ctrldat, sftpower = finalsft, mergesimilar = TRUE, mergecut = 0.25){
  
  library(WGCNA)
  
  datExpr <- datres[-334]
  probes <- names(datExpr)
  
  #Calculate TOM
  adjacency <- adjacency(datExpr, power = sftpower)
  TOM <- TOMsimilarity(adjacency)
  dimnames(TOM) <- list(probes, probes)
  
  dissTOM <- 1 - TOM
  
  # Call the hierarchical clustering function
  geneTree <- hclust(as.dist(dissTOM), method = "average")
  
  #Module identification using dynamic tree cut
  dynamicMods <- cutreeDynamic(dendro = geneTree, distM = dissTOM,
                               deepSplit = 2, pamRespectsDendro = FALSE,
                               minClusterSize = 30)
  
  dynamicColors <- labels2colors(dynamicMods)
  oricolors <- dynamicColors
  finalColors <- dynamicColors
  
  # Calculate eigengenes
  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  
  if(mergesimilar == TRUE){
    
    #Calculate dissimilarity of module eigengenes
    MEDiss <- 1-cor(MEs)
    
    MEDissThres <- mergecut
    
    # Call an automatic merging function
    merge <- mergeCloseModules(datExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
    finalColors <- merge$colors
    
  }
  
  plotDendroAndColors(geneTree, finalColors, 
                      c("Tree cut"), 
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05)
  
  if(mergesimilar == TRUE){
    res <- list(tom = TOM, melist = MEList, mergedmelist = merge, dynamicColors = oricolors, 
                mergedColors = finalColors)
  }else{
    res <- list(tom = TOM, melist = MEList, dynamicColors = oricolors)
  }
  
  return(res)
}

ctrlmain <- makemain()
casemain <- makemain(datres = preedat)

#Own network organization####
ctrlcolors <- ctrlmain$mergedColors
names(ctrlcolors) <- row.names(ctrlmain$tom)

casecolors <- casemain$mergedColors
names(casecolors) <- row.names(casemain$tom)


orgnetwork <- function(TOM = ctrlmain$tom, netdata = ctrldat, colors = ctrlcolors, fixcolor = NULL, 
                       quantilecut = topquantile, abscut = 0.5, 
                       largenodesize = largevexsize, largeedgesize = largeedgsize, nodesize = 5, 
                       edgesize = 2, 
                       disfun = ctrldisfun, intnum = internum, pdfprefix = 'ctrlquant', plot = TRUE, 
                       plotonly = TRUE, calconnec = NULL){
  
  datExpr <- netdata[-334]
  
  probes <- names(datExpr)
  
  modTOM <- TOM
  modProbes <- probes
  
  dimnames(modTOM) <- list(modProbes, modProbes)
  
  if(!is.null(abscut)){
    weightcutoff <- abscut
  }else{
    library(corpcor)
    tomvalues <- sm2vec(modTOM)
    weightcutoff <- as.numeric(quantile(tomvalues, quantilecut))
  }
  
  unicolors <- unique(colors)
  modulecolors <- unicolors[unicolors != 'grey']
  
  i <- 1
  
  modulenodes <- data.frame(nodeName = character(), 
                            nodeColor = character(), 
                            stringsAsFactors = FALSE)
  for(i in 1:length(modulecolors)){
    
    modulecolor <- modulecolors[i]
    modulegenes <- names(colors)[colors == modulecolor]
    subTOM <- modTOM[modulegenes, modulegenes]
    
    if(nrow(subTOM) <= 0){
      next()
    }
    
    subcyt <- exportNetworkToCytoscape(subTOM, 
                                       edgeFile = NULL, nodeFile = NULL, 
                                       weighted = TRUE, 
                                       threshold = weightcutoff, 
                                       nodeNames = row.names(subTOM))
    
    subnodedata <- subcyt$nodeData
    subnodedata <- subnodedata[c(1)]
    if(nrow(subnodedata) == 0){
      next()
    }
    subnodedata$nodeColor <- modulecolor
    subnodedata$nodeName <- as.character(subnodedata$nodeName)
    
    modulenodes <- rbind(modulenodes, subnodedata)
    
    
  }
  
  if(is.null(edgesize)){
    edgeweight <- TRUE
  }else{
    edgeweight <- FALSE
  }
  
  cyt <- exportNetworkToCytoscape(modTOM, 
                                  edgeFile = NULL, 
                                  nodeFile = NULL, 
                                  weighted = edgeweight, 
                                  threshold = weightcutoff, 
                                  nodeNames = modProbes)
  edgedata <- cyt$edgeData
  nodedata <- cyt$nodeData
  edgedata <- edgedata[c(1, 2, 3)]
  nodedata <- nodedata[c(1)]
  
  edgedata$fromNode <- as.character(edgedata$fromNode)
  edgedata$toNode <- as.character(edgedata$toNode)
  nodedata$nodeName <- as.character(nodedata$nodeName)
  
  greynodes <- names(colors)[colors == 'grey']
  
  if(length(greynodes) > 0){
    candnodes <- unique(c(greynodes, modulenodes$nodeName))
    edgedata <- subset(edgedata, (fromNode %in% candnodes) & (toNode %in% candnodes))
    nodedata <- subset(nodedata, nodeName %in% unique(c(edgedata$fromNode, edgedata$toNode)))
  }
  
  
  
  library(igraph)
  
  inet <- graph_from_data_frame(d = edgedata, directed = FALSE, vertices = nodedata)
  
  vertexattr <- vertex.attributes(inet)
  edgeattr <- edge.attributes(inet)
  
  if(is.null(calconnec)){
    
    calconn <- function(edgedat = edgedata){
      dat1 <- edgedat
      dat2 <- edgedat[c('toNode', 'fromNode', 'weight')]
      names(dat1) <- names(dat2) <- c('node1', 'node2', 'weight')
      dat <- rbind(dat1, dat2)
      dat <- unique(dat)
      row.names(dat) <- 1:nrow(dat)
      dat <- dat[c('node1', 'weight')]
      names(dat)[1] <- 'node'
      
      library(plyr)
      
      calsum <- function(block){
        nodename <- unique(block[,1])
        numblock <- block[c(2)]
        numsum <- colSums(numblock)
        resblock <- data.frame(node = nodename, conn = numsum, stringsAsFactors = FALSE)
        row.names(resblock) <- 1:nrow(resblock)
        return(resblock)
        
      }
      
      dat <- dat[order(dat$node),]
      row.names(dat) <- 1:nrow(dat)
      conres <- ddply(.data = dat, .variables = c('node'), .fun = calsum)
      
      return(conres)
    }
    
    nodeconn <- calconn()
    connmax <- max(nodeconn$conn)
    nodeconn$normconn <- nodeconn$conn/connmax
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
    
  }else{
    nodeconn <- calconnec
    
    normconnvals <- nodeconn$normconn
    names(normconnvals) <- nodeconn$node
  }
  
  
  
  asscolors <- function(wholelabels = colors, nodedat = nodedata){
    
    library(scales)
    
    netlabels <- wholelabels[nodedat$nodeName]
    sigcolors <- names(table(netlabels))[order(-table(netlabels))]
    sigcolors  <- sigcolors[sigcolors != 'grey']
    newcolors <- hue_pal()(length(sigcolors))
    netlabels[netlabels == 'grey'] <- '#FFFFFF'
    
    for(i in 1:length(sigcolors)){
      sigcolor <- sigcolors[i]
      netlabels[netlabels == sigcolor] <- newcolors[i]
      
    }
    
    return(netlabels)
    
  }
  
  if(!is.null(fixcolor)){
    labels <- fixcolor[nodedata$nodeName]
  }else{
    if('grey' %in% unique(colors) & length(unique(colors)) == 1){
      labels <- rep('#FFFFFF', nrow(nodedata))
    }else{
      labels <- asscolors()
    }
  }
  
  colormapping <- colors[names(labels)]
  colormapping <- data.frame(ori = colormapping, new = labels, stringsAsFactors = FALSE)
  colormapping <- unique(colormapping)
  row.names(colormapping) <- 1:nrow(colormapping)
  
  vertexattr$nodeColors <- as.vector(labels)
  vsize <- log10(normconnvals)
  
  esize <- log10(edgeattr$weight)
  
  library(scales)
  vsize <- rescale(vsize, to = c(1, largenodesize))
  
  V(inet)$size <- vsize
  
  esize <- rescale(esize, to = c(1, largeedgsize))
  
  E(inet)$width <- esize
  
  if(!is.null(nodesize)){
    V(inet)$size <- nodesize
  }
  
  if(edgeweight == FALSE){
    E(inet)$width <- edgesize
  }
  
  
  set.seed(7)
  
  if(is.null(disfun)){
    netstyle <- layout_with_fr(inet, coords = NULL, dim = 2, niter = intnum, grid = c("nogrid"))
    
    saveRDS(netstyle, paste0(pdfprefix, '_interation', intnum, '_disfun.rds'))
  }else{
    netstyle <- disfun
  }
  
  if(plot == TRUE){
    
    pdf(paste0(pdfprefix, "wgcna.pdf"), height=10, width=10, useDingbats=FALSE)
    print(
      
      plot(inet, vertex.label = NA, vertex.color = vertexattr$nodeColors, 
           edge.color = '#000000', 
           layout = netstyle)
      
      
    )
    dev.off()
    
  }
  
  
  if(plotonly == FALSE){
    res <- list(inet = inet, nodenormconn = normconnvals, nodeconn = nodeconn, 
                nodecolor = labels, edgeweight = edgeattr$weight, colormap = colormapping)
  }else{
    res <- NULL
  }
  
  return(res)
  
  
}

ctrlquantnetres <- orgnetwork(disfun = NULL, plotonly = FALSE, pdfprefix = 'ctrlquant', edgesize = NULL, 
                              abscut = absolutecut, plot = FALSE)

preequantnetres <- orgnetwork(TOM = casemain$tom, netdata = preedat, colors = casecolors, 
                              disfun = NULL, pdfprefix = 'preequant', plotonly = FALSE, edgesize = NULL, 
                              abscut = absolutecut, plot = FALSE)

#Export to Cytoscape####
getclass <- function(nodename = ctrlcyto$nodes$source){
  
  unknowngroupidx <- grep(pattern = 'Unknown ', x = nodename)
  normgroupidx <- setdiff(1:length(nodename), unknowngroupidx)
  
  normgroup <- nodename[normgroupidx]
  unknowngroup <- nodename[unknowngroupidx]
  
  normclasses <- gsub(pattern = ' .*$', replacement = '', x = normgroup)
  unknownclasses <- gsub(pattern = 'Unknown ', replacement = '', x = unknowngroup)
  unknownclasses <- gsub(pattern = ' .*$', replacement = '', x = unknownclasses)
  
  nodeclasses <- data.frame(source = c(normgroup, unknowngroup), 
                            classes = c(normclasses, unknownclasses), 
                            idx = c(normgroupidx, unknowngroupidx), stringsAsFactors = FALSE)
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses$labels <- nodeclasses$source
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = 'Unknown ', replacement = 'u', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses$labels[unknowngroupidx] <- gsub(pattern = ';.*$', replacement = '', 
                                              x = nodeclasses$labels[unknowngroupidx])
  nodeclasses <- nodeclasses[order(nodeclasses$idx),]
  nodeclasses <- nodeclasses[-grep(pattern = 'idx', x = names(nodeclasses))]
  
  return(nodeclasses)
  
}

ctrlclass <- getclass(nodename = names(ctrlquantnetres$nodecolor))
preeclass <- getclass(nodename = names(preequantnetres$nodecolor))



classcolors <- data.frame(classes = unique(c(ctrlclass$classes, preeclass$classes)), 
                          stringsAsFactors = FALSE)
classcolors$classcolors <- hue_pal()(nrow(classcolors))

exportcyto <- function(res = ctrlquantnetres, writefile = TRUE, 
                       prefix = 'Control_adj_screen', 
                       classcolor = classcolors){
  
  nodes <- data.frame(source = names(res$nodenormconn), normconn = as.vector(res$nodenormconn), 
                      color = as.vector(res$nodecolor), stringsAsFactors = FALSE)
  edges <- as_edgelist(res$inet)
  edges <- as.data.frame(edges, stringsAsFactors = FALSE)
  names(edges) <- c('source', 'target')
  edges$weight <- res$edgeweight
  
  sourcecolor <- nodes[c('source', 'color')]
  targetcolor <- sourcecolor
  names(targetcolor) <- c('target', 'targetcolor')
  names(sourcecolor) <- c('source', 'sourcecolor')
  
  edges <- merge(edges, sourcecolor, by = 'source')
  edges <- merge(edges, targetcolor, by = 'target')
  edges$finalcolor <- '#FFFFFF'
  edges$finalcolor[edges$sourcecolor == edges$targetcolor] <- 
    edges$sourcecolor[edges$sourcecolor == edges$targetcolor]
  
  nodeclass <- getclass(nodename = nodes$source)
  nodes <- merge(nodes, nodeclass, by = c('source'))
  nodes <- merge(nodes, classcolors, by = c('classes'))
  nodes <- nodes[c('source', 'normconn', 'color', 'labels', 'classes', 'classcolors')]
  
  sourceclass <- nodeclass[c('source', 'labels')]
  names(sourceclass) <- c('source', 'sourcelabel')
  targetclass <- sourceclass
  names(targetclass) <- c('target', 'targetlabel')
  edges <- merge(targetclass, edges, by = c('target'))
  edges <- merge(sourceclass, edges, by = c('source'))
  edges <- edges[c('sourcelabel', 'targetlabel', 'weight', 'finalcolor')]
  names(edges) <- c('source', 'target', 'weight', 'finalcolor')
  
  nodes <- nodes[c('labels', 'normconn', 'source', 'classes', 'classcolors', 'color')]
  names(nodes) <- c('source', 'normconn', 'oriname', 'classes', 'classcolors', 'modulecolor')
  
  cytores <- list(nodes = nodes, edges = edges)
  
  if(writefile == TRUE){
    
    write.table(nodes, paste0(prefix, '_WGCNA_nodes.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
    write.table(edges, paste0(prefix, '_WGCNA_edges.txt'), sep = '\t', row.names = FALSE, 
                quote = FALSE)
  }
  
  return(cytores)
  
}

ctrlcyto <- exportcyto(prefix = 'Control_limmaadj_screen_0.8', 
                       writefile = TRUE)
preecyto <- exportcyto(res = preequantnetres, 
                       prefix = 'Preterm_limmaadj_screen_0.8', 
                       writefile = TRUE)



ctrlCARlipids <- ctrlcyto$nodes[grep(pattern = 'CAR\\(', x = ctrlcyto$nodes$source),]
preeCARlipids <- preecyto$nodes[grep(pattern = 'CAR\\(', x = preecyto$nodes$source),]

allCARlipids <- unique(ctrlCARlipids$source, preeCARlipids$source)

#metabolite class enrichment#####
allclasses <- getclass(nodename = colnames(newdat)[1:(ncol(newdat)-1)])

colordic <- c('red', 'cyan')
names(colordic) <- c('#F8766D', '#00BFC4')

moduleenrich <- function(groupnet = preecyto$nodes, backclass = allclasses, 
                         colordict = colordic){
  
  modulenames <- unique(groupnet$modulecolor)
  
  
  classenrich <- function(module = moduleset, back = backclass){
    modulecounts <- table(module$classes)
    backcounts <- table(back$classes)
    modulesum <- sum(modulecounts)
    backsum <- sum(backcounts)
    for(j in 1:length(modulecounts)){
      classname <- names(modulecounts)[j]
      moduleclasscount <- modulecounts[classname]
      backclasscount <- backcounts[classname]
      
      a11 <- moduleclasscount
      a12 <- backclasscount
      a21 <- modulesum - moduleclasscount
      a22 <- backsum - backclasscount
      fishermat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisherres <- fisher.test(fishermat, alternative = 'greater')
      fisherp <- fisherres$p.value
      
      if(j == 1){
        fisherps <- fisherp
      }else{
        fisherps <- c(fisherps, fisherp)
      }
      
      
    }
    
    names(fisherps) <- names(modulecounts)
    
    return(fisherps)
    
  }
  
  mod_phe_link <- function(modcolors = moduleset, modulelogps = moduleps, textsize = 50, 
                           modulecolorname = figuremodname){
    
    modclasscolors <- modcolors[c('classes', 'classcolors')]
    modclasscolors <- unique(modclasscolors)
    names(modclasscolors)[1] <- 'classname'
    modulelogps <- merge(modulelogps, modclasscolors, by = c('classname'))
    modulelogps <- modulelogps[order(modulelogps$logp),]
    modulelogps$classname <- factor(modulelogps$classname, levels = modulelogps$classname, 
                                    ordered = TRUE)
    modulelogps$classcolors <- factor(modulelogps$classcolors, levels = modulelogps$classcolors, 
                                      ordered = TRUE)
    row.names(modulelogps) <- 1:nrow(modulelogps)
    
    library(ggplot2)
    
    p <- ggplot(modulelogps, aes(x = classname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = modulelogps$classcolors) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in module ', modulecolorname)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = textsize)) + 
        theme(axis.text.y = element_text(size = textsize)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
  }
  
  for(i in 1:length(modulenames)){
    modulename <- modulenames[i]
    figuremodname <- as.vector(colordict[modulename])
    
    moduleset <- subset(groupnet, modulecolor == modulename)
    
    moduleps <- classenrich()
    
    moduleps <- data.frame(classname = names(moduleps), fisherp = as.vector(moduleps), 
                           stringsAsFactors = FALSE)
    moduleps$logp <- -log2(moduleps$fisherp)
    moduleps <- moduleps[order(-moduleps$logp),]
    moduleps <- subset(moduleps, fisherp < 0.05)
    
    if(nrow(moduleps) < 1){
      next()
    }
    
    mod_phe_link()
    
    print(moduleps)
  }
  
}

moduleenrich(groupnet = ctrlcyto$nodes)

colordic <- c('red')
names(colordic) <- c('#F8766D')

moduleenrich(groupnet = preecyto$nodes, colordict = colordic)

relatetopquantilemodules <- function(ctrlquantnet = ctrlquantnetres, casequantnet = preequantnetres, 
                                     titlesufix = tag, 
                                     ctrlpdcolornames = 
                                       c('red'), 
                                     preepdcolornames = 
                                       c('red', 'cyan')){
  
  vctrl <- V(ctrlquantnet$inet)
  vcase <- V(casequantnet$inet)
  allgenes <- unique(c(names(vctrl), names(vcase)))
  
  ctrlnodecolors <- ctrlquantnet$nodecolor
  casenodecolors <- casequantnet$nodecolor
  
  ctrlmodulecount <- table(ctrlnodecolors)
  casemodulecount <- table(casenodecolors)
  
  ctrlmodulecount <- ctrlmodulecount[ctrlquantnet$colormap$new]
  casemodulecount <- casemodulecount[casequantnet$colormap$new]
  
  names(ctrlmodulecount) <- ctrlpdcolornames
  names(casemodulecount) <- preepdcolornames
  
  
  ctrluninodecolors <- unique(ctrlnodecolors)
  caseuninodecolors <- unique(casenodecolors)
  
  pTable <- matrix(0, nrow = length(ctrluninodecolors), ncol = length(caseuninodecolors))
  rownames(pTable) <- ctrluninodecolors
  colnames(pTable) <- caseuninodecolors
  countTable <- pTable
  
  for(i in 1:length(ctrluninodecolors)){
    ctrluninodecolor <- ctrluninodecolors[i]
    for(j in 1:length(caseuninodecolors)){
      caseuninodecolor <- caseuninodecolors[j]
      
      ctrlset <- names(ctrlnodecolors)[ctrlnodecolors == ctrluninodecolor]
      caseset <- names(casenodecolors)[casenodecolors == caseuninodecolor]
      unionset <- union(ctrlset, caseset)
      interset <- intersect(ctrlset, caseset)
      
      mat <- c(length(setdiff(allgenes, unionset)), 
               length(setdiff(ctrlset, caseset)), 
               length(setdiff(caseset, ctrlset)), 
               length(interset))
      mat <- matrix(mat, nrow = 2)
      
      fisherp <- fisher.test(mat, alternative = "greater")$p.value
      overlap <- length(interset)
      
      pTable[i, j] <- fisherp
      countTable[i, j] <- overlap
      
    }
    
  }
  
  logpTable <- -log10(pTable)
  #Truncate p values smaller than 10^{-50} to 10^{-50}
  logpTable[is.infinite(logpTable)] <- 1.3*max(logpTable[is.finite(logpTable)]);
  logpTable[logpTable > 50 ] <- 50
  
  # Marginal counts (really module sizes)
  ctrlTotals <- apply(countTable, 1, sum)
  caseTotals <- apply(countTable, 2, sum)
  
  ctrlname <- sub(pattern = '^.*to ', replacement = '', x = titlesufix)
  casename <- sub(pattern = ' to.*$', replacement = '', x = titlesufix)
  
  row.names(logpTable) <- row.names(pTable) <- row.names(countTable) <- ctrlpdcolornames
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- preepdcolornames
  
  ctrlmodulecount <- ctrlmodulecount[order(-as.vector(ctrlmodulecount))]
  casemodulecount <- casemodulecount[order(-as.vector(casemodulecount))]
  
  logpTable <- logpTable[names(ctrlmodulecount), names(casemodulecount)]
  pTable <- pTable[names(ctrlmodulecount), names(casemodulecount)]
  countTable <- countTable[names(ctrlmodulecount), names(casemodulecount)]
  
  logpTable <- matrix(logpTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  pTable <- matrix(pTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  countTable <- matrix(countTable, nrow = length(names(ctrlmodulecount)), byrow = FALSE)
  
  colnames(logpTable) <- colnames(pTable) <- colnames(countTable) <- names(casemodulecount)
  rownames(logpTable) <- rownames(pTable) <- rownames(countTable) <- names(ctrlmodulecount)
  
  pTable[pTable >= 0.05] <- 1
  
  textMat <- paste(countTable, "\n(", signif(pTable, 3), ")", sep = "")
  textMat <- gsub(pattern = '\n(1)', replacement = '', x = textMat, fixed = TRUE)
  
  par(mfrow = c(1,1))
  par(cex = 1.0)
  par(mar = c(8, 12, 2.7, 1) + 0.3)
  
  labeledHeatmap(Matrix = logpTable, 
                 xLabels = paste(" ", colnames(logpTable)), 
                 yLabels = paste(" ", rownames(logpTable)),
                 colorLabels = TRUE,
                 xSymbols = paste0(casename, colnames(logpTable), ": ", as.vector(casemodulecount)),
                 ySymbols = paste0(ctrlname, rownames(logpTable), ": ", as.vector(ctrlmodulecount)),
                 textMatrix = textMat,
                 colors = greenWhiteRed(100)[50:100],
                 main = paste0('Correspondence of modules (', titlesufix, ')'),
                 cex.text = 1.5, cex.lab = 1.0, setStdMargins = FALSE)
  
  
}


relatetopquantilemodules(ctrlpdcolornames = c('white','red', 'cyan'), 
                         preepdcolornames = c('white','red'), 
                         titlesufix = 'Case to Control')

sub.logpTable = logpTable[2:3,2,drop=F]

rownames(sub.logpTable) = c( "red","cyan")
colnames(sub.logpTable) = c("red")
par(mar = c(12, 12, 2.7, 1))
labeledHeatmap(Matrix = sub.logpTable, 
               xLabels = paste(" ", colnames(sub.logpTable)), 
               yLabels = paste(" ", rownames(sub.logpTable)),
               colorLabels = TRUE,
               xSymbols = paste0(casename, colnames(sub.logpTable), ": ", as.vector(casemodulecount)[1]),
               ySymbols = paste0(ctrlname, rownames(sub.logpTable), ": ", as.vector(ctrlmodulecount)[c(1,2)]),
               #textMatrix = textMat[-c(1,2,3,4)],
               colors = greenWhiteRed(100)[50:100],xLabelsAngle = 30,xColorWidth = .1,
               main = paste0('Correspondence of modules (', titlesufix, ')'),
               cex.text = 1.5, cex.lab = 1.0, setStdMargins = FALSE)

## cor bettern clinical fators and metabolites
#Correlation matrix#
library(ltm)
pd.meta.cor = data.frame(Preterm = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Preterm)}),
                         BMI = apply(meta[-1], 2, function(x) {cor(x, pd$BMI)}),
                         Income = apply(meta[-1], 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Age = apply(meta[-1], 2, function(x) {cor(x, pd$Income)}), 
                         Alcohol =  apply(meta[-1], 2, function(x) {cor(x, pd$Income, method = "spearman")}),
                         Smoker =  apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Smoking)}),
                         SGA = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$SGA)}),
                         LGA = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$LGA)}),
                         BabyLength = apply(meta[-1], 2, function(x) {cor(x, pd$BabyLength)}),
                         Gender = apply(meta[-1], 2, function(x) {biserial.cor(x, pd$Gender)})
                         )

pheatmap(t(pd.meta.cor), 
         color = colorRampPalette(colors = c('green', 'black', 'red'))(100), 
         show_rownames = T, show_colnames = T, 
         scale = 'none', fontsize_row = 15)

d = dist(pd.meta.cor)
test = hclust(d)
memb = cutree(test, k=3)
datanno <- data.frame(group = memb)
datanno$group <- factor(datanno$group, levels = c('1','2','3'), ordered = TRUE)
pdf("by.clinical.vs.meta.corrmat.pdf",width = 8,height = 5)
pheatmap(t(pd.meta.cor), annotation_col = datanno, 
         color = colorRampPalette(colors = c('blue', 'black', 'yellow'))(100), 
         show_colnames = F, cluster_cols = test, cluster_rows = F,
         fontsize = 15,
         scale = 'row', show_rownames = T)
dev.off()
c1 = colnames(meta)[memb==1]
c2 = colnames(meta)[memb==2]
c3 = colnames(meta)[memb==3]
##cluster enrich analysis:

lipidcolor <- sort(cutree(test, k=3))
mode(lipidcolor) <- 'character'

tlipidgroup <- data.frame(Cluster=factor(lipidcolor, labels = paste0('Cluster', seq(3))))
library(scales)

lipidspeciecolor <- data.frame(lipid = c('FA','OLEIC','PC','Cer-AS', 'Cer_NDS', 'Cer_NS',
                                         'DAG','PE', 'PG', 'PS','SM', 'TAG','CAR'), 
                               lipidcolor = c('#F8766D', '#E7851E', '#FF689E','#D09400', '#B2A100','#89AC00', 
                                              '#45B500', '#00BC51', '#00C087', '#00C0B2', 
                                              '#00BCD6', '#00B3F2','#FF61C7'), 
                               stringsAsFactors = FALSE)
lipidgroupenrich <- function(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor){
  
  getlipidspecies <- function(lipidnames = row.names(lipidclusteranno)){
    
    processed <- gsub(pattern = 'Unknown ', replacement = '', x = lipidnames)
    processed <- gsub(pattern = ' .*$', replacement = '', x = processed)
    processed <- gsub(pattern = '\\(.*$', replacement = '', x = processed)
    
    return(processed)
    
  }
  
  lipidclusteranno$species <- getlipidspecies()
  
  enrichres <- function(targets = sub$species, background = lipidclusteranno$species){
    
    targetssum <- length(targets)
    backsum <- length(background)
    
    uniquetargets <- unique(targets)
    for(i in 1:length(uniquetargets)){
      uniquetarget <- uniquetargets[i]
      a11 <- sum(targets == uniquetarget)
      a12 <- sum(background == uniquetarget)
      a21 <- targetssum - a11
      a22 <- backsum - a12
      mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
      fisheres <- fisher.test(mat)
      fisherp <- fisheres$p.value
      ratio <- (a11/(a11 + a21))/(a12/(a12 + a22))
      if(i == 1){
        fisherps <- fisherp
        ratios <- ratio
      }else{
        fisherps <- c(fisherps, fisherp)
        ratios <- c(ratios, ratio)
      }
    }
    
    res <- data.frame(lipid = uniquetargets, fisherp = fisherps, ratio = ratios, stringsAsFactors = FALSE)
    
    return(res)
    
  }
  
  groupnames <- unique(lipidclusteranno$Cluster)
  
  i <- 1
  
  for(i in 1:length(groupnames)){
    
    groupname <- groupnames[i]
    sub <- lipidclusteranno[lipidclusteranno$Cluster == groupname,]
    
    subspecies <- unique(sub$species)
    
    enrichresult <- enrichres(targets = sub$species, background = lipidclusteranno$species)
    enrichresult$cluster <- groupname
    
    if(i == 1){
      enrichresults <- enrichresult
    }else{
      enrichresults <- rbind(enrichresults, enrichresult)
    }
    
  }
  
  enrichresults <- unique(enrichresults)
  
  plotdat <- subset(enrichresults, fisherp < 0.05)
  plotdat <- merge(plotdat, lipidspeciecolor, by = c('lipid'))
  plotdat <- plotdat[order(plotdat$cluster, -plotdat$fisherp),]
  plotdat$logp <- -log2(plotdat$fisherp)
  
  i <- 1
  for(i in 1:nrow(plotdat)){
    line <- plotdat[i,]
    sub <- subset(lipidclusteranno, Cluster == line$cluster)
    a11 <- sum(sub$species == line$lipid)
    a12 <- sum(lipidclusteranno$species == line$lipid)
    a21 <- nrow(sub) - a11
    a22 <- nrow(lipidclusteranno) - a21
    mat <- matrix(c(a11, a12, a21, a22), nrow = 2, byrow = TRUE)
    
    greaterp <- fisher.test(mat, alternative = 'greater')$p.value
    lessp <- fisher.test(mat, alternative = 'less')$p.value
    
    if(greaterp < lessp){
      line$dir <- 'UP'
    }else{
      line$dir <- 'DN'
    }
    
    if(i == 1){
      lines <- line
    }else{
      lines <- rbind(lines, line)
    }
    
  }
  
  plotdat$logp[lines$dir == 'DN'] <- -plotdat$logp[lines$dir == 'DN']
  
  library(ggplot2)
  
  clusternames <- unique(plotdat$cluster)
  
  i <- 1
  for(i in 1:length(clusternames)){
    clustername <- clusternames[i]
    plotsub <- subset(plotdat, cluster == clustername)
    plotsub$lipidname <- factor(plotsub$lipid, levels = plotsub$lipid, ordered = TRUE)
    
    p <- ggplot(plotsub, aes(x = lipidname, y = logp))
    print(
      p + geom_bar(stat = 'identity', fill = plotsub$lipidcolor) + 
        xlab('') + ylab('') + 
        ggtitle(paste0('Significantly enriched lipids in ', clustername)) + 
        theme_bw() + 
        theme(panel.grid = element_blank()) + 
        theme(axis.text.x = element_text(angle = 90, size = 30)) + 
        theme(axis.text.y = element_text(size = 30)) + 
        theme(plot.title = element_text(size = 25)) + 
        theme(panel.border = element_blank(), axis.line = element_line()) + 
        scale_y_continuous(position = 'right') + 
        coord_flip()
    )
    
    
  }
  
  return(plotdat)
  
  
}


lipidenrichres <- lipidgroupenrich(lipidclusteranno = tlipidgroup, lipidspeciecolor = lipidspeciecolor)
FA_inc3 = meta[,c("Label",c3[grep("FA", c3)])]
#plot box plot to check c3 FA is higher in preterm
library(reshape2)
library(ggpubr)
indat = melt(FA_inc3)
ggplot(indat, aes(x = variable, y = value))+ 
  geom_boxplot(aes(fill = Label)) + 
  scale_fill_viridis_d()
p <- ggboxplot(indat, x = "variable", y = "value",
               color = "Label", palette = "jco")
#  Add p-value
p + stat_compare_means() + stat_compare_means(method = "wilcox")+ 
  stat_compare_means( aes(label = ..p.signif..), 
                       label.x = 1, label.y = 10)
p = ggboxplot(indat, x = "Label", y = "value",
          color = "Label", palette = "jco",
               facet.by = "variable", short.panel.labs = FALSE)
# Use only p.format as label. Remove method name.
p + stat_compare_means(label = "p.format")
