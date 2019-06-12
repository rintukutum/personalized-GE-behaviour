rm(list=ls())
load('./data/processed/DC-MTB-12/dc.mtb12.freezed.RData')
load('./data/raw/annotation/DC-MTB-12/dc.mtb12.geneSymbol.RData')
probeIDs <- dc.mtb12.geneSymbol$Probe_Id
dc.temp <- dc.mtb12.freezed
dc.temp$R1$`MTB-infected` <- dc.temp$R1$`MTB-infected`[probeIDs,]
dc.temp$R1$`non-infected` <- dc.temp$R1$`non-infected`[probeIDs,]

dc.temp$R2$`MTB-infected` <- dc.temp$R2$`MTB-infected`[probeIDs,]
dc.temp$R2$`non-infected` <- dc.temp$R2$`non-infected`[probeIDs,]

rm(dc.mtb12.freezed)
probeID.anno <- dc.mtb12.geneSymbol[probeIDs,c('Probe_Id','Symbol','Entrez_Gene_ID')]
rownames(probeID.anno) <- 1:nrow(probeID.anno)
dc.mtb12.working <- list(
  data = dc.temp,
  annotation = probeID.anno
)
save(dc.mtb12.working,
     file = './data/processed/DC-MTB-12/dc.mtb12.working.RData'
)
rm(list=ls())
load('./data/processed/DC-MTB-12/dc.mtb12.working.RData')

R1.inf <- dc.mtb12.working$data$R1$`MTB-infected`
R1.noninf <- dc.mtb12.working$data$R1$`non-infected`

delta <- R1.inf - R1.noninf

responsivity.data <- cbind(R1.noninf, delta)
resV.cor <- apply(
  responsivity.data,
  1,
  function(x){
    x <- as.numeric(x)
    h.x <- length(x)/2
    out <- cor(x[1:h.x],x[(h.x+1):length(x)])
    }
)

response.data <- cbind(R1.inf, delta)
reS.cor <- apply(
  response.data,
  1,
  function(x){
    x <- as.numeric(x)
    h.x <- length(x)/2
    out <- cor(x[1:h.x],x[(h.x+1):length(x)])
  }
)
dir.create('./figures/',showWarnings = FALSE)
dir.create('./figures/DC-MTB-12',showWarnings = FALSE)
png('./figures/DC-MTB-12/01-histogram-DC-MTB-12__R1.png',width = 1200,height = 800,
    res = 150)
par(mfrow=c(1,2))
hist(resV.cor,main = 'Responsivity genes',xlab = 'Correlation')
hist(reS.cor,main = 'Response genes',xlab = 'Correlation')
dev.off()
stats.mtb <- c(table(abs(resV.cor) >= 0.8)['TRUE'],
  table(abs(reS.cor) >= 0.8)['TRUE']
)
names(stats.mtb) <- c('responsivity','response')
png('./figures/DC-MTB-12/01-Statistics-DC-MTB-12__R1.png',width = 1200,height = 800,
    res = 150)
par(mfrow=c(1,1))
barplot(stats.mtb,ylab = '# probes(Genes)',main = 'DC-MTB-2012 | R1')
dev.off()

#####
which.max(abs(resV.cor)) # LAMP3
which.max(abs(reS.cor)) # CXCR7

#####
extractCorAnno <- function(cor,threshold=0.8){
   probes.cor <- cor[abs(cor) >= threshold]
   idx <- dc.mtb12.working$annotation$Probe_Id %in% names(probes.cor)
   genes.cor <- dc.mtb12.working$annotation[idx,]
   genes.cor$cor <- as.numeric(probes.cor[genes.cor$Probe_Id])
   return(genes.cor)
}

R1.responsivity.genes <- extractCorAnno(
   cor = resV.cor
)
R1.response.genes <- extractCorAnno(
   cor = reS.cor
)

#####
R2.inf <- dc.mtb12.working$data$R2$`MTB-infected`
R2.noninf <- dc.mtb12.working$data$R2$`non-infected`

R2.delta <- R2.inf - R2.noninf

R2.responsivity.data <- cbind(R2.noninf, R2.delta)
R2.resV.cor <- apply(
  R2.responsivity.data,
  1,
  function(x){
    x <- as.numeric(x)
    h.x <- length(x)/2
    out <- cor(x[1:h.x],x[(h.x+1):length(x)])
  }
)

R2.response.data <- cbind(R2.inf, R2.delta)
R2.reS.cor <- apply(
  R2.response.data,
  1,
  function(x){
    x <- as.numeric(x)
    h.x <- length(x)/2
    out <- cor(x[1:h.x],x[(h.x+1):length(x)])
  }
)


png('./figures/DC-MTB-12/01-histogram-DC-MTB-12__R2.png',width = 1200,height = 800,
    res = 150)
par(mfrow=c(1,2))
hist(R2.resV.cor,main = 'Responsivity genes',xlab = 'Correlation')
hist(R2.reS.cor,main = 'Response genes',xlab = 'Correlation')
dev.off()
R2.stats.mtb <- c(table(abs(R2.resV.cor) >= 0.8)['TRUE'],
               table(abs(R2.reS.cor) >= 0.8)['TRUE']
)
names(stats.mtb) <- c('responsivity','response')
png('./figures/DC-MTB-12/01-Statistics-DC-MTB-12__R2.png',width = 1200,height = 800,
    res = 150)
par(mfrow=c(1,1))
barplot(R2.stats.mtb,ylab = '# probes(Genes)',main = 'DC-MTB-2012 | R2')
dev.off()

which.max(abs(R2.resV.cor)) # LAMP3 # ILMN_2170814
which.max(abs(R2.reS.cor)) # CXCR7 # ILMN_1798360

R2.responsivity.genes <- extractCorAnno(
   cor = R2.resV.cor
)
R2.response.genes <- extractCorAnno(
   cor = R2.reS.cor
)
dc.mtb.12.profile <- list(
   responsivity = rbind(
     data.frame(
         R1.responsivity.genes,
         replicate = 'R1',
         stringsAsFactors = FALSE
     ),
     data.frame(
         R2.responsivity.genes,
         replicate = 'R2',
         stringsAsFactors = FALSE
     )
     ),
     response = rbind(
     data.frame(
         R1.response.genes,
         replicate = 'R1',
         stringsAsFactors = FALSE
     ),
     data.frame(
         R2.response.genes,
         replicate = 'R2',
         stringsAsFactors = FALSE
     )
     )
)
dir.create('./data/results/',showWarnings=FALSE)
dir.create('./data/results/DC-MTB-12/',showWarnings=FALSE)


save(
  dc.mtb.12.profile,
  file = './data/results/DC-MTB-12/dc.mtb.12.profile.RData'
)

png('./figures/DC-MTB-12/01-Replicate-results-LAMP3.png',
    width = 1200,height = 800,res = 150)
par(mfrow=c(1,2))
plot(y = delta['ILMN_2170814',],
     x = R1.noninf['ILMN_2170814',],
     ylab='delta',
     xlab = 'baseline',
     main = 'R1 | ILMN_2170814, LAMP3')

plot(y = R2.delta['ILMN_2170814',],
     x = R2.noninf['ILMN_2170814',],
     ylab='delta',
     xlab = 'baseline',
     main = 'R2 | ILMN_2170814, LAMP3')

dev.off()
#-------------------
#     Gene ontology
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

if (!requireNamespace("clusterProfiler", quietly = TRUE))
  BiocManager::install("clusterProfiler")

library("clusterProfiler")
filter.kegg <- function(kegg.output){
  idx <- kegg.output@result$pvalue <= 0.05
  kegg.output@result <- kegg.output@result[idx,]
  return(kegg.output)
}

load('./data/results/DC-MTB-12/dc.mtb.12.profile.RData')

responsivity.kegg <- filter.kegg(enrichKEGG(
  gene = unique(dc.mtb.12.profile$responsivity$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
))

response.kegg <- filter.kegg(enrichKEGG(
  gene = unique(dc.mtb.12.profile$response$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
))

# 121 entrez id are common
intersect(unique(dc.mtb.12.profile$response$Entrez_Gene_ID),
          unique(dc.mtb.12.profile$responsivity$Entrez_Gene_ID))

R.resV <- plyr::dlply(
  .data = dc.mtb.12.profile$responsivity,
  .variables = 'replicate')

R1.resV.kegg <- filter.kegg(enrichKEGG(
  gene = unique(R.resV$R1$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
))
#----------
# 16 pathways
dim(R1.resV.kegg@result)

R2.resV.kegg <- filter.kegg(enrichKEGG(
  gene = unique(R.resV$R2$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
  )
)
#----------
# 40 pathways
dim(R2.resV.kegg@result)
#----------
# 7 KEGG pathways common
intersect(
  R1.resV.kegg@result$Description,
  R2.resV.kegg@result$Description
)

setdiff(
  R1.resV.kegg@result$Description,
  R2.resV.kegg@result$Description
)

setdiff(
  R2.resV.kegg@result$Description,
  R1.resV.kegg@result$Description
)

R.res <- plyr::dlply(
  .data = dc.mtb.12.profile$response,
  .variables = 'replicate')

R1.res.kegg <- filter.kegg(enrichKEGG(
  gene = unique(R.res$R1$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
))
#----------
# 72 pathways
dim(R1.res.kegg@result)

R2.res.kegg <- filter.kegg(clusterProfiler::enrichKEGG(
  gene = unique(R.res$R2$Entrez_Gene_ID),
  organism = 'hsa',keyType = "ncbi-geneid"
))
#----------
# 84 pathways
dim(R2.res.kegg@result)
#----------
# 57 KEGG pathways common
intersect(
  R1.res.kegg@result$Description,
  R2.res.kegg@result$Description
)


responsivity.GObp <- enrichGO(
  gene = unique(dc.mtb.12.profile$responsivity$Entrez_Gene_ID),
  OrgDb = 'org.Hs.eg.db',
  ont = 'BP'
)
de_symbols <- read.table('clipboard',stringsAsFactors = FALSE)

intersect(
  unique(dc.mtb.12.profile$responsivity$Symbol),
  unique(de_symbols$V1)
)

intersect(
  unique(dc.mtb.12.profile$response$Symbol),
  unique(de_symbols$V1)
)

setdiff(
  c(
    unique(dc.mtb.12.profile$responsivity$Symbol),
    unique(dc.mtb.12.profile$response$Symbol)
  ),
  unique(de_symbols$V1)
)
