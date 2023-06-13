library( ChIPseeker )
library( Signac )
library( Seurat )
library(TxDb.Mmusculus.UCSC.mm10.knownGene)
txdb <- TxDb.Mmusculus.UCSC.mm10.knownGene

obj = readRDS('NAc.Cleaned.FinalObj.rds')
DefaultAssay(obj) = 'macs2_peaks'

peaks = granges(obj)

peakAnno <- annotatePeak(peaks, tssRegion=c(-1000, 100),
                         TxDb=txdb, annoDb="org.Mm.eg.db")

df = as.data.frame( peakAnno )

df$peak = paste( df$seqnames , df$start , df$end , sep = '-' )
 
df = merge( rownames(obj) , df , by.x = 1 , by.y = 'peak' , all.x = T )
colnames(df)[1] = 'peak'
rownames(df) = df$peak
df = df[ rownames(obj) , ]

write.table(df , quote=F , sep='\t' , row.names=F , file = 'NAc.AllMACS2Peaks.ChIPseeker.txt')
