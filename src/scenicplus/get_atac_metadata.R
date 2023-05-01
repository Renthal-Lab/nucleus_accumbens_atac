# Get metadata of nucleus accumbens ATAC-seq data

atac_rds = "/path/to/NAc.snATAC.ReMACS2.Motifs.ChromVAR.rds"
output_dir = '/path/to/output_dir/'
  
atac_rds_data <- readRDS(atac_rds)

write.csv(atac_rds_data@meta.data, paste0(output_dir, "data/nucleus_accumbens_metadata.csv"))
