library(DEWSeq)
library(IHW)
library(tidyverse)
library(data.table)

cat("Load Data\n")
countData <- fread(snakemake@input[['matrix']], sep = "\t")
annotationData <- fread(snakemake@input[['annotation']], sep = "\t")
gene <- data.frame(geneid = factor(unique(sapply(strsplit(colnames(countData[,-1]), "_"), function(x) x[1]))))
CUTOFF <- snakemake@params[['padj']]
CUTOFF2 <- snakemake@params[['logfc']]
setwd(file.path(snakemake@params[['outdir']]))

cat("Calculate sum\n")
countData2 <- countData %>% mutate(!!paste(gene[1,], "_293T_UVC_R1", sep = "") := countData %>% dplyr::select(paste(gene[1,], "_293T_UVC_R1_S1", sep = ""):paste(gene[1,], "_293T_UVC_R1_S3", sep = "")) %>% rowSums, 
                                          !!paste(gene[1,], "_293T_UVC_R2", sep = "") := countData %>%  dplyr::select(paste(gene[1,], "_293T_UVC_R2_S1", sep = ""):paste(gene[1,], "_293T_UVC_R2_S3", sep = "")) %>% rowSums)

countData2 <- countData2 %>% dplyr::select(unique_id, paste(gene[1,], "_293T_UVC_R1", sep = ""),paste(gene[1,], "_293T_UVC_R2", sep = ""),
                                           paste(gene[1,], "_293T_noUV_R1", sep = ""),paste(gene[1,], "_293T_noUV_R2", sep = ""))

cat("Prepare colData\n")
count <- countData2
colData <- data.frame(row.names = colnames(count[,-1]), 
                      type = factor(sapply(strsplit(colnames(count[,-1]), "_"), function(x) x[3])),
                      cell = factor(sapply(strsplit(colnames(count[,-1]), "_"), function(x) x[2])),
                      rep = factor(sapply(strsplit(colnames(count[,-1]), "_"), function(x) x[4])))
colData$type <- factor(colData$type)

cat("Creating DEWSeq object\n")
ddw <- DESeqDataSetFromSlidingWindows(countData  = count,
                                      colData    = colData,
                                      annotObj   = annotationData,
                                      tidy       = TRUE,
                                      design     = ~type)

cat("Apply sizefactor and filtering\n")
ddw <- estimateSizeFactors(ddw)
keep <- rowSums(counts(ddw)) >= 10
ddw <- ddw[keep,]

cat("Performing DE analysis\n")
ddw <- estimateDispersions(ddw, fitType = "local", quiet = TRUE)
ddw <- nbinomWaldTest(ddw)

cat("Saving RDS objects\n")
saveRDS(ddw, file = paste(getwd(), paste(gene[1,], "DEWseq_res_sum.rds", sep = "_"), sep = "/"))

cat("Getting the results\n")
section <- levels(colData$type)

get_results <- function(dew, section1, section2) {
  resultWindows <- resultsDEWSeq(dew, contrast = c("type", section1, section2), tidy = TRUE) %>% as_tibble
  resultWindows[,"padj"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, data = resultWindows, alpha = 0.05, nfolds = 10))
  resultWindows <- resultWindows %>% mutate(significant = resultWindows$padj < CUTOFF)
  sign.windows <- resultWindows %>% filter(significant) %>% filter(log2FoldChange > CUTOFF2) %>% arrange(desc(log2FoldChange))
  return(sign.windows)
}

cat("Windows: Subzone vs input\n")
sign.wdw.final.i <- get_results(ddw, section[2], section[1])

cat("Regions: collapse unique regions between S1, S2,and S3\n")
sign.region.final.i <- extractRegions(windowRes  = sign.wdw.final.i, padjCol = "padj", padjThresh = CUTOFF, log2FoldChangeThresh = CUTOFF2) %>% as_tibble

# Write to textfile
cat("Saving the results\n")
write.table(sign.wdw.final.i, file = paste(getwd(), paste(gene[1,], "DEWSeq_sign_windows_sumvsinput.txt", sep = "_"), sep="/"), sep = "\t", row.names = FALSE, quote = FALSE)
write.table(sign.region.final.i, file = paste(getwd(), paste(gene[1,], "DEWSeq_sign_regions_sumvsinput.txt", sep = "_"), sep="/"), sep = "\t", row.names = FALSE, quote = FALSE)



