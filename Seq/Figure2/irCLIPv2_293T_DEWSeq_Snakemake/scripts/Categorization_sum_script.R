library(DEWSeq)
library(IHW)
library(tidyverse)
library(data.table)

cat("Load Data\n")
w <- read.delim(snakemake@input$window, header = TRUE)
r <- read.delim(snakemake@input$region, header = TRUE)
ddw <- readRDS(snakemake@input$ddw)
gene <- sapply(strsplit(rownames(colData(ddw))[1], "_"), function(x) x[1])
setwd(file.path(snakemake@params$outdir))

cat("Prepare datasets\n")
r$regionID <- paste(r$chromosome, r$region_begin, r$region_end, r$strand, r$regionStartId, sep = "_")
r2 <- r %>% mutate(unique_id = str_split(unique_ids, ", ")) %>% unnest(unique_id)
countsn <- as.data.frame(counts(ddw, normalized=TRUE))
countsn$unique_id <- rownames(countsn)
rownames(countsn) <- NULL
r2 <- merge(countsn, r2, by = "unique_id")

cat("Categorize region according to slope value\n")
r3 <- r2 %>% group_by(regionID) %>% summarise(across(paste(gene, "_293T_UVC_R1_S1", sep = ""):paste(gene, "_293T_noUV_R2", sep = ""), ~ sum(.x)))
r3 <- as.data.frame(r3)
rownames(r3) <- r3$regionID
r3$regionID <- NULL
countsn <- r3[,c(1,4,2,5,3,6)]
mean <- data.frame(apply(array(as.matrix(countsn), c(nrow(countsn),2, ncol(countsn)/2)),3, rowMeans))
rownames(mean) <- rownames(countsn)
column <- c("avg_S1", "avg_S2", "avg_S3")
colnames(mean) <- column
mean <- log10(mean)
mean[mean == -Inf] <- 0

slope  <-  function(x){
  if(all(is.na(x)))
    return(NA)
  else
    return(coef(lm(x~I(1:3)))[2])
}

mean$lmslope <- apply(mean, 1,slope)
mean$regionID <- rownames(mean)
rownames(mean) <- NULL

mean.l <- subset(mean, lmslope < -sd(mean$lmslope, na.rm = TRUE))
mean.l$category <- "low"
mean.m <- subset(mean, lmslope > -sd(mean$lmslope, na.rm = TRUE) & lmslope < sd(mean$lmslope, na.rm = TRUE))
mean.m$category <- "med"
mean.h <- subset(mean, lmslope > sd(mean$lmslope, na.rm = TRUE))
mean.h$category <- "high"
cat <-  rbind(mean.l, mean.m, mean.h)

cat("Determine windows with highest fold change inside regions\n")
section <- levels(ddw$section)

get_results <- function(dew, section1, section2) {
  resultWindows <- resultsDEWSeq(dew, contrast = c("section", section1, section2), tidy = TRUE) %>% as_tibble
  resultWindows[,"padj"] <- adj_pvalues(ihw(pSlidingWindows ~ baseMean, data = resultWindows, alpha = 0.05, nfolds = 10))
  resultWindows <- resultWindows %>% mutate(significant = resultWindows$padj < 0.05)
  return(resultWindows)
}

cat("Windows: Subzone vs input\n")
S1_I <- get_results(ddw, section[2], section[1])
S1_I$section <- "S1_I"
S2_I <- get_results(ddw, section[3], section[1])
S2_I$section <- "S2_I"
S3_I <- get_results(ddw, section[4], section[1])
S3_I$section <- "S3_I"

r4 <- r2
r4 <- merge(r4, S3_I[,c(6,15)], by = "unique_id", all.x = TRUE)
r4 <- r4 %>% dplyr::rename( logFC_S3 = log2FoldChange )
r4$logFC_S3[is.na(r4$logFC_S3)] <- 0
r4 <- merge(r4, S2_I[,c(6,15)], by = "unique_id", all.x = TRUE)
r4 <- r4 %>% dplyr::rename( logFC_S2 = log2FoldChange )
r4$logFC_S2[is.na(r4$logFC_S2)] <- 0
r4 <- merge(r4, S1_I[,c(6,15)], by = "unique_id", all.x = TRUE)
r4 <- r4 %>% dplyr::rename( logFC_S1 = log2FoldChange )
r4$logFC_S1[is.na(r4$logFC_S1)] <- 0
r4 <- r4 %>% group_by(regionID) %>% mutate(max_S1_id = unique_id[which.max(logFC_S1)], max_S1_fc = logFC_S1[which.max(logFC_S1)],
                                           max_S2_id = unique_id[which.max(logFC_S2)], max_S2_fc = logFC_S2[which.max(logFC_S2)],
                                           max_S3_id = unique_id[which.max(logFC_S3)], max_S3_fc = logFC_S3[which.max(logFC_S3)]) %>% ungroup()
max <- unique(r4 %>% dplyr::select(regionID, max_S1_id:max_S3_fc))

cat("Create a table with slope categorization and max window in each region for each slice\n")
r <- merge(r %>% dplyr::select(chromosome:region_length, gene_id:gene_region,unique_ids,regionID), cat, by = "regionID")
r <- merge(r, max, by = "regionID")

w <- w %>% dplyr::select(chromosome:end, strand:unique_id)
w <- unique(w)

r_f <- left_join(r, w, by = c("max_S1_id" = "unique_id"))
r_f <- r_f %>% dplyr::rename(all_of(c(chr_S1 = "chromosome.y", begin_S1 = "begin", end_S1 = "end", strand_S1 = "strand.y" )))
r_f <- left_join(r_f, w, by = c("max_S2_id" = "unique_id"))
r_f <- r_f %>% dplyr::rename(all_of(c(chr_S2 = "chromosome", begin_S2 = "begin", end_S2 = "end", strand_S2 = "strand" )))
r_f <- left_join(r_f, w, by = c("max_S3_id" = "unique_id"))
r_f <- r_f %>% dplyr::rename(all_of(c(chromosome = "chromosome.x", strand = "strand.x", chr_S3 = "chromosome", begin_S3 = "begin", end_S3 = "end", strand_S3 = "strand" )))

cat("Saving the results\n")
write.table(r_f, file = paste(getwd(), paste(gene, "DEWSeq_sign_regions_sum_categorization.txt", sep = "_"), sep="/"), sep = "\t", row.names = FALSE, quote = FALSE)






