library("DEXSeq")
library("stringi")
library("tidyverse")
library('BiocParallel')
library('parallel')

BPPARAM = MulticoreParam(workers=6)

args = commandArgs(trailingOnly=TRUE)

# Optparse wouldn't install in the conda environment.
COUNTS = args[1]  # subread_fc.out counts file.
SAMPLE_TABLE = args[2]  # column of sample names match bam files and column of groups.
GTF = args[3]  # Can't just be a GTF. Has to be reformatted.
GROUP_A = args[4]  # Should match the "group" values in SAMPLE_TABLE.
GROUP_B = args[5]

DEX_RESULTS = paste0(dirname(COUNTS), "/dexseq")
dir.create(DEX_RESULTS, showWarnings = F, recursive = T)

SIMPLIFIED_COUNTS = paste0(COUNTS, '.simplified_col_names')

if (length(args)<5) {
  stop("At least 5 arguments must be supplied.", call.=FALSE)
}

tables_from_dexseq = function(results_folder){
    # No, I'm not very good with R, that's true.
    
    res = readRDS(paste0(DEX_RESULTS, "/dxd_results.R"))

    # Make data.frame to print results.
    q = res@listData

    fc_col = names(q)[grepl(pattern =  "log2", x = names(q))]
    to_include = strsplit(fc_col, "_")[[1]][2:3]
    
    a = names(q)[grepl(pattern =  GROUP_A, x = names(q))]
    a_col = a[!grepl(pattern = 'log2', a)]
    b = names(q)[grepl(pattern =  GROUP_B, x = names(q))]
    b_col = b[!grepl(pattern = 'log2', b)]
    
    # There is some correct way to do this, which is not this.
    df_to_write = cbind(q$groupID, q$featureID, q$dispersion, q$exonBaseMean, 
                        as.character(q$transcripts), as.character(q$genomicData@seqnames), 
                        as.character(q$genomicData@ranges@start), 
                        as.character(q$genomicData@ranges@start + q$genomicData@ranges@width),
                        as.character(q$genomicData@strand),
                        unlist(q[a_col]), unlist(q[b_col]), unlist(q[fc_col]), q$padj)
    df_to_write = data.frame(df_to_write)

    colnames(df_to_write) = c("groupID", "featureID", "dispersion", "exonBaseMean",
                              "transcripts", "chrm", 
                              "start",
                              "end",
                              "strand",
                              GROUP_A, GROUP_B, "log2fc", "padj")

    # If we truncated ENSG (removed the "ENSG" and leading zeroes), recreate the full ID.
    #recover_ensg = function(s) {ifelse(
    #    grepl(pattern = "ENSG", x = s), paste(s), paste("ENSG", str_pad(s, 11, pad="0"), sep="") )}

    #df_to_write %>% mutate(ENSG=recover_ensg(groupID))

    # Print results as a table.
    write.table(df_to_write, file=paste0(DEX_RESULTS, "/DEXSeq_results.txt"), sep='\t', quote = F, row.names = F)

    return(0)
}

source("Subread_to_DEXSeq/load_SubreadOutput.R")

# Read subread counts. check.names=F to prevent converting '/' and '_' to '.'.
fc = read.csv(COUNTS, sep='\t', comment.char = '#', check.names=F)

# Simplify column names to the basename without the STAR _Aligned suffix.
bamName = function(x) {ifelse(grepl('bam', x), gsub("(.*)_Aligned.out.bam", "\\1", basename(x)), x) }
colnames(fc) = lapply(colnames(fc), bamName)

# Write the file with simplified column names.
con <- file(SIMPLIFIED_COUNTS, 'w')
writeLines("# Header line required for DEXSeqDataSetFromFeatureCounts.R", con = con)
write.table(fc, file = con, sep = '\t', quote = FALSE, row.names=FALSE)
close(con)

# Function to perform DEXSeq given counts and gtf filenames, plus a dataframe of conditions.
call_dexseq = function(subread_fname, gtf_fname, sampleData) {

    # Create output directory.
    results_folder = DEX_RESULTS #paste0(top, "/outs/R/dexseq/", comparison_name, "/")

    #> formuladispersion <- count ~ sample + (exon + type) * condition
    #> pasillaExons <- estimateDispersions(pasillaExons, formula = formuladispersion)
    
    sampleData$batch <- as.factor(sampleData$batch)
    # Call DEXSeq.
    dxd <- DEXSeqDataSetFromFeatureCounts(
        subread_fname, flattenedfile = gtf_fname, sampleData = sampleData,
        design = ~ sample + (exon + batch) * condition)
    # default: design = ~sample + exon + condition:exon,

    saveRDS(dxd, file=paste0(results_folder, "/dxd_object.R"))

    #sh_samples = paste0('W', 1:6)
    #fc_out$counts = fc_out$counts[,sh_samples]
    #dGeneId = fc_out$annotation$GeneID
    #dFeatureID = paste(fc_out$annotation$Start, fc_out$annotation$End, sep='.')
    #design <- formula( ~ sample + exon + condition:exon )
    #dxd = DEXSeqDataSet(fc_out$counts, sampleData, design, dFeatureID, dGeneId, normalized=F)
    #dxd = DEXSeqDataSet(countData, sampleData, design, featureID, groupID, normalized=FALSE)

    formuladispersion <- ~ sample + (exon + batch) * condition
    # The next four function calls each take hours.
    dxd = estimateSizeFactors(dxd)
    dxd = estimateDispersions(dxd, formula = formuladispersion, BPPARAM=BPPARAM)
    
    # The formula0, formula1 arguments to testForDEU are not supported in v. 1.32.
    # Null hypothesis formula.
    formula0 <- ~ sample + batch * exon + condition
    # Formula where the exon identity matters (rejecting the null hypothesis).
    formula1 <- ~ sample + batch * exon + condition * I(exon == exonID)
    
    print(paste0("design(dxd)=", design(dxd)))
    
    # For v. 1.32:
    #defaults for testForDEU(): fullModel = design(object),  reducedModel = ~ sample + exon, 
    # Should reducedModel be specified as ~ sample + batch * exon ?
    dxd = testForDEU(dxd, BPPARAM=BPPARAM)
    
    #dxd = testForDEU(dxd, fullModel = formula1, reducedModel = formula0, BPPARAM=BPPARAM)
    dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

    saveRDS(dxd, file=paste0(results_folder, "/dxd_object.R"))

    # Get results object.
    res = DEXSeqResults(dxd)

    saveRDS(res, file=paste0(results_folder, "/dxd_results.R"))

    tables_from_dexseq(results_folder)

    return(0)
}

# Function to perform DEXSeq given counts and gtf filenames, plus a dataframe of conditions.
call_dexseq_no_batch = function(subread_fname, gtf_fname, sampleData) {

    # Create output directory.
    results_folder = DEX_RESULTS #paste0(top, "/outs/R/dexseq/", comparison_name, "/")

    #> formuladispersion <- count ~ sample + (exon + type) * condition
    #> pasillaExons <- estimateDispersions(pasillaExons, formula = formuladispersion)
    
    # Call DEXSeq.
    dxd <- DEXSeqDataSetFromFeatureCounts(
        subread_fname, flattenedfile = gtf_fname, sampleData = sampleData,
        design = ~ sample + exon * condition)
    # default: design = ~sample + exon + condition:exon,

    saveRDS(dxd, file=paste0(results_folder, "/dxd_object.R"))

    #sh_samples = paste0('W', 1:6)
    #fc_out$counts = fc_out$counts[,sh_samples]
    #dGeneId = fc_out$annotation$GeneID
    #dFeatureID = paste(fc_out$annotation$Start, fc_out$annotation$End, sep='.')
    #design <- formula( ~ sample + exon + condition:exon )
    #dxd = DEXSeqDataSet(fc_out$counts, sampleData, design, dFeatureID, dGeneId, normalized=F)
    #dxd = DEXSeqDataSet(countData, sampleData, design, featureID, groupID, normalized=FALSE)

    formuladispersion <- ~ sample + exon * condition
    # The next four function calls each take hours.
    dxd = estimateSizeFactors(dxd)
    dxd = estimateDispersions(dxd, formula = formuladispersion, BPPARAM=BPPARAM)
    
    # The formula0, formula1 arguments to testForDEU are not supported in v. 1.32.
    # Null hypothesis formula.
    formula0 <- ~ sample + exon + condition
    # Formula where the exon identity matters (rejecting the null hypothesis).
    formula1 <- ~ sample + exon + condition * I(exon == exonID)
    
    print(paste0("design(dxd)=", design(dxd)))
    
    # For v. 1.32:
    #defaults for testForDEU(): fullModel = design(object),  reducedModel = ~ sample + exon, 
    # Should reducedModel be specified as ~ sample + batch * exon ?
    dxd = testForDEU(dxd, BPPARAM=BPPARAM)
    
    #dxd = testForDEU(dxd, fullModel = formula1, reducedModel = formula0, BPPARAM=BPPARAM)
    dxd = estimateExonFoldChanges(dxd, fitExpToVar="condition", BPPARAM=BPPARAM)

    saveRDS(dxd, file=paste0(results_folder, "/dxd_object.R"))

    # Get results object.
    res = DEXSeqResults(dxd)

    saveRDS(res, file=paste0(results_folder, "/dxd_results.R"))

    tables_from_dexseq(results_folder)

    return(0)
}

read_sampleTable = function() {
    # Fix this mess.
    
    # Make the sampleTable data.frame to pass to DESeq.
    samples = read.csv(SAMPLE_TABLE, sep='\t')
    rownames(samples) = samples$sample
#    samples$condition = samples$group
    
    batch_from_name = function(x) { str_extract(x, "(Batch_\\d+)") }
    samples$batch = factor(unlist(lapply(samples$sample, batch_from_name)))
    
    sampleTable = data.frame(
        condition = samples$group, batch = factor(samples$batch),
        sample=samples$sample
    )
    rownames(sampleTable) = rownames(samples)
    return(sampleTable)
}

####################################################################################
# Run the comparison.

df = read_sampleTable()
sampleData = df %>% filter(condition==GROUP_A | condition==GROUP_B)

# Write a smaller counts file with just the shRNA groups.
subread_fname = paste0(SIMPLIFIED_COUNTS, "_subset_to_", GROUP_A, "_vs_", GROUP_B)
con <- file(subread_fname, 'w')
writeLines("# Header line required for DEXSeqDataSetFromFeatureCounts.R", con = con)
write.table(
    fc[,c("Geneid", "Chr", "Start", "End", "Strand", "Length", row.names(sampleData))],
    file = con, sep = '\t', quote = FALSE, row.names=FALSE)
close(con)

# Call the function to run DEXSeq.
ifelse(
    length(unique(sampleData$batch))==1,
    call_dexseq_no_batch(subread_fname, GTF, sampleData),
    call_dexseq(subread_fname, GTF, sampleData)
)
