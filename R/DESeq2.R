# Create heatmaps given a geneset list (adapted from DESeq2.Rmd)
setwd("/Users/kaihu/Projects/Analytical_projects/Lin_Zhou_RNA-Seq/OneStopRNAseq/temp_res")

library(BiocManager)
library(DESeq2)
library(EnhancedVolcano)
library(ggplot2)
library(ggpubr)
library(ggrepel)
library(ashr)
library(MASS)
library(WriteXLS)
library(plyr)
library(gdata)
library(dplyr)
library(RColorBrewer)
library(pheatmap)
library(PoiClaClu)  
library(gridExtra)  
library(grid)
library(berryFunctions)
library(readr)
library(textshape)

# Rscript -e "rmarkdown::render(         './DESeq2/DESeq2.Rmd',         params=list(             max_fdr=0.05,             min_lfc=1,             cookscutoff=TRUE,             indfilter=FALSE,             countFile='../feature_count/counts.strict.txt',             annoFile='https://raw.githubusercontent.com/radio1988/OneStopRNAseq/master/snakemake/envs/anno_tables/mm10.gencode.vm25.te_included.anno.txt',             metaFile='../meta/meta.xlsx',             contrastFile='../meta/contrast.de.xlsx',             blackSamples='sampleXXXXX,sampleYYY'             ),         output_file='DESeq2.html')" > DESeq2/DESeq2.log 2>&1 ;D=DESeq2; rm -f $D/$D.zip && [ -d $D ] && zip -rq  $D/$D.zip $D/ &>> DESeq2/DESeq2.log;

max_fdr <- 0.05
min_lfc <- 1
indfilter <- FALSE
cookscutoff <- TRUE
countFile <- "feature_count/counts.strict.txt"
# annoFile <- "https://raw.githubusercontent.com/radio1988/OneStopRNAseq/master/snakemake/envs/anno_tables/mm10.gencode.vm25.te_included.anno.txt"
annoFile <- "https://raw.githubusercontent.com/radio1988/OneStopRNAseq/master/snakemake/workflow/resources/anno_tables/mm10.gencode.vm25.te_included.anno.txt"
metaFile <- "meta/meta.xlsx"
contrastFile <- "meta/contrast.de.xlsx"
blackSamples <- 'sampleXXXXX,sampleYYY'

## Functions used:
source("scripts/utils.R")
source("scripts/plots.R")

## Read in data: from feature count output
meta.df <- readExcel(metaFile)
meta.df <- meta.df[!(meta.df[, 1] %in% blackSamples), ]
contrast.df <- readExcel(contrastFile)
df <- read_delim(countFile, delim = "\t", comment = '#') 

df <- df[, !(colnames(df) %in% blackSamples)]
df <- column_to_rownames(df, 'Geneid')
colnames(df) <- gsub("\\.bam$", "", colnames(df))
colnames(df) <- gsub("sorted_reads.", "", colnames(df))
colnames(df) <- gsub("mapped_reads.", "", colnames(df))
df <- df %>% mutate_at(vars(6:ncol(df)), as.integer)
cts <- as.matrix(df[, 6:ncol(df)])
cts <- cts[, match(meta.df$SAMPLE_LABEL, colnames(cts))] # match order in meta

## Read in annotation file:
anno <- getAnnotation(annoFile)
anno <- anno[!duplicated(anno[[1]]), ] # assuming first column as ENSEMBL ID
anno$Name <- toupper(anno$Name)

## Filter out lowly expressed genes by rowMeans
cts <- cts[rowMeans(cts) >= 2, ] # reduce rows from 24621 to 18564 if usig rowMeans >= 2 vs rowSums(cts) >= 10
df <- df[match(row.names(cts), row.names(df)), ]
group_count <- plyr::count(meta.df, 'GROUP_LABEL')

match_order <- match(colnames(cts), meta.df$SAMPLE_LABEL)
match_order <- match_order[!is.na(match_order)]
meta.df <- meta.df[match_order, ] # match row order in meta.df and remove those that don't exist in colnames(cts)

# match meta.df row order
tpm <- calculateTPM(cts, df$Length) %>% data.frame()

## Design matrix:
sample <- factor(meta.df$SAMPLE_LABEL)
batch <- factor(meta.df$BATCH)
group <- factor(meta.df$GROUP_LABEL)
coldata <- data.frame(row.names = colnames(cts), sample, batch, group)

## Model fitting:
if (length(levels(factor(meta.df$BATCH))) > 1) {
	dds <- DESeqDataSetFromMatrix(countData = cts, 
																colData = coldata, 
																design = ~  0 + group + batch) 
} else {
	dds <- DESeqDataSetFromMatrix(countData = cts, 
																colData = coldata, 
																design = ~  0 + group)
}
dds <- DESeq(dds)

## QC Plots
### Dispersion plot:
plotDispEsts(dds, main = "Dispersion plot")

### Sample PCA plot:
plotQC_PCA(dds, ntop = 5000)

## DEG analysis:
# i <- 1
for (i in 1:ncol(contrast.df)) {
	name1 <- parse_name(contrast.df[1, i]) # would be treatment
	name2 <- parse_name(contrast.df[2, i]) # would be control
	poss <- match(name1$names, gsub("group", "", resultsNames(dds)))
	negs <- match(name2$names, gsub("group", "", resultsNames(dds)))
	contrast <- rep(0, length(resultsNames(dds)))
	contrast[poss] <- 1/length(poss)
	contrast[negs] <- -1/length(negs)
	# double check the contrast:
	print(data.frame(resNames=gsub("group", "", resultsNames(dds)), 
									 contrast=contrast))
	
	res <- lfcShrink(dds, contrast = contrast, type = 'ashr')
	res2 <- results(dds, contrast = contrast, independentFilter = indfilter, cooksCutoff = cookscutoff)
	norm_exp = tpm
	
	res.df <- as.data.frame(res)
	names(res.df)[2] <- "log2FoldChange_shrunken"
	names(res.df)[3] <- "lfcSE_shrunken"
	
	res2.df <- as.data.frame(res2)
	names(res2.df)[2] <- "log2FoldChange_raw"
	names(res2.df)[3] <- "lfcSE_raw"
	res2.df <- res2.df[, c(2,3)]
	
	resdata <- merge(res.df, res2.df, by=0, sort=F, all.x=T)
	resdata <- merge(anno, resdata, by.x=1, by.y=1, sort=F, all.y=T)
	resdata <- merge(resdata, norm_exp, by.x=1, by.y=0, all.x=T, sort=F)
	name <- paste(name1$name, name2$name, sep = "_vs_")
	output_excel <- paste0("DEG_list/", name, "_deseq2.xlsx")
	WriteXLS(x = resdata,
					 ExcelFileName = output_excel,
					 row.names = F, SheetNames = 'sheet1', na = 'NA')
	
	## For GSEA analysis:
	rnk <- subset(resdata, select = c("Name","log2FoldChange_shrunken"))
	colnames(rnk) <- c("# Name","log2FoldChange_shrunken")
	rnk <- rnk[order(rnk$log2FoldChange_shrunken), ]
	output_rnk <- paste0("RNK_files/", name, ".rnk")
	write.table(rnk, output_rnk, row.names = F, quote = F, sep = "\t")
	
	## Plot: P-value distribution
	ggplot(data.frame(res), aes(x=pvalue)) +
		geom_histogram(color="darkblue", fill="lightblue")
	output_pvalue <- paste0("Pvalue_distribution/", name, ".pdf")
	ggsave(output_pvalue)
	
	## Plot: MA-plot and Volcano plot
	output_ma <- paste0("MA_plots/", name, ".shrunken_lfc.pdf")
	maplot(res, anno, output_ma)
	output_volcano <- paste0("Volcano_plots/", name, ".shrunken_lfc.pdf")
	volcanoplot(res, anno, lfcthresh = min_lfc, sigthresh = max_fdr,
							textcx=.8,  name = output_volcano)
	
	## Plot: DEG
	sig_idx <- resdata$padj < max_fdr & abs(resdata$log2FoldChange_shrunken) > min_lfc
	sig_idx[is.na(sig_idx)] <- FALSE
	resdata.sig <- resdata[sig_idx, ]
	
	n1 <- dim(resdata.sig)[2]
	n2 <- dim(norm_exp)[2]
	zscore.df <- zscore(resdata.sig[, (n1 - n2 + 1):n1])
	rownames(zscore.df) <- resdata.sig[, 1]
	output_heatmap <- paste0("DEG_heatmap/", name, ".heatmap")
	label_rows <- anno[anno$Gene %in% rownames(zscore.df), "Name"]
	Heatmap(zscore.df, nclass = 1, coldata = coldata,
					fname = output_heatmap, show_rownames = FALSE,
					labels_row = label_rows,
					main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr ))
	
	## Plot: DEG using mean value only
	zscore_mean.df <- get_zscore_mean(zscore.df, meta.df$GROUP_LABEL)
	output_heatmap_mean <- paste0("DEG_heatmap_mean/", name, ".heatmap")
	
	Heatmap(zscore_mean.df, nclass = 1, coldata = NULL, # use coldata = NULL to skip adding sample names
					fname = output_heatmap_mean, show_rownames = FALSE,
					labels_row = label_rows,
					main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr ))
	
	## Plot: DEG using mean value and with rowname
	output_heatmap_mean_rowname <- paste0("DEG_heatmap_mean_rowname/", name, ".heatmap")
		### match ENSEMBL to SYMBOL
	label_rows <- anno[anno$Gene %in% rownames(zscore_mean.df), "Name"]
	Heatmap(zscore_mean.df, show_rownames = TRUE,
					labels_row = label_rows, labels_n = 50,
					fname = output_heatmap_mean_rowname,
					main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr))
}

## Plot: heatmap for custom gene sets
geneset_names <- c("MHC_NK", "NEW_SENESCENE", "SENESCENCE_SASP")
for (geneset_name in geneset_names) {
	# geneset_name <- "MHC_NK"
	message("processing geneset: ", geneset_name)
	genesets <- read_gmt(gmt_filename = paste0("../gmt/gmt/", geneset_name, ".gmt"))
	for (collection_name in names(genesets)) {
		# collection_name <- "GO_MHC_CLASS_II_PROTEIN_COMPLEX"
		tem_res <- resdata[resdata$Name %in% genesets[[collection_name]], ]
		n1 <- dim(tem_res)[2]
		n2 <- dim(norm_exp)[2]
		zscore.df <- zscore(tem_res[, (n1 - n2 + 1):n1])
		rownames(zscore.df) <- tem_res[, 1]
		message("collection: ", collection_name)
		
		## plot: geneset
		output_heatmap_geneset <- paste0("Geneset_heatmap/", geneset_name, "/", collection_name, ".heatmap")
		labels_row <- anno[anno$Gene %in% rownames(zscore.df), "Name"]
		Heatmap(zscore.df, nclass = 1, coldata = coldata,
						fname = output_heatmap_geneset, show_rownames = FALSE,
						labels_row = labels_row,
						main = collection_name)
	
		## plot: geneset using mean value
		zscore_mean.df <- get_zscore_mean(zscore.df, meta.df$GROUP_LABEL)
		output_heatmap_geneset_mean <- paste0("Geneset_heatmap_mean/", geneset_name, "/", collection_name, ".heatmap")
		Heatmap(zscore_mean.df, nclass = 1, coldata = NULL, # use coldata = NULL to skip adding sample names
						fname = output_heatmap_geneset_mean, show_rownames = FALSE,
						labels_row = labels_row,
						main = collection_name)
		
		## plot: geneset using mean value with rowname
		output_heatmap_geneset_mean_rowname <- paste0("Geneset_heatmap_mean_rowname/", geneset_name, "/", collection_name, ".heatmap")
		Heatmap(zscore_mean.df, nclass = 1, coldata = NULL, # use coldata = NULL to skip adding sample names
						fname = output_heatmap_geneset_mean_rowname, show_rownames = TRUE,
						labels_row = labels_row,
						main = collection_name, labels_n = 55)
		
	}
}
