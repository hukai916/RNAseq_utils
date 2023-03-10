volcanoplot <- function(res, anno, name='name', 
												lfcthresh=2, sigthresh=0.05, 
												labelsig=FALSE, xlim=100, ylim=1000, textcx=1) {
	# If Annotation file is wrong, No annotation match
	res0 <- res
	res<-merge(data.frame(res0), anno, by.x = 0, by.y=1, all.x=T, sort=F)
	if (sum(anno[, 1] %in% row.names(res0)) < 1) {
		warning(
			c("\nThe annotation file and the count filt does have no match in gene-id:\n", 
				"count table gene-id: ", head(row.names(res0), 1), 
				"\nanno table gene-id: ", anno[1:1, 1], 
				"gene-id rather than gene-name used for Volcano-plot"
			))
		res$Name <- res[, 1]
	}
	
	# remove NA
	res$padj[is.na(res$padj)] <- 1
	# set lim on x, y
	res$padj[res$padj < 10^(-ylim) & !is.na(res$padj)] <- 10^(-ylim) # y-axis top value 50
	res$log2FoldChange[res$log2FoldChange > xlim] <- xlim
	res$log2FoldChange[res$log2FoldChange < -xlim] <- -xlim
	# show num_pos num_neg
	pos <- subset(res, padj<sigthresh & log2FoldChange>lfcthresh)
	neg <- subset(res, padj<sigthresh & log2FoldChange< -lfcthresh)
	pos.n <- dim(pos)[1]
	neg.n <- dim(neg)[1]
	
	EnhancedVolcano(res,
									lab = res$Name,
									#selectLab = as.character(res$Name[which(res$padj<labelcut)]), # mark top genes
									#selectLab = c("FOS", "LDHA"), # mark selected genes
									x = 'log2FoldChange',
									y = 'padj',
									title = name,
									subtitle = paste("Up:", pos.n, ", Down:", neg.n, sep = ""),
									xlab = bquote(~Log[2]~ "Fold Change"),
									ylab = bquote(~-Log[10]~italic(FDR)),
									pCutoff = max_fdr,
									FCcutoff = min_lfc,
									cutoffLineType = 'twodash',
									cutoffLineWidth = 0.8,
									legendLabels = c('NS', expression(Log[2]~FC),
																	 "FDR", expression(FDR~and~Log[2]~FC)),
									caption = paste0('Total = ', nrow(res), ' genes'),
									legendPosition = 'right',
									legendLabSize = 10,
									axisLabSize = 10,
									legendIconSize = 3.0)
	ggsave(paste(name, "pdf", sep="."), width=8, height=6)
}

maplot <- function(res, anno, fname){
	maplotdat <- merge(anno, data.frame(res), by.x=1, by.y=0)
	ggmaplot(maplotdat, main = name,genenames = as.vector(maplotdat$Name),
					 fdr = max_fdr, fc = 2^min_lfc, size = 0.4,
					 xlab = "Log2 Mean Expression",
					 ylab = "Log2 Fold Change",
					 palette = c("#B31B21", "#1465AC", "darkgray"),
					 legend = "top", top = 20,
					 font.label = c("bold", 11),label.rectangle = TRUE,
					 font.legend = "bold",
					 font.main = "bold",
					 ggtheme = ggplot2::theme_minimal())
	ggsave(fname)
}

plotQC_PCA <- function(dds, fname = NULL, ntop = 500, height = 8, width = 8) {
	vsd <- varianceStabilizingTransformation(dds)
	batch <- vsd$batch
	sample <- vsd$sample
	pcaData <- plotPCA(vsd, intgroup = 'group', returnData = TRUE, ntop = ntop)
	percentVar <- round(100 * attr(pcaData, 'percentVar'), 1)
	if (!(length(levels(batch)) > 1)) {
		batch <- NULL
	}
	
	res <- ggplot(pcaData, aes(PC1, PC2, color = group, shape = batch)) +
		geom_point(size = 3) +
		xlab(paste0("PC1: ", percentVar[1], "% variance")) +
		ylab(paste0("PC2: ", percentVar[2], "% variance")) +
		geom_label_repel(aes(label = sample),
										 box.padding = 0.35,
										 point.padding = 1,
										 segment.color = 'grey50',
										 segment.alpha = 0.5,
										 show.legend = FALSE) + # if TRUE, legend not right
		theme_classic()
	if (!is.null(fname)) {
		ggsave(fname, res, height = height, width = width)
	} else {
		print(res)
	}
}

Heatmap <- function(df, fname="heatmap", main = "title",
										coldata = NULL, # if NULL, use colname from df to annotate
										show_rownames = FALSE,
										labels_row = rownames(df),
										labels_n = 50,
										cellheight = NA,
										nclass = NA,
										fontsize_row = 10) {
	# df <- zscore_mean.df
	# coldata <- NULL
	# nclass <- NA
	# main <- "test"
	# fname <- "test"
	# labels_row <- rownames(df)
	# show_rownames <- FALSE
	# cellheight <- NA
	# labels_n <- 50
	
	# df <- zscore.df
	# nclass = 1
	# coldata = coldata # use coldata = NULL to skip adding sample names
	# fname = output_heatmap_geneset_mean_rowname
	# show_rownames = TRUE
	# labels_row = labels_row
	# main = collection_name

	if (dim(df)[1] < 1) {
		message("Too few rows to plot heatmap!")
		output_notice <- paste0(fname,".no_gene.txt")
		file.create(output_notice)
		return(1)
	}
	
	if (is.na(nclass)) {
		nclass <- 1
	}

	if (dim(df)[1] > 1) {
		cluster_rows <- TRUE
	} else {
		cluster_rows <- FALSE
	}
	
	# for annotating cols: if is.null(coldata), set column name according to colnames(df)
	if (is.null(coldata)) {
		my_sample_col <- data.frame(group = colnames(df))
		row.names(my_sample_col) <- colnames(df)
	} else {
		my_sample_col <- data.frame(group = coldata[, 'group'])
		row.names(my_sample_col) <- coldata[, 'sample']
	}
	
	if (labels_n == "all") {
		labels_n <- length(labels_row)
	}
	
	# first-round pheatmap to retrieve the clustering information
	p <- pheatmap(df, 
								annotation_col = my_sample_col,
								main = main,
								cluster_cols = F, 
								cluster_rows = cluster_rows,
								border_color = NA,
								cutree_rows = nclass, # split into nclass part horizontally for better visualization
								show_rownames = show_rownames,
								labels_row = labels_row,
								cellheight = cellheight,
								# fontsize_row = fontsize_row,
								filename = NA)
	
	# sparse the labels: default to display a max of 50 labels
	# print(labels_row)
	# print(p$tree_row$order)
	# print(dim(df))
	if (length(labels_row) > 1) { # if length(labels_row) == 1, p$tree_row$order does not exist
		labels_row_sparse <- sparse_label(labels = labels_row, labels_order = p$tree_row$order, max_n = labels_n)
	} else {
		labels_row_sparse <- labels_row
	}

	# save zscores used to create heatmap to Excel
	output_zscore <- paste0(fname,".zscores.xlsx")
	tem_df <- as.data.frame(df)
	tem_df$GENE_ID <- rownames(tem_df)
	tem_df$GENE_SYMBOL <- labels_row
	tem_df <- select(tem_df, c(GENE_ID, GENE_SYMBOL), everything())
	WriteXLS(tem_df, row.names = FALSE, output_zscore)
	
	# output final plots
	p2 <- pheatmap(df, annotation_col = my_sample_col,
								 main = main, border_color = NA, 
								 cluster_cols = FALSE, cluster_rows = cluster_rows,
								 cutree_rows = nclass,
								 show_rownames = show_rownames,
								 labels_row = labels_row_sparse,
								 filename = paste0(fname, ".pdf"))
	
	p3 <- pheatmap(df, annotation_col = my_sample_col,
								 main = main, border_color = NA,
								 cluster_cols = TRUE, cluster_rows = cluster_rows,
								 cutree_rows = nclass,
								 show_rownames = show_rownames,
								 labels_row = labels_row_sparse,
								 filename = paste0(fname, ".cluster_col.pdf"))
}
