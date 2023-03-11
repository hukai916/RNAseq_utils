fix_hyphen <- function(x){
	return(gsub("-", ".", x))
}

readExcel <- function(fname){
	df <- readxl::read_xlsx( fname, na='NA', sheet=1)  # na term important
	df <- data.frame(df)
	return(df)
}

getAnnotation <- function(urlpath) {
	tmp <- tempfile()
	download.file(urlpath, destfile = tmp, method = 'auto')
	return(read.table(tmp, sep="\t", header = TRUE, quote=""))
}

get_density <- function(x, y, ...) {
	dens <- MASS::kde2d(x, y, ...)
	ix <- findInterval(x, dens$x)
	iy <- findInterval(y, dens$y)
	ii <- cbind(ix, iy)
	return(dens$z[ii])
}

calculateTPM <- function(counts, len) {
	# michael's version https://support.bioconductor.org/p/91218/
	x <- counts/len
	return(t(t(x)*1e6/colSums(x)))
}

calculateFPKM <- function(counts, len) {
	x <- counts %>% t(t(counts)*1e6/colSums(counts))
	return (x/len*1e3)
}

zscore <- function(matrix) {
	return( t(scale(t(matrix))))
}

get_zscore_mean <- function(zscore.df, group_label) {
	# assuming the colname(zscore.df) match the order in group_label
	# group_label <- meta.df$GROUP_LABEL
	df <- t(zscore.df) %>% as.data.frame()
	df$Name <- group_label
	df <- df %>% group_by(Name) %>% summarise_each(mean) ## Caution: summarizse_each() reorderds the rows!
	df <- df[match(unique(group_label), df$Name), ] %>% select(-Name) ## revert the row orders back
	zscore_mean.df <- t(as.matrix(df))
	colnames(zscore_mean.df) <- unique(group_label)
	zscore_mean.df
}

parse_name <- function(name) {
	name <- gsub(" ", "", name)
	name <- gsub(";$", "", name)
	names <- strsplit(name, ";") [[1]]
	name <- gsub(";", ".", name)
	name <- name %>% str_replace("-", ".") %>% str_replace("/", ".") # since DESeq2 auto corrected the names
	names <- names %>% str_replace("-", ".") %>% str_replace("/", ".")
	return (list(name=name, names=names))
}

sparse_label <- function(labels = labels, labels_order = NULL, max_n = 50) {
	if (is.null(labels_order)) {
		labels_order <- seq(1:length(labels))
	}
	if (length(labels) != length(labels_order)) {
		stop("labels length must match labels_order!")
	}
	idx <- round(seq(1, length(labels_order), length.out = max_n))
	idx2 <- labels_order[idx]
	labels[setdiff(seq_along(labels), idx2)] <- ""
	labels
}
	


rename_num_vector_by_order <- function(l) {
	# l have to be a vector of numbers
	# output vector of roman numbers ordered by appearance in the input vector
	# e.g. c(2,3,3,2,1) -> c(I, II, II, I, III)
	# test rename_num_vector_by_order(c(2,3,3,2,1))
	u <- unique(l)
	n=0
	for (i in u){
		n = n+1; 
		l <- replace(l, l==i, as.character(as.roman(n)))
	}
	return(l)
}

process_deseq_res <- function(res="lfcshrink.res", res2="results.res", name='name', anno='anno.df', norm_exp="tpm.df") {
	## Summary
	print(name)
	print("\n>>> Summary using FDR cut-off only (LFC not used)")
	summary(res, alpha=max_fdr)
	
	print("\n>>> Summary using both FDR and LFC_shrunken cut-off")
	sig_idx <- res$padj<max_fdr & abs(res$log2FoldChange) > min_lfc
	sig_idx[is.na(sig_idx)] <- FALSE
	res_sig <- res[sig_idx,]
	print(table(sig_idx))
	
	up_idx <- res$padj<max_fdr & res$log2FoldChange > min_lfc
	up_idx[is.na(up_idx)] <- FALSE
	res_sig <- res[up_idx,]
	print(table(up_idx))
	
	down_idx <- res$padj<max_fdr & res$log2FoldChange < -min_lfc
	down_idx[is.na(down_idx)] <- FALSE
	res_sig <- res[down_idx,]
	print(table(down_idx))
	
	# Prep
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
	head(resdata)
	sig_idx <- resdata$padj<max_fdr & abs(resdata$log2FoldChange_shrunken) > min_lfc # important to put this line right before output sig.xlsx
	sig_idx[is.na(sig_idx)] <- FALSE
	resdata.sig <- resdata[sig_idx,]
	head(resdata.sig)
	
	## Write results
	if (!TEST){
		WriteXLS(x = resdata,
						 ExcelFileName = paste(name, 'deseq2.xlsx', sep = '.'),
						 row.names = F, SheetNames = 'sheet1', na = 'NA')  # for user
		
		WriteXLS(x = resdata.sig,
						 ExcelFileName = paste(name, 'deseq2.sig.FDR', max_fdr, 
						 											'LFC', min_lfc, 'xlsx', sep = '.'),
						 row.names = F, SheetNames = 'sheet1', na = 'NA')  # for user
	}
	
	## For GSEA
	rnk <- subset(resdata, select = c("Name","log2FoldChange_shrunken"))
	colnames(rnk) <- c("# Name","log2FoldChange_shrunken")
	rnk <- rnk[order(rnk$log2FoldChange_shrunken), ]
	rnk[, 1] <- toupper(rnk[, 1])
	write.table(rnk, 
							paste('rnk/', name, '.rnk', sep = ''), 
							row.names = F, quote = F, sep='\t')
	
	##  Plots
	ggplot(data.frame(res), aes(x=pvalue))+
		geom_histogram(color="darkblue", fill="lightblue")
	ggsave(paste0(name, '.pvalue.pdf'))
	
	ggplot(data.frame(res), aes(x=padj))+
		geom_histogram(color="darkblue", fill="lightblue")
	ggsave(paste0(name, '.fdr.pdf'))
	
	maplot(res, anno, paste0(name, '.maplot.shrunken_lfc.pdf'))
	maplot(res2, anno, paste0(name, '.maplot.raw_lfc.pdf'))
	
	volcanoplot(res, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
							textcx=.8,  name= paste(name, "LFC_shrunken", sep="."))
	volcanoplot(res2, anno, lfcthresh=min_lfc, sigthresh=max_fdr,
							textcx=.8,name= paste(name, "LFC_raw", sep="."))
	
	n1 <- dim(resdata.sig)[2]
	n2 <- dim(norm_exp)[2]
	zscore.df <- zscore(resdata.sig[, (n1-n2+1):n1])
	rownames(zscore.df) <- resdata.sig[,1]
	colnames(zscore.df) <- gsub (":TPM", "", colnames(zscore.df))
	Heatmap(zscore.df, nclass = 2,
					fname = paste(name, "heatmap", sep="."),
					main = paste(name, "LFC >", min_lfc, "FDR <", max_fdr ))
}

read_gmt <- function(gmt_filename = "../gmt/gmt/MHC_NK.gmt") {
	lines <- readLines(gmt_filename)
	list_of_vectors <- list()
	for (line in lines) {
		components <- strsplit(line, "\\s+")[[1]]
		vector_name <- components[1]
		vector_elements <- components[-c(1, 2)]
		list_of_vectors[[vector_name]] <- as.vector(vector_elements)
	}
	list_of_vectors
}

