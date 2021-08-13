library(edgeR)
library(PCAtools) # loads tools for making heatmaps

sample_info <- read.csv("sample_info.csv", row.names = 1)
raw_counts <- read.csv("filtered_counts.csv", row.names = 1, stringsAsFactors = FALSE)

# Remove any rows from 'sample_info' not corresponding to samples present in 'raw_counts'
sample_info <- sample_info[rownames(sample_info) %in% colnames(raw_counts),]
# Ensure the columns of 'raw_counts' match the rows of 'sample_info', and are in the same order
raw_counts <- raw_counts[,rownames(sample_info)]
nrow(raw_counts)

#
dgList_NASH <- DGEList(
  counts = raw_counts, # counts data
  samples = sample_info, # sample data 
  group=sample_info$NASH # specify experimental groups
)

dgList_NASH <- dgList_NASH[, which(dgList_NASH$samples$control != "Y")]
dgList_NASH <- dgList_NASH[, which(dgList_NASH$samples$remove != "TRUE")]
dgList_NASH <- dgList_NASH[, which(dgList_NASH$samples$lib.size > 1000000)]
# 693 sample left

dgList_NASH_unfiltered<-dgList_NASH # we will make a copy of the complete dataset
nrow(dgList_NASH_unfiltered)

keep <- filterByExpr(dgList_NASH_unfiltered, group=dgList_NASH_unfiltered$sample$group)
dgList_NASH <- dgList_NASH_unfiltered[keep,]
nrow(dgList_NASH)
ncol(dgList_NASH)

help(biplot)

counts_for_pca_NASH <-cpm(dgList_NASH$counts,log=TRUE,prior.count=1)
pca_output_NASH <- pca(counts_for_pca_NASH, metadata = dgList_NASH$samples)
png("NASHPCA.png")
biplot(pca_output_NASH, colby = 'group', pointSize = 1, drawConnectors = 'False', legendPosition = 'bottom', title = 'PCA plot of NASH samples')
dev.off()

dgList_NASH <- calcNormFactors(dgList_NASH)
design <- model.matrix(~ 0+ group, data = dgList_NASH$samples)
dgGlm <- estimateDisp(dgList_NASH, design, robust = TRUE)

fit <- glmQLFit(dgGlm, design, robust = TRUE)
contrast_name<-'groupNASH-groupNAFLD'
contrast_matrix <- makeContrasts(contrasts=contrast_name,levels=design)
de <- glmQLFTest(fit, contrast=contrast_matrix)
top_genes <- topTags(de, n=20)
top_genes

NASH_result <- topTags(de, n=nrow(dgList_NASH), sort.by='none')$table
write.csv(NASH_result, file='differential_NASH_result.csv')
fdr_threshold <- 0.05
fc_threshold <- 2
diffexp_genes <- rownames(NASH_result)[abs(NASH_result$logFC) >= log2(fc_threshold) & NASH_result$FDR <= fdr_threshold ]
print(paste(length(diffexp_genes), 'genes are differentially expressed at a fold change of at least',fc_threshold, 'and a maximum FDR of',
            fdr_threshold))

NASH_result[which(NASH_result$logFC >= 1 & NASH_result$FDR <= 0.05),'Significance'] <- 'Up'
NASH_result[which(NASH_result$FDR <= 0.05 & NASH_result$logFC <= -1),'Significance'] <- 'Down'
NASH_result[!(NASH_result$Significance %in% c('Up', 'Down')),'Significance'] <- 'None'
NASH_result[(NASH_result$Significance %in% c('Up', 'Down')),]
write.csv(NASH_result, file='differential_NASH_result.csv')

volcano_plot<- function(NASH_result_table, fc_threshold=2,
                        fdr_threshold=0.05,label_threshold=9){
  genes<-NASH_result
  ggplot(data=genes, aes(logFC, -log10(FDR))) +
    geom_point(aes(col=Significance),size=0.2)+
    geom_text_repel(aes(label=ifelse(abs(logFC)>=label_threshold,as.character(row.names(genes)),'')))+ 
    scale_colour_manual(name = 'Significance', values = setNames(c('blue','grey','red'),c('Down', 'None','Up'))) +
    geom_vline(xintercept=c(log2(fc_threshold),log2(fc_threshold)*-1), linetype= "dashed")+
    geom_hline(yintercept= 1.3, linetype= "dashed")+
    labs(title="Volcano plot from group NASH vs NAFLD",x = 'log2 Fold Change')
}

png("v.png")
volcano_plot(NASH_result, fc_threshold,fdr_threshold,8)
dev.off()