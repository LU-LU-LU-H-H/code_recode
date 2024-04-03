library(ggplot2)
library(ggrepel)

# Rscript 1.work.R ../supp/DEG/F___E venn_mRNA/all/F__E___G__E/F__E___G__E.gene.txt 6.F__E.pdf
##print("""Rscript volcano_with_name.R    ../supp/DEG/F___E    gene_list_file    out_file""")

logFC <- 0.58
PValue <- 0.05
args <- commandArgs(TRUE)
f_d <- args[1] # Gene_expr_with_pval.xls
gene_f <- args[2]
out_f <- args[3]

f <- paste0(f_d, "/Gene_expr_with_pval.xls")
up <- paste0(f_d, "/sig_UP_genes.xls")
down <- paste0(f_d, "/sig_Down_genes.xls")

UpNames <- read.table(up, sep='\t', header=T)$Geneid
DownNames <- read.table(down, sep='\t', header=T)$Geneid
out <- read.table(f, sep='\t', header=T)
rownames(out) <- out$Geneid
gene <- read.table(gene_f, header=F)

colors <- rep("grey80", nrow(out))
##colors[(rownames(out) %in% UpNames) & (rownames(out) %in% gene$V1)] <- "red"
##colors[(rownames(out) %in% DownNames) & (rownames(out) %in% gene$V1)] <- "blue"
colors[(rownames(out) %in% UpNames) ] <- "red"
colors[(rownames(out) %in% DownNames) ] <- "blue"
colors[(rownames(out) %in% gene$V1)] <- "black"

tmp <- out[match(gene$V1 ,out$Geneid),]
print(length(tmp))

px <- out$logFC
px <- px[is.finite(px)]
py <- -log10(out$PValue)
py <- py[is.finite(py)]

pc <- max(out[c(UpNames, DownNames),]$PValue)
pdf(out_f) # 'volcano.pdf'
#ggplot(out, aes(logFC, -log10(PValue))) + geom_point(colour=colors) + xlim(min(px),max(px)) + ylim(0,max(py)) + geom_hline(yintercept=-log10(pc), linetype = "dashed") + geom_vline(xintercept=min(out[UpNames,]$logFC), linetype = "dashed") + geom_vline(xintercept=max(out[DownNames,]$logFC), linetype = "dashed") + geom_text_repel(data=out[match(gene$V1 ,out$Geneid),], aes(label=Geneid),  box.padding = unit(0.35, "lines"),point.padding = unit(0.3, "lines")) + labs(x='log2(fold change)') + labs(y='-log10(p-value)') + labs(title='Volcano Plot') + theme(plot.title = element_text(hjust = 0.5))
ggplot(out, aes(logFC, -log10(PValue))) + 
geom_point(colour=colors) + 
xlim(min(px),max(px)) + 
ylim(0,max(py)) + 
geom_hline(yintercept=-log10(PValue), linetype = "dashed") + 
geom_vline(xintercept=-logFC, linetype = "dashed") + 
geom_vline(xintercept=logFC, linetype = "dashed") + 
geom_text_repel(data=out[match(gene$V1 ,out$Geneid),], aes(label=Geneid),  box.padding = unit(0.8, "lines"),point.padding = unit(0.3, "lines")) + 
labs(x='log2(fold change)') + 
labs(y='-log10(p-value)') + 
labs(title='Volcano Plot') + 
theme(plot.title = element_text(hjust = 0.5),
panel.background = element_rect(fill="white",colour = "black"),
)
dev.off()

