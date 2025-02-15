
library(Seurat)
library(ggplot2)
library(ggsignif)
library(ggrastr)

Cluster_colors <- c(
	'#A4CEE2',
	'#1F77B5',
	'#2B6E38', 
	'#BBD58D', 
	'#842758', 
	'#E1BB6D', 
	'#C3B1CC', 
	'#673F95', 
	'#010000', 
	'#E99C9A', 
	'#E0A630', 
	'#D33C28', 
	'#992223', 
	'#A76031')


names(Cluster_colors) <- c(
	'HSC', 
	'LMPP', 
	'GMP', 
	'Granulocyte',
	'Monocytes', 
	'DendriticCell', 
	'CLP', 
	'pro-B', 
	'T', 
	'MEP', 
	'MK_Prog', 
	'EarlyErythroid',
	'LateErythroid', 
	'Basophil')


get_volcano_comparison <- function(df, plot_title=NULL, genes_2_plot = gene_list) {
    if(!'gene' %in% colnames(df)){df$gene <- rownames(df)}
    df$color <- ifelse(df$p_val_adj > 0.05, "Gray", ifelse(df$avg_log2FC > 0, "Red", "Blue"))
    
    p <- ggplot(df, aes(x = avg_log2FC, y = -log(p_val_adj), color=color)) +
        geom_point_rast(alpha = 0.9, size = 0.8) +
        ggrepel::geom_label_repel(data = df[df$gene %in% genes_2_plot & df$p_val_adj<0.05, ], 
                                  aes(x = avg_log2FC, y = -log(p_val_adj),label = gene),
                                  color='black', segment.colour = "black",fontface = 'plain', size = 1,
                                  min.segment.length = unit(0, 'lines'), max.overlaps = 75,
                                  nudge_y = 0.1, family = "Helvetica") + 
        scale_color_manual(values = c("Red" = "#8C3F35", "Blue" = "#7C93DA", "Gray" = "#BEBEBE")) +
        labs(x = "Average log2FC",
             y = "-log(p_val_adj)",
             color = "Gene Color") +
        theme_minimal() + 
        theme(legend.position='none', 
              axis.text = element_text(size = 5, family='Helvetica'),
              axis.title = element_text(size = 7, family='Helvetica'),
              title = element_text(size = 7, family='Helvetica'))
    if(is.null(plot_title)){
        p <- p + labs(title = plot_title)
    }
    return(p)
}



########################################################################################################################
# EXPLORATORY ANALYSIS FOR THE CSNK1A1
##########################################

phenotypes <- c('del(5q)', 'non-del(5q)', 'Healthy')
color_map <- c('#003049', '#d62828', '#f77f00')
names(color_map) <- phenotypes

elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
norm_data_elder <- elder_data@assays$RNA@data
norm_data_elder <- t(as.data.frame(norm_data_elder['CSNK1A1', , drop=FALSE]))
plotter_elder <- FetchData(elder_data, c('Cluster_names', 'Sample'))
plotter_elder$cell_5q <- 'Healthy'
plotter_elder <- merge(plotter_elder, norm_data_elder, by=0)

sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'

plotter_5q <- FetchData(sc_data_MDS5q_subseted, c('cell_5q', 'Cluster_names', 'Sample'))
nomr_data_5q <- sc_data_MDS5q_subseted@assays$RNA@data
nomr_data_5q <- t(as.data.frame(nomr_data_5q['CSNK1A1', , drop=FALSE]))
plotter_5q <- merge(plotter_5q, nomr_data_5q, by=0)


plotter <- rbind(plotter_5q, plotter_elder)
plotter$cell_5q <- ifelse(plotter$cell_5q == 'normal', 'non-del(5q)', 
                   ifelse(plotter$cell_5q == 'del5q', 'del(5q)', 'Healthy'))
plotter$cell_5q <- factor(plotter$cell_5q , levels=c('Healthy', 'non-del(5q)',  'del(5q)'))

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/CSNK1A1_Expression.pdf')
cowplot::plot_grid(
        ggplot(plotter, aes(x= cell_5q, y= CSNK1A1, fill=cell_5q)) + 
        geom_violin(alpha=0.7) +
        scale_fill_manual(values=color_map) + 
        theme_classic() + theme(legend.position = 'none') +
        geom_signif(comparisons = list(
            c("del(5q)", "non-del(5q)"), 
            c("Healthy", "non-del(5q)"), 
            c("del(5q)", "Healthy")), 
            step_increase = 0.1,
            map_signif_level=TRUE)
    ,
        ggplot(plotter, aes(x= cell_5q, y= CSNK1A1, fill=cell_5q)) + 
        geom_violin(alpha=0.7) +
        scale_fill_manual(values=color_map) + 
        theme_classic() + theme(legend.position = 'none') +
        ylim(c(0, 5)) +
        geom_signif(comparisons = list(c("del(5q)", "non-del(5q)")),   
                map_signif_level=TRUE) + facet_wrap(~Sample, nrow=2)
    ,nrow=2)
dev.off()




########################################################################################################################
# SERIOUS ANALYSIS FOR THE CSNK1A1
########################################################################################################################
library(Seurat)
library(Libra)
library(ggplot2)

get_prop_tables_per_sample <- function(Sobj){
    print(
        prop.table(table(Sobj$Sample, Sobj$cell_5q), margin=1)*100
    )
}

post_data <-  readRDS('/home/tereshkova/data/gserranos/MDS/Data/POST_Samples_Annotated_final.rds')
all_5q_depleted_cells_COPYKAT <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/CASPER/all_5q_depleted_cells_CASPER_AND_COPYKAT_POST.rds')
post_data$cell_5q <- ifelse(colnames(post_data) %in% all_5q_depleted_cells_COPYKAT, 'del5q', 'normal')
#  There are less cells as some clusters have disapeared based on the annotation
# elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')
MDS_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
pre_post_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/pre_post_5q_Annotated_final.rds')
pre_post_data_del <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/SelectedCells5q_PrePost_CASPER_COPYKAT.rds')
pre_post_data$cell_5q <- ifelse(colnames(pre_post_data) %in% pre_post_data_del, 'del5q', 'normal')
elder_data <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/elder_Samples_Annotated_final.rds')



# nonDel(5q) vs Del(5q) MDs
get_prop_tables_per_sample(MDS_data)

Idents(MDS_data) <- 'cell_5q'
MDS_normal     <- subset(x = MDS_data, subset = cell_5q == "normal")
expression_MDS_normal<- as.data.frame(MDS_normal@assays$RNA@counts)
MDS_normal$Status <- 'non-del(5q)'

MDS_5q     <- subset(x = MDS_data, subset = cell_5q == "del5q")
expression_MDS_5q <- as.data.frame(MDS_5q@assays$RNA@counts)
MDS_5q$Status <- 'del(5q)'

genes_2_keep <- intersect(rownames(expression_MDS_5q), rownames(expression_MDS_normal))
'CSNK1A1' %in% genes_2_keep


all_data <- merge(expression_MDS_5q[genes_2_keep,], expression_MDS_normal[genes_2_keep,], by=0)
rownames(all_data) <- all_data$Row.names
all_data$Row.names <- NULL

ncol(all_data) == ncol(expression_MDS_5q) + ncol(expression_MDS_normal)

metadata <- setNames(rbind(MDS_5q[[c('Status', 'Sample', 'Cluster_names')]], 
                           MDS_normal[[c('Status', 'Sample', 'Cluster_names')]]), 
                           c('label', 'replicate', 'cell_type'))

metadata$label <- factor(metadata$label, levels=c('non-del(5q)', 'del(5q)'))
table(metadata$replicate, metadata$label)
table(metadata$cell_type, metadata$label)
ncol(all_data) == nrow(metadata)

DE = run_de(as.matrix(all_data), meta = metadata, min_reps=1)
de_results <- as.data.frame(DE)
de_results[de_results$gene =='CSNK1A1',]




########################################################################################################################
# Comparison of the DE results
########################################################################################################################
library(Seurat)
library(Libra)
library(RColorBrewer)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(ggrepel)
library(colorRamp2)
# library(enrichR)

# 1 the number of shared DE genes accross all comparisons
sc_data_MDS5q_subseted <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/5qSamples_Annotated_final.rds')
expression <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@counts)
metadata <-	 setNames(as.data.frame(sc_data_MDS5q_subseted$Cluster_names), 'cell_type')
metadata$label <- sc_data_MDS5q_subseted$cell_5q
metadata$replicate <- sc_data_MDS5q_subseted$Sample

metadata$label <- factor(metadata$label, levels = c('del5q', 'normal'))
metadata <- metadata[rownames(metadata) %in% colnames(expression), ]

DE = run_de(as.matrix(expression), meta = metadata)
DE <- as.data.frame(DE)
DE <- setNames(DE, c( "cell_type", "gene","avg_log2FC", "p_val", "p_val_adj", "de_family", "de_method", "de_type"))
DE_plotter <- DE[DE$gene %in% c('PRSS21', 'ZCCHC10', 'FIGN', 'SIL1',  'CCL5', 'MAP3K7CL',  'IGKC'), 
                c('cell_type', 'gene', 'avg_log2FC', 'p_val_adj')]


gene_list <- c('PRSS21', 'ZCCHC10', 'FIGN', 'SIL1',  'CCL5', 'MAP3K7CL',  'IGKC')

results_wilcox <- list()
Idents(sc_data_MDS5q_subseted) <- 'Cluster_names'
DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
for(cell_type in sort(unique(metadata$cell_type))){
    print(cell_type)
    tmp <- subset(sc_data_MDS5q_subseted, idents = cell_type)
    Idents(tmp) <- 'cell_5q'
    results_wilcox[[cell_type]] <- FindMarkers(tmp, 
                                    ident.1 = 'del5q', 
                                    ident.2 = 'normal',
                                    min.cells.group = 1, 
                                    min.cells.feature = 1,
                                    min.pct = 0,
                                    logfc.threshold = 0,
                                    only.pos = FALSE)
}

# saveRDS(results_wilcox, '/home/tereshkova/data/gserranos/MDS/Data/Revision/DE_results_Wilcox_NoFilters.rds')
# saveRDS(DE, '/home/tereshkova/data/gserranos/MDS/Data/Revision/DE_results_edgeR.rds')

results_wilcox <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/Revision/DE_results_Wilcox_NoFilters.rds')
DE <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/Revision/DE_results_edgeR.rds')

plotter_list <- list()
for(cell_type in sort(unique(metadata$cell_type))){
    plotter_list[[paste0('Wilcox: ', cell_type)]] <-  get_volcano_comparison(
                        results_wilcox[[cell_type]], 
                        plot_title=paste0('Wilcox: ', cell_type))
    plotter_list[[paste0('edgeR: ', cell_type)]] <-  get_volcano_comparison(
                        DE[DE$cell_type == cell_type,], 
                        plot_title=paste0('edgeR: ', cell_type))
}

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/DE_Comparisons.pdf', height = 15)
cowplot::plot_grid(plotlist=plotter_list, ncol=4)
dev.off()

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/DE_Number_perCT.pdf')
plotter <- data.frame(celltype=character(), Ngenes_Wilcox=integer(), Ngenes_edgeR=integer())
for(cell_type in sort(unique(metadata$cell_type))){
    tmp <- results_wilcox[[cell_type]]
    plotter <- rbind(plotter, 
                    data.frame( celltype = cell_type, 
                                Ngenes_Wilcox = nrow(tmp[ tmp$p_val_adj<0.05,]),
                                Ngenes_edgeR = nrow(DE[DE$cell_type == cell_type & 
                                                    DE$p_val_adj<0.05,])))
}
ggplot(plotter, aes(x=Ngenes_edgeR, y=Ngenes_Wilcox, label=celltype)) + 
geom_point() + geom_label_repel() + theme_classic()
dev.off()

plotter <- data.frame(gene=character(), cell_type=character(), 
                      de_method=character(), p_val_adj=numeric(), 
                      avg_log2FC=numeric())
for(cell_type in sort(unique(metadata$cell_type))){
    print(cell_type)
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_wilcox$gene <- rownames(tmp_wilcox)
    tmp_wilcox$de_method <- 'Wilcox'
    tmp_wilcox$cell_type <- cell_type
    plotter <- rbind(plotter, tmp_wilcox[, c('gene', 'cell_type', 'de_method', 
                                           'p_val_adj', 'avg_log2FC')])
    tmp_edgeR <-  DE[DE$cell_type == cell_type,]
    plotter <- rbind(plotter, tmp_edgeR[, c('gene', 'cell_type', 'de_method', 
                                           'p_val_adj', 'avg_log2FC')])
}

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Pval_logFC_Histograms.pdf', height = 15)
ggplot(plotter, aes(x = avg_log2FC, fill = de_method)) + 
geom_histogram(alpha = 0.8, position="dodge") + 
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
theme_classic() + facet_wrap(~cell_type, , scales = "free", nrow=8)

ggplot(plotter, aes(x = p_val_adj, fill = de_method)) + 
geom_histogram(alpha = 0.8, position="dodge") + 
scale_fill_manual(values=c("#999999", "#E69F00", "#56B4E9")) +
theme_classic() + facet_wrap(~cell_type, , scales = "free", nrow=8)
dev.off()

compare_log2FC <- function(cell_type){
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_wilcox$gene <- rownames(tmp_wilcox)
    genes_in_comon <- intersect(
                            rownames(tmp_wilcox), 
                            DE[DE$cell_type == cell_type,'gene'])
    plotter <- data.frame(gene = genes_in_comon)
    plotter <- merge(plotter, setNames(tmp_wilcox[,c('gene', 'avg_log2FC')], c(c('gene', 'Wilcox_avg_log2FC'))), 
                    by='gene', all.x=TRUE)
    plotter <- merge(plotter, setNames(DE[DE$cell_type == cell_type,c('gene', 'avg_log2FC')], c(c('gene', 'edgeR_avg_log2FC'))), 
                    by='gene', all.x=TRUE)
   p <-  ggplot(plotter, aes(x=Wilcox_avg_log2FC, y=edgeR_avg_log2FC)) + 
    geom_point(shape=21, fill='gray', color='black', alpha = 0.7) + 
    geom_abline() + theme_bw() + 
    labs(title= cell_type, 
         subtitle=paste0('Spearman: ', 
         signif(cor(plotter$Wilcox_avg_log2FC, 
            plotter$edgeR_avg_log2FC, 
            method='spearman'), digits=3)), 
            x='wilcox_log2FC',
            y='edgeR_log2FC')
    return(p)
}

plot_list <- list()
for(cell_type in sort(unique(metadata$cell_type))){
   plot_list[[cell_type]] <- compare_log2FC(cell_type)
}
pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/LogFC_Correlation.pdf', height = 10)
cowplot::plot_grid(plotlist=plot_list, ncol=4)
dev.off()


pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/NDEG_per_NumCell.pdf')
plotter <- data.frame(cell_type=character(), Ngenes=integer())
for(cell_type in sort(unique(metadata$cell_type))){
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_wilcox <- tmp_wilcox[tmp_wilcox$p_val_adj < 0.05,]
    plotter <- rbind(plotter, data.frame(cell_type=cell_type, Ngenes=nrow(tmp_wilcox)))
}
plotter <- merge(plotter, setNames(as.data.frame(as.matrix(table(metadata$cell_type))), c('Ncells')),
                by.x='cell_type', by.y=0)
ggplot(plotter, aes(x=Ncells, y=Ngenes, color = cell_type, label = cell_type)) + 
geom_smooth(method="lm", aes(group=1))+ 
geom_point(size =3) + theme_classic() + ggrepel::geom_label_repel() + 
scale_color_manual(values=Cluster_colors) + theme(legend.position='none') +
labs(title='Correlation between NDEG and Ncells', subtitle=paste0('spearman: ', cor(plotter$Ngenes, plotter$Ncells, method='spearman')))

dev.off()


pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/FDR_logFC_DEG.pdf', height=5)
plotter <- data.frame(gene=character(), cell_type=character(), 
                      de_method=character(), p_val_adj=numeric(), 
                      avg_log2FC=numeric())
for(cell_type in sort(unique(metadata$cell_type))){
    print(cell_type)
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_wilcox$gene <- rownames(tmp_wilcox)
    tmp_wilcox$de_method <- 'Wilcox'
    tmp_wilcox$cell_type <- cell_type
    plotter <- rbind(plotter, tmp_wilcox[, c('gene', 'cell_type', 'de_method', 
                                           'p_val_adj', 'avg_log2FC')])
    tmp_edgeR <-  DE[DE$cell_type == cell_type,]
    plotter <- rbind(plotter, tmp_edgeR[, c('gene', 'cell_type', 'de_method', 
                                           'p_val_adj', 'avg_log2FC')])
}

plotter <- plotter[plotter$p_val_adj < 0.05,]
plotter$Direction <- ifelse(plotter$avg_log2FC >0, 'UpDEG', 'DownDEG')
cowplot::plot_grid(
ggplot(plotter, aes(x=abs(avg_log2FC), y=de_method, fill=de_method)) + 
geom_boxplot(width=0.3) + theme_classic() + scale_fill_manual(values=c( "#E69F00", "#999999")),
ggplot(plotter, aes(x=-log10(p_val_adj), y=de_method, fill=de_method)) + 
geom_boxplot(width=0.3) + theme_classic() + scale_fill_manual(values=c( "#E69F00", "#999999")),
nrow=2)
dev.off()


pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/PCT_expr_DEG.pdf')
plotter <- data.frame(gene=character(), cell_type=character(), 
                      pct.1=numeric(), pct.2=numeric(),
                      p_val_adj=numeric(), 
                      avg_log2FC=numeric())
for(cell_type in sort(unique(metadata$cell_type))){
    print(cell_type)
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_wilcox$gene <- rownames(tmp_wilcox)
    tmp_wilcox$cell_type <- cell_type
    plotter <- rbind(plotter, tmp_wilcox[, c('gene', 'cell_type', 'pct.1', 'pct.2', 
                                           'p_val_adj', 'avg_log2FC')])
}
plotter$Ratio <- plotter$pct.1/plotter$pct.2
ggplot(plotter, aes(x = gene, y = Ratio)) + 
coord_flip() + 
geom_bar(stat = "identity", position = "identity", width = 0.525) +
theme_classic()+facet_wrap(~cell_type, scales='free')
dev.off()



# tmp1 <- AggregateExpression(sc_data_MDS5q_subseted, 
#                             return.seurat = F, 
#                             slot = "counts", 
#                             assays = "RNA", 
#                             group.by = c("Cluster_names"))

# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf')
# for(cell_type in sort(unique(metadata$cell_type))){
#     tmp <- results_wilcox[[cell_type]]
#     tmp1$RNA[, cell_type, drop=FALSE]
#     plotter <- merge(tmp, tmp1$RNA[, cell_type, drop=FALSE], by=0)
#     ggplot(plotter, aes(y=.data[[cell_type]], x=p_val_adj)) + geom_point(alpha = 0.8)
# }
# dev.off()



# calculating the coeficient of variation
# Assuming data_replicates is a matrix where each row represents a gene and each column represents a replicate
aggregated_CellType     <- AggregateExpression(sc_data_MDS5q_subseted, 
                            return.seurat = F, 
                            slot = "counts", 
                            assays = "RNA", 
                            group.by = c("Cluster_names"))

aggregated_CellType_Sample <- AggregateExpression(sc_data_MDS5q_subseted, 
                            return.seurat = F, 
                            slot = "counts", 
                            assays = "RNA", 
                            group.by = c("Sample", "Cluster_names"))

aggregated_CellType_Sample_5q <- AggregateExpression(sc_data_MDS5q_subseted, 
                            return.seurat = F, 
                            slot = "counts", 
                            assays = "RNA", 
                            group.by = c("Sample", "Cluster_names", 'cell_5q'))




sc_norm_data <- as.data.frame(sc_data_MDS5q_subseted@assays$RNA@data)
metadata <- FetchData(sc_data_MDS5q_subseted, c('Cluster_names', 'cell_5q'))

results_cv_pseudo_CT     <- data.frame(gene=rownames(aggregated_Samples))
results_cv_sc_CT         <- data.frame(gene=results_cv_pseudo_CT$gene)
results_cv_pseudo_CT_Dis <- data.frame(gene=rownames(aggregated_Samples))
results_cv_sc_CT_Dis     <- data.frame(gene=results_cv_pseudo_CT$gene)

for (celltype in unique(sc_data_MDS5q_subseted$Cluster_names)){
    print(celltype)
    # pseudobulk
    # filter genes with 0 variance
    aggregated_Samples <- aggregated_CellType_Sample$RNA
    sds = matrixStats::rowSds(aggregated_Samples[, grep(celltype, colnames(aggregated_Samples), value=T)])
    aggregated_Samples <- aggregated_Samples[sds>0,]
    aggregated_Samples <- edgeR::cpm(aggregated_Samples)
     # calculate variances for each gene
    gene_vars = matrixStats::rowSds(aggregated_Samples)
    # also calculate mean expression for each gene
    gene_means = rowMeans(aggregated_Samples)
    cov = gene_vars / gene_means
    results_cv_pseudo_CT <- merge(results_cv_pseudo_CT, setNames(as.data.frame(cov), celltype), by.x='gene', by.y=0)
    # single cell
    cells_2_keep <- rownames(metadata[metadata$Cluster_names == celltype, ])
    tmp <- sc_norm_data[results_cv_pseudo_CT$gene, c(cells_2_keep)]
    gene_vars = matrixStats::rowSds(as.matrix(tmp))
    # also calculate mean expression for each gene
    gene_means = rowMeans(tmp)
    cov = gene_vars / gene_means
    results_cv_sc_CT <- merge(results_cv_sc_CT, setNames(as.data.frame(cov), celltype), by.x='gene', by.y=0)
     for (disease in unique(sc_data_MDS5q_subseted$cell_5q)){
        print(paste0(celltype, '--', disease))
        # pseudobulk
        # filter genes with 0 variance
        samples_2_keep <- grep(disease, grep(celltype, colnames(aggregated_Samples), value=T), value=T)
        sds = matrixStats::rowSds(aggregated_Samples[, grep(celltype, colnames(aggregated_Samples), value=T)])
        aggregated_Samples <- aggregated_Samples[sds>0,]
        aggregated_Samples <- edgeR::cpm(aggregated_Samples)
        # calculate variances for each gene
        gene_vars = matrixStats::rowSds(aggregated_Samples)
        # also calculate mean expression for each gene
        gene_means = rowMeans(aggregated_Samples)
        cov = gene_vars / gene_means
        results_cv_pseudo_CT_Dis <- merge(results_cv_pseudo_CT_Dis, setNames(as.data.frame(cov), paste0(celltype, '_', disease)), by.x='gene', by.y=0)
        # single cell
        cells_2_keep <- rownames(metadata[metadata$Cluster_names == celltype & metadata$cell_5q == disease, ])
        tmp <- sc_norm_data[results_cv_pseudo_CT_Dis$gene, c(cells_2_keep)]
        gene_vars = matrixStats::rowSds(as.matrix(tmp))
        # also calculate mean expression for each gene
        gene_means = rowMeans(tmp)
        cov = gene_vars / gene_means
        results_cv_sc_CT_Dis <- merge(results_cv_sc_CT_Dis, setNames(as.data.frame(cov), paste0(celltype, '_', disease)), by.x='gene', by.y=0)

    }
}

results_cv_sc_CT <- reshape2::melt(results_cv_sc_CT)
results_cv_sc_CT$method <- 'singleCell'
results_cv_pseudo_CT <- reshape2::melt(results_cv_pseudo_CT)
results_cv_pseudo_CT$method <- 'pseudobulk'
plotter_CT <- rbind(results_cv_pseudo_CT, results_cv_sc_CT)

results_cv_sc_CT_Dis <- reshape2::melt(results_cv_sc_CT_Dis)
results_cv_sc_CT_Dis$method <- 'singleCell'
results_cv_pseudo_CT_Dis <- reshape2::melt(results_cv_pseudo_CT_Dis)
results_cv_pseudo_CT_Dis$method <- 'pseudobulk'
plotter_CT_Dis <- rbind(results_cv_pseudo_CT_Dis, results_cv_sc_CT_Dis)

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/CoefficientOfVariation.pdf', width = 15, height=5)
ggplot(plotter_CT, aes(x= variable, y=value, fill = method)) + 
geom_boxplot(outlier.alpha = 0.1, outlier.size=1) + scale_fill_manual(values=c( "#E69F00", "#999999")) +
theme_classic() + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + 
labs(title= 'Coeficient of variation per celltype', y='Coeficient of variation', x='Cell types')
ggplot(plotter_CT_Dis, aes(x= variable, y=value, fill = method)) + 
geom_boxplot(outlier.alpha = 0.1, outlier.size=1) + scale_fill_manual(values=c( "#E69F00", "#999999")) +
theme_classic() + theme(axis.text.x = element_text(angle=90, vjust = 0.5, hjust=1)) + 
labs(title= 'Coeficient of variation per celltype and Phenotype', y='Coeficient of variation', x='Cell types per Phenotype')
dev.off()





results_cv_sc <- reshape2::melt(results_cv_sc)
results_cv_sc$method <- 'singleCell'
results_cv_pseudo <- reshape2::melt(results_cv_pseudo)
results_cv_pseudo$method <- 'pseudobulk'
plotter <- rbind(results_cv_pseudo, results_cv_sc)

pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf')
ggplot(plotter, aes(x= variable, y=value, fill = method)) + 
geom_boxplot() + scale_fill_manual(values=c( "#E69F00", "#999999")) +
theme_classic() + theme(axis.text.x = element_text(angle=90)) + 
labs(title= 'Coeficient of variation per celltype', y='Coeficient of variation', x='Cell types')
dev.off()




pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/PCA_Pseudobulk.pdf')
pca_result <- prcomp(t(edgeR::cpm(as.data.frame(aggregated_Samples))))
variance_explained <- summary(pca_result)$importance[2, ] * 100
# Plot PCA
pca_plotter <- as.data.frame(pca_result$x)
pca_plotter$Sample <- stringr::str_extract(rownames(pca_plotter), '^[^_]+(?=_)')
pca_plotter$Celltype  <- stringr::str_extract(rownames(pca_plotter), '(?<=_)[\\w-]+(?=_)')
pca_plotter$Disease <- stringr::str_extract(rownames(pca_plotter), '[^_]+$')
tmp <- pca_plotter
ggplot(pca_plotter, aes(x = PC1, y = PC2, fill=Celltype)) +
			geom_point(alpha=0.8, shape=21, color='black', size = 5) + theme_classic() + theme(legend.position = "none") +
            scale_fill_manual(values=Cluster_colors)+
			labs(title = "PCA of TMM Normalized Samples", 
				x = paste0("PC1: ", round(variance_explained[1], 2), "%"), 
				y = paste0("PC2: ", round(variance_explained[2], 2), "%")
			)

ggplot(pca_plotter, aes(x = PC1, y = PC2, fill=Sample)) +
			geom_point(alpha=0.8, shape=21, color='black', size = 5) + theme_classic() + theme(legend.position = "none") +
			labs(title = "PCA of TMM Normalized Samples", 
				x = paste0("PC1: ", round(variance_explained[1], 2), "%"), 
				y = paste0("PC2: ", round(variance_explained[2], 2), "%")
			) + facet_wrap(~Sample)

ggplot(pca_plotter, aes(x = PC1, y = PC2)) +
			geom_point(data=transform(pca_plotter, Celltype=NULL), alpha=0.5, size = 2, color='grey') + 
			geom_point(alpha=0.8, shape=21, color='black', size = 5, aes(fill=Celltype)) + 
            theme_classic() + theme(legend.position = "none") +
            scale_fill_manual(values=Cluster_colors)+
			labs(title = "PCA of TMM Normalized Samples", 
				x = paste0("PC1: ", round(variance_explained[1], 2), "%"), 
				y = paste0("PC2: ", round(variance_explained[2], 2), "%")
			)+ facet_wrap(~Celltype)

ggplot(pca_plotter, aes(x = PC1, y = PC2, fill=Celltype)) +
			geom_point(alpha=0.8, shape=21, color='black', size = 5) + theme_classic() + theme(legend.position = "none") +
            scale_fill_manual(values=Cluster_colors)+
			labs(title = "PCA of TMM Normalized Samples", 
				x = paste0("PC1: ", round(variance_explained[1], 2), "%"), 
				y = paste0("PC2: ", round(variance_explained[2], 2), "%")
			)+ facet_wrap(Celltype~Disease)

ggplot(pca_plotter, aes(x = PC1, y = PC2, shape=Disease)) +
			geom_point(data=pca_plotter[, -113], alpha=0.5, size = 2, color='grey') + 
			geom_point(alpha=0.8, color='black', size = 5, aes(fill=Celltype)) + 
            theme_classic() + theme(legend.position = "bottom") +
            scale_shape_manual(values=c(21, 24))+
            scale_fill_manual(values=Cluster_colors, guide = "none")+
			labs(title = "PCA of TMM Normalized Samples", 
				x = paste0("PC1: ", round(variance_explained[1], 2), "%"), 
				y = paste0("PC2: ", round(variance_explained[2], 2), "%")
			)+ facet_wrap(~Celltype)
dev.off()



pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Ngenes_DE.pdf', height = 10)
plotter <- data.frame(cell_type=character(), method=character(), number=integer())
for(cell_type in sort(unique(metadata$cell_type))){
    tmp_wilcox <- results_wilcox[[cell_type]]
    tmp_edger  <- DE[DE$cell_type == cell_type,]
    tmp_wilcox <- tmp_wilcox[abs(tmp_wilcox$avg_log2FC) > 0.2 & tmp_wilcox$p_val_adj <0.05,]
    tmp_edger <- tmp_edger[abs(tmp_edger$avg_log2FC) > 1 & tmp_edger$p_val_adj <0.05,]
    plotter <- rbind(plotter, data.frame(cell_type=cell_type, method='Wilcox', number=nrow(tmp_wilcox)))
    plotter <- rbind(plotter, data.frame(cell_type=cell_type, method='edgeR', number=nrow(tmp_edger)))
}







########################################################################################################################
# DE Analysis Wilcox 
########################################################################################################################
# ensembl <- read.table('/home/tereshkova/data/gserranos/MDS/Data/Annotation/MDS_5qdeletion.tsv', 
# 						sep='\t', header=TRUE)
ensembl <- read.table('/home/serrang/Projects/mds5qRev/mart_export.txt', 
						sep='\t', header=TRUE)
ensembl <- ensembl[ensembl$Gene.name != '',]

sc_data_MDS5q_subseted <- readRDS('/home/serrang/Projects/mds5qRev/5qSamples_Annotated_final.rds')
Idents(sc_data_MDS5q_subseted) <- 'cell_5q'

# 'SMD34459'<-'Patient_1'
# 'SMD35109'<-'Patient_2'
# 'SMD35303'<-'Patient_3'
# 'SMD37209'<-'Patient_4'

# All samples together
DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
wilcox_normal_Vs_del5q <- FindMarkers(sc_data_MDS5q_subseted, 
                                ident.1 = 'del5q', 
                                ident.2 = 'normal',
                                logfc.threshold = 0,
                                only.pos = FALSE)

# saveRDS(wilcox_normal_Vs_del5q, '/home/tereshkova/data/gserranos/MDS/Data/wilcox_normal_Vs_del5q_allcells.rds')
wilcox_normal_Vs_del5q <- readRDS('/home/serrang/Projects/mds5qRev/wilcox_normal_Vs_del5q_allcells.rds')
wilcox_normal_Vs_del5q$gene <- rownames(wilcox_normal_Vs_del5q)
pdf('/home/serrang/Projects/mds5qRev/Plots/wilcox_normal_Vs_del5q.pdf')
get_volcano_comparison(wilcox_normal_Vs_del5q, 'wilcox_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )
dev.off()

q <- nrow(wilcox_normal_Vs_del5q[wilcox_normal_Vs_del5q$p_val_adj < 0.05 & wilcox_normal_Vs_del5q$avg_log2FC <0 & wilcox_normal_Vs_del5q$gene %in% ensembl$Gene.name, ])
m <- nrow(wilcox_normal_Vs_del5q[wilcox_normal_Vs_del5q$p_val_adj < 0.05 & wilcox_normal_Vs_del5q$avg_log2FC <0, ])
n <- nrow(wilcox_normal_Vs_del5q[wilcox_normal_Vs_del5q$p_val_adj < 0.05 & wilcox_normal_Vs_del5q$avg_log2FC >0, ])
k <- length(ensembl$Gene.name)
a <- phyper(q-1, m, n, k, lower.tail = FALSE, log.p = F)


# Per sample
# 'SMD34459'<-'Patient_1'
# 'SMD35109'<-'Patient_2'
# 'SMD35303'<-'Patient_3'
# 'SMD37209'<-'Patient_4'

Idents(sc_data_MDS5q_subseted) <- 'Sample'
DE_Wilcox_Per_Sample <- list()
for (SMP in unique(sc_data_MDS5q_subseted$Sample)){
    print(SMP)
    tmp <- subset(sc_data_MDS5q_subseted, idents=SMP)
    DefaultAssay(tmp) <- 'RNA'
    Idents(tmp) <- 'cell_5q'
    # All samples together
    DE_Wilcox_Per_Sample[[SMP]] <- FindMarkers(tmp, 
                                    ident.1 = 'del5q', 
                                    ident.2 = 'normal',
                                    logfc.threshold = 0,
                                    only.pos = FALSE)
}

for (nm in names(DE_Wilcox_Per_Sample)){
    DE_Wilcox_Per_Sample[[nm]]$gene <- rownames(DE_Wilcox_Per_Sample[[nm]])
}

pdf('/home/serrang/Projects/mds5qRev/Plots/wilcox_normal_Vs_del5q_perSample.pdf')
cowplot::plot_grid(
    get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD34459']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name ),
    get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35109']], 'SMD35109_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name ),
    get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35303']], 'SMD35303_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name ),
    get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD37209']], 'SMD37209_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name ),
ncol=2)
dev.off()


# # RANDOM
# sc_data_MDS5q_subseted <- readRDS('./5qSamples_Annotated_final.rds')
# Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
# DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
# cells_1 <- sample(colnames(sc_data_MDS5q_subseted), length(colnames(sc_data_MDS5q_subseted))/2)
# cells_2 <- colnames(sc_data_MDS5q_subseted)[!colnames(sc_data_MDS5q_subseted) %in% cells_1]
# wilcox_Random <- FindMarkers(sc_data_MDS5q_subseted, 
#                                 ident.1 = cells_1, 
#                                 ident.2 = cells_2,
#                                 logfc.threshold = 0,
#                                 only.pos = FALSE)

# # saveRDS(wilcox_Random, '/home/tereshkova/data/gserranos/MDS/Data/Revision/wilcox_Random_allcells.rds')
# # wilcox_Random <- readRDS('/home/serrang/Projects/mds5qRev/wilcox_Random_allcells.rds')

# wilcox_Random$gene <- rownames(wilcox_Random)
# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/wilcox_Random.pdf')
# get_volcano_comparison(wilcox_Random, 'wilcox_Random', genes_2_plot =ensembl$Gene.name )
# dev.off()



assign_rdm <- function(barcode){
    return( sample(c('normal', 'del5q'), 1))
}

wilcox_Random_list <- list()
for (pct in seq(10, 100, 10)){
    print(paste0('=======', pct, '======='))
    DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
    Idents(sc_data_MDS5q_subseted) <- 'cell_5q'
    all_cells <- FetchData(sc_data_MDS5q_subseted, 'cell_5q')
    print('Original cells')
    print(table(all_cells$cell_5q))
    all_cells$cell_id <- rownames(all_cells)
    # change randomly the label for X percent of the cells
    cells_2_change <- sample(all_cells$cell_id, round(nrow(all_cells)*pct/100))
    # this code assigns normal or del5q randomly based on how the barcode ends
    sapply(cells_2_change, function(cell) {
        all_cells[all_cells$cell_id == cell, 'cell_5q'] <<- assign_rdm(cell)
    })
    # this code inverts the labels, it not assign them randomly
    # all_cells$cell_5q <- ifelse(all_cells$cell_id %in% cells_2_change,
    #                             ifelse(all_cells$cell_5q=='del5q','normal','del5q'),
    #                             all_cells$cell_5q)
    cells_1 <- all_cells[all_cells == 'del5q', 'cell_id']
    cells_2 <- all_cells[all_cells == 'normal', 'cell_id']
    print('Shuffled cells')
    print(table(all_cells$cell_5q))
    wilcox_Random_list[[paste0(pct)]] <- FindMarkers(sc_data_MDS5q_subseted, 
                                ident.1 = cells_1, 
                                ident.2 = cells_2,
                                logfc.threshold = 0,
                                only.pos = FALSE)
}

plots <- lapply(wilcox_Random_list, function(x) get_volcano_comparison(x, genes_2_plot =ensembl$Gene.name ))
# pdf('./Plots/wilcox_Random.pdf', 20, 20)
# cowplot::plot_grid(plotlist=plots, labels= names(wilcox_Random_list),ncol=2)
# dev.off()


pdf('/home/serrang/Projects/mds5qRev/Plots/SupFig2.pdf', width=5, height=7)
cowplot::plot_grid(
    cowplot::plot_grid(
        get_volcano_comparison(wilcox_normal_Vs_del5q, 'wilcox_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )+ 
        labs(title='All samples', subtitle='del(5q) cells Vs non-del(5q) cells'),
        cowplot::plot_grid(
            get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD34459']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )+ labs(title='Patient 1'),
            get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35109']], 'SMD35109_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )+ labs(title='Patient 2'),
            get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35303']], 'SMD35303_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )+ labs(title='Patient 3'),
            get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD37209']], 'SMD37209_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )+ labs(title='Patient 4'),
        ncol=2),
    ncol=2, labels=c('A', 'B')),
    cowplot::plot_grid(
        plots[['30']] + xlim(-3, 1.5) + ylim(0, 600)+ labs(subtitle='30 % of the labels randomly shuffled'), 
        plots[['60']] + xlim(-3, 1.5) + ylim(0, 600)+ labs(subtitle='60 % of the labels randomly shuffled'), 
        plots[['80']] + xlim(-3, 1.5) + ylim(0, 600)+ labs(subtitle='80 % of the labels randomly shuffled'), 
        plots[['100']]+ xlim(-3, 1.5) + ylim(0, 600)+ labs(subtitle='100 % of the labels randomly shuffled'),
        ncol=2, nrow=2),
nrow=2, labels=c('A', 'C'))
dev.off()


data_per_patient <- rbind(rbind(rbind(
                        data.frame(get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD34459']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )$data, Sample='Patient 1'), 
                        data.frame(get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35109']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )$data, Sample='Patient 2')),
                        data.frame(get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD35303']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )$data, Sample='Patient 3')),
                        data.frame(get_volcano_comparison(DE_Wilcox_Per_Sample[['SMD37209']], 'SMD34459_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )$data, Sample='Patient 4'))

shuffled_data <-  rbind(rbind(rbind(
                  data.frame(plots[['30']]$data, Shuffled_percent = 30),
                  data.frame(plots[['60']]$data, Shuffled_percent = 60)),
                  data.frame(plots[['80']]$data, Shuffled_percent = 80)),
                  data.frame(plots[['100']]$data, Shuffled_percent = 100))

data_to_file <- list()
data_to_file[['FigS2-A']] <- get_volcano_comparison(wilcox_normal_Vs_del5q, 'wilcox_normal_Vs_del5q', genes_2_plot =ensembl$Gene.name )$data
data_to_file[['FigS2-B']] <- data_per_patient
data_to_file[['FigS2-C']] <- shuffled_data


WriteXLS::WriteXLS(data_to_file, ExcelFileName='/home/serrang/Projects/mds5qRev/FigS2.xlsx', 
SheetNames = names(data_to_file),  col.names=TRUE, row.names=TRUE, BoldHeaderRow=TRUE)



# results_wilcox <- readRDS('/home/tereshkova/data/gserranos/MDS/Data/Revision/DE_results_Wilcox_NoFilters.rds')
# plotter_list <- list()
# for (celltype in names(results_wilcox)){
#     plotter_list[[celltype]] <- get_volcano_comparison(results_wilcox[[celltype]],celltype, genes_2_plot =ensembl$Gene.name )
# }
# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf', height=10)
# cowplot::plot_grid(plotlist= plotter_list, ncol=4)
# dev.off()



# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf', width=15)
# Idents(sc_data_MDS5q_subseted) <- 'Cluster_names'
# VlnPlot(object = sc_data_MDS5q_subseted, features = 'FBXO38', split.by = 'cell_5q')
# dev.off()


# # Idents(sc_data_MDS5q_subseted) <- 'Sample'
# # tmp <- subset(sc_data_MDS5q_subseted, idents= 'SMD34459')
# # Idents(tmp) <- 'cell_5q'
# # tmp <- FindMarkers(tmp, 
# #                                 ident.1 = 'del5q', 
# #                                 ident.2 = 'normal',
# #                                 logfc.threshold = 0,
# #                                 only.pos = FALSE)

# # pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf', width=15)
# # get_volcano_comparison(tmp,'SMD34459', genes_2_plot =ensembl$Gene.name )
# # dev.off()



# tmp <- AggregateExpression(sc_data_MDS5q_subseted, 
#                             return.seurat = F, 
#                             slot = "counts", 
#                             assays = "RNA", 
#                             group.by = c("Cluster_names", 'cell_5q'))

# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Violins_Representatn.pdf', width=25)
# DefaultAssay(sc_data_MDS5q_subseted) <- 'RNA'
# Idents(sc_data_MDS5q_subseted) <- 'Cluster_names'
# VlnPlot(object = sc_data_MDS5q_subseted, features = c('RPS14', 'CD74', 'HINT1', 'CAMLG', 'TGFBI', 'UQCRQ'), split.by = 'cell_5q')
# dev.off()

# plotter <- reshape2::melt(tmp$RNA)
# plotter$Cell_type <- stringr::str_extract(plotter$Var2, '^[^_]+')
# plotter$Genotype  <- stringr::str_extract(plotter$Var2, '[^_]+$')
# plotter$gene <- ifelse(as.character(plotter$Var1) %in% c('RPS14', 'CD74', 'HINT1', 'CAMLG', 'TGFBI', 'UQCRQ'), 
#                         as.character(plotter$Var1), '')

# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Test.pdf', width=15)
# ggplot(plotter, aes(x= Cell_type, y= log10(value), fill=Genotype)) + 
# geom_boxplot(position = position_dodge(0.5), width=0.2) + 
# ggrepel::geom_text_repel(data = plotter[plotter$gene!= '', ], 
#                         mapping=aes(x= Cell_type, y= log10(value),label = gene), 
#                         position = position_dodge(0.5),color='black', 
#                         segment.colour = "black", fontface = 'plain', size = 2,
#                         min.segment.length = unit(0, 'lines'), 
#                         max.overlaps = 180, show.legend=FALSE, fill='white') + 
# theme_bw() + theme(legend.position='bottom', strip.background = element_blank(),
#   strip.text.x = element_blank()) +
# facet_wrap(~Cell_type, scales='free', nrow = 3)
# dev.off()



# # PSeudobulk per sample and wilcox del5q Vs normal

# aggregated_CellType_Sample_5q     <- AggregateExpression(sc_data_MDS5q_subseted, 
#                             return.seurat = T, 
#                             slot = "counts", 
#                             assays = "RNA", 
#                             group.by = c("Sample", 'cell_5q'))

# aggregated_CellType_Sample_5q$cell_5q <- stringr::str_extract(colnames(aggregated_CellType_Sample_5q), '[^_]+$')
# Idents(aggregated_CellType_Sample_5q) <- 'cell_5q'
# results <- FindMarkers(aggregated_CellType_Sample_5q,  ident.1 = 'del5q', 
#                                                 ident.2 = 'normal',
#                                                 logfc.threshold = 0,
#                                                 only.pos = FALSE)
# pdf('/home/serrang/Projects/mds5qRev/Plots/Revision/Psuedobulk_perSample_5qVsNornal.pdf')
# results$gene<- rownames(results)
# get_volcano_comparison(results, 'Pseudobulk per sample', genes_2_plot =ensembl$Gene.name )
# dev.off()
