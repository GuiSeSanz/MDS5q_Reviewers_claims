library(ggplot2)

clean_file <- function(sample_name){
	cmd <- paste0("grep '^5\t' /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/", sample_name, "_norm.bam.mpileup.GermLine_snp > /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/", sample_name, "_norm.bam.mpileup.GermLine_snp_clean")
	system(cmd)
	cmd <- paste0("grep '^5\t' /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/", sample_name, "_del5q.bam.mpileup.GermLine_snp > /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/", sample_name, "_del5q.bam.mpileup.GermLine_snp_clean")
	system(cmd)
}

plot_vars <- function(data, threshold=0, smp_name = SAMPLE){
	tmp <- data[data$Pos1 + data$Pos2 >threshold, ]
	p <- ggplot(tmp, aes(x=Pos, y=Var, color=Geno ))+ 
		scale_color_manual(values=c('norm' = "#999999", 'del5' = "#E69F00"))+
		geom_point(alpha=0.8, size = 0.8) + theme_classic() + 
		guides(color = guide_legend(override.aes = list(size=5, alpha = 1))) +
		facet_wrap(~Geno, nrow=2) + xlim(min(tmp$Pos), max(tmp$Pos)) + 
		geom_rect(data=rects, inherit.aes=FALSE, aes(xmin=start, xmax=end, ymin=min(tmp$Var),
						ymax=max(tmp$Var), group=group), color="transparent", fill="steelblue", alpha=0.2) + 
		theme(legend.position='none')
	if(threshold>0){
		p <- p+ labs(y=paste0('Sum of Reads > ', threshold)) + theme(axis.title.x=element_blank())
	}else{
		p <- p+ labs(subtitle=paste0('All Data for sample ', smp_name))
	}
	return(p)
}

plot_list <- list()
plot_list_desity <- list()
plot_GermLine <- list()
for (SAMPLE in c('SMD37209', 'SMD35303', 'SMD34459', 'SMD35109')){
	print(SAMPLE)
	# Keep only the coordinates, no header no tail
	# system("grep '^5' /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/SMD37209_norm.bam.mpileup.GermLine_snp > /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/SMD37209_norm.bam.mpileup.GermLine_snp_clean")
	# system("grep '^5' /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/SMD37209_del5q.bam.mpileup.GermLine_snp > /home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/SMD37209_del5q.bam.mpileup.GermLine_snp_clean")
	clean_file(SAMPLE)
	norm_data <- read.table(paste0('/home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/',SAMPLE,'_norm.bam.mpileup.GermLine_snp_clean'), sep='\t')
	del5_data <- read.table(paste0('/home/serrang/Projects/mds5qRev/Rev/BAM_Bogaloo/',SAMPLE,'_del5q.bam.mpileup.GermLine_snp_clean'), sep='\t')

	norm_data <- setNames(norm_data[, c(1, 2, 5, 6, 7)], c('Chr', 'Pos', 'Pos1', 'Pos2', 'Var'))
	del5_data <- setNames(del5_data[, c(1, 2, 5, 6, 7)], c('Chr', 'Pos', 'Pos1', 'Pos2', 'Var'))
	norm_data$Var <- as.numeric(sub('%', '', norm_data$Var))
	del5_data$Var <- as.numeric(sub('%', '', del5_data$Var))
	norm_data$Geno <- 'norm'
	del5_data$Geno <- 'del5'

	all_data<-rbind(norm_data, del5_data)
	all_data$Var <- as.numeric(sub('%', '', all_data$Var))
	# 74084068       q13.3
	# 145005789      q31.3
	if(SAMPLE %in% c('SMD35303', 'SMD34459', 'SMD35109')){
		start <- 74258689
		end <- 149430733
	}else if(SAMPLE == 'SMD37209'){
		start <- 74258689
		end <- 168616996
	}

	all_data$locus <- ifelse(all_data$Pos <= start, 'Pre', 
				  ifelse(all_data$Pos>= start & all_data$Pos <= end, 'Target', 'Post'))
	all_data$locus <- factor(all_data$locus , levels=c('Pre', 'Target', 'Post'))
	rects <- data.frame(start=start, end=end, group='5q')

	pdf(paste0('/home/serrang/Projects/mds5qRev/Plots/VarsCan_',SAMPLE,'.pdf'))
	print(cowplot::plot_grid(plot_vars(all_data), 
			cowplot::plot_grid(
				plot_vars(all_data, 8),
				plot_vars(all_data, 20),
				plot_vars(all_data, 50),
				plot_vars(all_data, 100),
				ncol=2, nrow=2),
		nrow=2, rel_heights=c(0.4,0.6)))
	print(cowplot::plot_grid(plot_vars(all_data[all_data$Pos %in% norm_data$Pos,]), 
			cowplot::plot_grid(
				plot_vars(all_data[all_data$Pos %in% norm_data$Pos,], 8),
				plot_vars(all_data[all_data$Pos %in% norm_data$Pos,], 20),
				plot_vars(all_data[all_data$Pos %in% norm_data$Pos,], 50),
				plot_vars(all_data[all_data$Pos %in% norm_data$Pos,], 100),
				ncol=2, nrow=2),
	nrow=2, rel_heights=c(0.4,0.6)))
	dev.off()
	tmp <- all_data[all_data$Pos %in% norm_data$Pos,]
	tmp <- tmp[tmp$Pos1 + tmp$Pos2 >20, ]
	pval <- ks.test(tmp[tmp$locus =='Target' & tmp$Geno=='del5', 'Var'], tmp[tmp$locus =='Target' & tmp$Geno=='norm', 'Var'])$p.value
	if(pval == 0){
		pval <- 'p-value < 2.2e-16'
	}
	plot <- plot_vars(all_data[all_data$Pos %in% norm_data$Pos,], 8) + 
	labs(title=paste0('Sample ', SAMPLE), subtitle=paste0(' Ks_pval: ', pval))
	plot_list[[SAMPLE]] <- plot

	pos_2_plot_1 <- all_data[all_data$Geno == 'norm' & all_data$Pos1 + all_data$Pos2 >=8 , ]
	pos_2_plot_1 <- pos_2_plot_1[pos_2_plot_1$Var >= 40 & pos_2_plot_1$Var <= 60, ]

	pos_2_plot_2 <- all_data[all_data$Geno == 'norm' & all_data$Pos1 + all_data$Pos2 >=8, ]
	pos_2_plot_2 <- pos_2_plot_2[pos_2_plot_2$Var >= 80 & pos_2_plot_2$Var <= 100, ]
	pos_2_plot <- rbind(pos_2_plot_1, pos_2_plot_2)


	pos_2_plot_3 <- all_data[all_data$Pos %in% pos_2_plot$Pos & all_data$Geno == 'del5',]
	pos_2_plot <- rbind(pos_2_plot, pos_2_plot_3)

	pval <- ks.test(pos_2_plot_3$Var, pos_2_plot_2$Var)$p.value
	if(pval == 0){
		pval <- 'p-value < 2.2e-16'
	}
	plot_GermLine[[SAMPLE]] <- plot_vars(pos_2_plot, 8) + labs(title=paste0('Sample ', SAMPLE), subtitle=paste0(' Ks_pval: ', pval))

	plot <- ggplot(pos_2_plot, aes(x=Var, color=Geno ,  fill=Geno))+ 
	scale_color_manual(values=c('norm' = "#999999", 'del5' = "#E69F00"))+
	scale_fill_manual(values=c('norm' = "#999999", 'del5' = "#E69F00"))+
	theme_classic() + stat_density(adjust = 1/8)+ theme(legend.position='none')+
	xlim(0, 100) + facet_wrap(~Geno, nrow=2)  + 
	labs(title=paste0('Sample ', SAMPLE), subtitle=paste0(' Ks_pval: ', pval))
	plot_list_desity[[SAMPLE]] <- plot
}

pdf('/home/serrang/Projects/mds5qRev/Plots/plot_GermLine.pdf')
cowplot::plot_grid(plotlist=plot_GermLine, ncol=2)
dev.off()

pdf(paste0('/home/serrang/Projects/mds5qRev/Plots/VarsCan_All_samples_Ks.pdf'))
cowplot::plot_grid(plotlist=plot_list, ncol=2)
dev.off()


pdf('/home/serrang/Projects/mds5qRev/Plots/Hist_Varscan.pdf')
cowplot::plot_grid(plotlist=plot_list_desity, ncol=2)
dev.off()



ks.test(pos_2_plot_3$Var, pos_2_plot_2$Var)

# pdf(paste0('/home/serrang/Projects/mds5qRev/Plots/Hist.pdf'))
# tmp <- all_data[all_data$Pos %in% norm_data$Pos,]
# tmp <- tmp[tmp$Pos1 + tmp$Pos2 >20, ]
# ggplot(tmp, aes(x= locus,y=Var, fill=Geno))+geom_point(alpha=0.3) 
# # , scales='free')
# dev.off()
tmp <- all_data[all_data$Pos %in% norm_data$Pos,]
tmp <- tmp[tmp$Pos1 + tmp$Pos2 >20, ]
ks.test(tmp[tmp$locus =='Pre'    & tmp$Geno=='del5', 'Var'], tmp[tmp$locus =='Pre'    & tmp$Geno=='norm', 'Var'])
ks.test(tmp[tmp$locus =='Target' & tmp$Geno=='del5', 'Var'], tmp[tmp$locus =='Target' & tmp$Geno=='norm', 'Var'])
ks.test(tmp[tmp$locus =='Post'   & tmp$Geno=='del5', 'Var'], tmp[tmp$locus =='Post'   & tmp$Geno=='norm', 'Var'])







