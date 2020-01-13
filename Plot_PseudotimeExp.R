###source('/users/jwang3/RetReg/Programs/Plot_PseudotimeExp.R')

Plot_PseudotimeExp <- function(cds1, File1, PlotType1=NULL, Gene1=NULL, Spec1='', Group1='', Group2='', show_backbone=T, RedMeth1='Component', Legend1='top', nClusters=1, width=12, height=10, axis.text.size=24, axis.title.size=24, plot.title.size=10, Col1=NULL, facet_wrap_nrow=3, logT1=T, show_rownames=F, use_gene_short_name=T, length.out=5, cell_size=0.75, Original.gene.order=T){
  library(VGAM); library(monocle)
  library(pheatmap)
    
  if(is.null(Col1)){ Col2 <- c(rgb(152/255,152/255,152/255), rgb(255/255,191/255,15/255), rgb(247/255,147/255,30/255), rgb(220/255,20/255,60/255), rgb(70/255,157/255,214/255), rgb(103/255,199/255,193/255), rgb(140/255,198/255,63/255), rgb(206/255,121/255,178/255), rgb(143/255,72/255,156/255), rgb(157/255,115/255,194/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(150/255,206/255,180/255), rgb(163/255,165/255,0/255), rgb(82/255,24/255,74/255), rgb(129/255,70/255,58/255))
    if(Group1=='protocol'){ Col2 <- c(rgb(169/255,169/255,169/255), rgb(103/255,199/255,193/255), rgb(3/255,161/255,198/255), rgb(97/255,156/255,255/255), rgb(0/255,114/255,189/255), rgb(140/255,198/255,63/255), rgb(192/255,193/255,48/255), rgb(247/255,147/255,30/255), rgb(220/255,20/255,60/255), rgb(230/255,134/255,201/255), rgb(157/255,115/255,194/255), rgb(143/255,72/255,156/255), rgb(163/255,165/255,0/255), rgb(129/255,70/255,58/255), rgb(237/255,221/255,82/255), rgb(255/255,195/255,0/255), rgb(255/255,141/255,26/255), rgb(255/255,87/255,51/255), rgb(199/255,0/255,57/255), rgb(144/255,12/255,62/255), rgb(82/255,24/255,74/255), rgb(60/255,61/255,107/255), rgb(42/255,122/255,155/255), rgb(0/255,186/255,173/255), rgb(86/255,199/255,133/255), rgb(173/255,212/255,93/255)) }
  }else{ Col2 <- Col1 }
  
  if(!is.null(Gene1)){
    source('/users/jwang3/RetReg/Programs/Converse_GeneIDSymbol.R')
    GeneInf1 <- read.table(paste0('/users/jwang3/RetReg/Public/',Spec1,'scRNA_genes.txt'), row.names=1, sep='\t')
    Gene2 <- Converse_GeneIDSymbol(Gene1, GeneInf1)
  }
  
  if(is.element('Trajectory', PlotType1)){
    uGroup1 <- unique(pData(cds1)[, Group1])
    if(length(uGroup1)==3){ Col2 <- Col2[c(1,3,5)] }
    pdf(file=paste0(File1,'_',Group1,'Traj.pdf'), width=width, height=height)
    print(plot_cell_trajectory(cds1, cell_size=1.3, cell_link_size=1, color_by = Group1, show_backbone=show_backbone) + scale_color_manual(values=Col2[1:length(uGroup1)]) + labs(x=paste(RedMeth1,'1'), y=paste(RedMeth1,'2')) + theme(legend.text=element_text(size=12), axis.line=element_line(size=1), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size), panel.border = element_blank(), legend.position=Legend1) )
    dev.off()
  }
  
  if(is.element('Trajectory_sample', PlotType1)){
    print('Plot trajectory for each sample')
    if(Group1==''){ stop('No Group1 information')
    }else if(Group2==''){ Group2 <- Group1 }
    
    pdf(file=paste0(File1,'_',Group1,'Traj_Sample.pdf'), width=width, height=height)
    print(plot_cell_trajectory(cds1, cell_size=0.6, cell_link_size=1, color_by = Group1) + facet_wrap(paste0('~',Group2), nrow = facet_wrap_nrow) + scale_color_manual(values=Col2[1:length(unique(pData(cds1)[, Group1]))]) + labs(x=paste(RedMeth1,'1'), y=paste(RedMeth1,'2')) + theme(legend.position = Legend1, strip.text=element_text(size=plot.title.size), axis.line=element_line(size=1), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size)))
    dev.off()
  }

  if(is.element('Trajectory_expression', PlotType1)){
    source('/users/jwang3/RetReg/Programs/Plot_MarkersExp.R')
    
    if(Gene1[1]=='' & is.element(Group1, colnames(pData(cds1)))){
      if(!is.numeric(pData(cds1)[1, Group1])){
        cds11 <- as.data.frame(cbind(t(cds1@reducedDimS), as.factor(pData(cds1)[, Group1])))
      }
      colnames(cds11) <- c(paste(RedMeth1,1:3), Group1)
      LengGroup1 <- length(unique(cds11[, Group1]))
    }else{
      if(logT1==T){ cds11 <- as.data.frame(t(rbind(cds1@reducedDimS, log2(exprs(cds1)[Gene2[, 'EnsemblID'], ]+1))))
      }else{ cds11 <- as.data.frame(t(rbind(cds1@reducedDimS, exprs(cds1)[Gene2[, 'EnsemblID'], ]))) }
      if(nrow(cds1@reducedDimS)==2){ colnames(cds11) <- c('Componet1', 'Componet2', Gene2[, 'Symbol'])
      }else if(nrow(cds1@reducedDimS)==3){ colnames(cds11) <- c('Componet1', 'Componet2', 'Componet3' , Gene2[, 'Symbol']) }
    }

    for(i in (nrow(cds1@reducedDimS)+1):ncol(cds11)){
      GeneSymb1 <- as.character(colnames(cds11)[i]); print(GeneSymb1)
      cds2 <- cds11[, c(1:(nrow(cds1@reducedDimS)), i)]
      ExpGrade1 <- seq(min(cds2[, nrow(cds1@reducedDimS)+1]), max(cds2[, nrow(cds1@reducedDimS)+1]), length.out=length.out)
      cds3 <- cbind(cds2[cds2[, nrow(cds1@reducedDimS)+1]==ExpGrade1[1], ], 1)
      colnames(cds3)[(ncol(cds3)-1):ncol(cds3)] <- c(GeneSymb1,'Group')
      for(i in 2:length(ExpGrade1)){
        cds21 <- cds2[cds2[, nrow(cds1@reducedDimS)+1]>ExpGrade1[i-1] & cds2[, nrow(cds1@reducedDimS)+1]<=ExpGrade1[i], ]
        if(nrow(cds21)>0){ 
          cds22 <- cbind(cds21, 2)
          colnames(cds22)[(ncol(cds22)-1):ncol(cds22)] <- c(GeneSymb1,'Group')
          cds3 <- rbind(cds3, cds22); 
      } }

      pdf(file=paste0(File1,'_traj_',GeneSymb1,'.pdf'), width=width, height=height)
      if(is.numeric(cds3[1,ncol(cds3)])){
        print(ggplot(cds3, aes(Componet1, Componet2)) + geom_point(aes(colour = cds3[, nrow(cds1@reducedDimS)+1], group=Group), size=0.8) + scale_colour_gradient(name=GeneSymb1, low = rgb(169/255,169/255,169/255), high = rgb(220/255,20/255,60/255)) + labs(x=paste(RedMeth1,'1'), y=paste(RedMeth1,'2')) + theme(legend.text=element_text(size=12), axis.line=element_line(size=0.75), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size), panel.border = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) )
      }else{
        print(ggplot(cds3, aes(Componet1, Componet2)) + geom_point(aes(colour = cds3[, nrow(cds1@reducedDimS)+1], group=Group), size=0.8) + scale_color_manual(values=Col2[1:LengGroup1]) + labs(x=paste(RedMeth1,'1'), y=paste(RedMeth1,'2')) + theme(legend.text=element_text(size=12), axis.line=element_line(size=0.75), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size), panel.border = element_blank()) + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank()) )
      }
      dev.off()
  } }
      
  if(is.element('Pseudotime', PlotType1)){
    pdf(file=paste0(File1,'_pseudotime.pdf'), width=width, height=height)
    print(plot_cell_trajectory(cds1, cell_size=1.3, cell_link_size=1, show_backbone=show_backbone) + labs(x=paste(RedMeth1,'1'), y=paste(RedMeth1,'2')) + theme(legend.text=element_text(size=12), axis.line=element_line(size=1), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size), panel.border = element_blank(), legend.position=Legend1) )
    dev.off() 
  }
  
  if(is.element('Pseudotime_sample_density', PlotType1)){
    print('Plot sample density based on pseudotime')
    pdf(file=paste0(File1,'_pseudotime_sample_density.pdf'), width=width, height=height)
    print(qplot(Pseudotime, data = pData(cds1), color = pData(cds1)[, Group1], geom ="density"))
    dev.off() 
  }
  
  if(is.element('Pseudotime_expression', PlotType1)){
    print('plot pseudotime expression of genes')
    cds1_subset <- cds1[Gene2[,'EnsemblID'], ]

    if(Original.gene.order==T){
      fData(cds1_subset)[, 'gene_short_name'] <- paste0(1:nrow(Gene2),'.',fData(cds1_subset)[, 'gene_short_name'])
    }
    pdf(file=paste0(File1,'_Pseudotime_expression.pdf'), width=width, height=height)
    print(plot_genes_in_pseudotime(cds1_subset, color_by = Group1, cell_size=cell_size) + scale_color_manual(values=Col2[1:length(unique(pData(cds1)[, Group1]))]) + theme(legend.position = Legend1, strip.text=element_text(size=plot.title.size), axis.line=element_line(size=1), axis.text = element_text(size=axis.text.size), axis.title = element_text(size=axis.title.size)))
    dev.off()
  }

  if(is.element('Pseudotime_heatmap', PlotType1)){
    print('Clustering Genes by Pseudotemporal Expression Pattern')
    fData(cds1)[, 'gene_short_name'] <- make.names(fData(cds1)[, 'gene_short_name'], unique=T)
    pdf(file=paste0(File1,'_G',nClusters,'_traj_heatmap.pdf'), width=width, height=height)    
    if(!is.null(Col1)){
      p1 <- plot_pseudotime_heatmap(cds1[Gene1, ], num_clusters = nClusters, cores = 3, show_rownames = show_rownames, use_gene_short_name=use_gene_short_name, return_heatmap=T, hmcols=Col1)
    }else{ p1 <- plot_pseudotime_heatmap(cds1[Gene1, ], num_clusters = nClusters, cores = 3, show_rownames = show_rownames, use_gene_short_name=use_gene_short_name, return_heatmap=T) }
    dev.off()
    
    Group1 <- cutree(p1$tree_row, nClusters);
    Gene2 <- cbind(Gene1, Group1); colnames(Gene2) <- c('ID','Group')
    Group2 <- c()
    for(i in 1:nClusters){
      Group2 <- c(Group2, rownames(Gene2[Gene2[,'Group']==i,])[1:2])
    }
    
    cds11 <- cds1[Group2, ]
    pdf(file=paste0(File1,'_G',nClusters,'_traj_heatmap_part.pdf'), width=width, height=height)
    plot_pseudotime_heatmap(cds11, cluster_rows=F, num_clusters = 1, cores = 3, show_rownames = T, return_heatmap=F)
    dev.off()
    
    return(Gene2)
  }
      
  if(is.element('Branched_pseudotime', PlotType1)){
    print('Analyzing Branches in Single-Cell Trajectories')
    BEAM_res <- BEAM(cds1, branch_point = 1, cores = 1)
    BEAM_res <- BEAM_res[order(BEAM_res$qval),]
    BEAM_res <- BEAM_res[,c("gene_short_name", "pval", "qval")]
    pdf(file=paste0(File1,'_traj_branched_heatmap.pdf'), width=width, height=height)
    plot_genes_branched_heatmap(cds1[row.names(subset(BEAM_res, qval < 0.01)),], branch_point = 1, num_clusters = 5, cores = 3, use_gene_short_name = T, show_rownames = F)
    dev.off()
    
    print('plot branched pseudotime expression')
    cds1_genes <- row.names(subset(fData(cds1), gene_short_name %in% c("Ccnd2", "Sftpb", "Pdpn")))
    plot_genes_branched_pseudotime(cds1[lung_genes,], branch_point = 1, color_by = "Time", ncol = 1)
  }
  
  if(is.element('Cell_clusters', PlotType1)){
    print('Plot cluster of single cells')
    pdf(file=paste0(File1,'_',Group1,'CellClusters.pdf'), width=width, height=height)
    print(plot_cell_clusters(cds1, color_by = Group1, cell_size = 0.5, show_group_id = T))
    dev.off()
  }
  
}
