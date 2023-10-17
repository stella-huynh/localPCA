
##  Set unique IDs to outlier regions

set_regionID <- function(tab_regions) {

  tab_regions$regID <- NA

  c <- 1
  for (n in 1:nrow(tab_regions)) {

    chr <- tab_regions[n,"chr"]
    nreg <- length(which(tab_regions$chr==chr)) #nb of outlier region on chromosome "chr"

    if(nreg==1) {

      tab_regions[n, "regID"] <- "reg0"

    } else {

      if(c==1) { c_chr <- chr
      } else if (c>1 & c_chr==chr) { c_chr <- chr
      } else if (c>1 & c_chr!=chr) { c_chr <- chr ; c <- 1
      }
      tab_regions[n, "regID"] <- paste0("reg",c)
      c <- c+1

    }

  }

  return(tab_regions)

}


##  Extract specific genomic region from VCF file

extractVCF <- function(vcf.in, chr, start, end, vcf.out, bcftools="bcftools", verbose=TRUE) {

  command_vcf <- sprintf("%s view -r %s:%i-%i -Ov -o %s %s", bcftools, chr, start, end, vcf.out, vcf.in)
  if(verbose) { print(command_vcf) }
  system(command_vcf)

}


##  Run PCA and generate plots

plotPCA <- function(vcf.fn, popFile, outDir, plots=TRUE, verbose=TRUE) {

  vcf.gds <- paste0(vcf.fn,".gds")
  SNPRelate::snpgdsVCF2GDS(vcf.fn, vcf.gds, method="biallelic.only", verbose=FALSE)
  if(verbose) { SNPRelate::snpgdsSummary(vcf.gds) }

  genofile <- SNPRelate::snpgdsOpen(vcf.gds)

  sample.id <- gdsfmt::read.gdsn(gdsfmt::index.gdsn(genofile, "sample.id"))
  pop <- read.delim(popFile, header=FALSE)
  if(ncol(pop)<3) { stop(paste("plotPCA() is expecting a 3-columns popFile (ind,pop,cols), but there only", ncol(pop), "columns is found.")) }
  colnames(pop) <- c("ind","pop","cols")

  pcaall <- SNPRelate::snpgdsPCA(genofile, autosome.only=F, remove.monosnp=F, verbose=FALSE)

  tab_pc <- data.frame(sample.id = pcaall$sample.id,
                       pop = pop$pop[match(sample.id, pop$ind)],
                       cols = pop$cols[match(sample.id, pop$ind)],
                       PC1 = pcaall$eigenvect[,1],
                       PC2 = pcaall$eigenvect[,2],
                       PC3 = pcaall$eigenvect[,3],
                       stringsAsFactors=FALSE
                       )
  #rownames(tab_pc) <- tab_pc$sample.id
  varPC1 <- round(pcaall$varprop[1]*100,2)
  varPC2 <- round(pcaall$varprop[2]*100,2)
  varPC3 <- round(pcaall$varprop[3]*100,2)

  tabfile <- paste0(outDir, basename(tools::file_path_sans_ext(vcf.fn)), "_coorpca.txt")
  write.table(tab_pc, file=tabfile, quote=F, sep="\t", row.names=F, col.names=T)

  gdsfmt::closefn.gds(genofile)


  plotfile <- paste0(outDir, basename(tools::file_path_sans_ext(vcf.fn)), "_PC1-2.pdf")
  p1 <- ggplot2::ggplot(tab_pc, aes(x=PC1, y=PC2, col=cols)) +
        ggplot2::scale_color_identity("Populations", guide="legend", labels=unique(tab_pc$pop), breaks=unique(tab_pc$cols)) +
        ggplot2::geom_point(size=8, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
        ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
                       axis.title = element_text(face="bold", size=16), axis.text = element_text(size=14),
                       legend.title = element_text(face="bold", size=16), legend.text = element_text(size=14))
  p2 <- ggplot2::ggplot(tab_pc, aes(x=PC1, y=PC2, col=cols, label=sample.id)) +
        ggplot2::scale_color_identity("Populations", guide="legend", labels=unique(tab_pc$pop), breaks=unique(tab_pc$cols)) +
        ggplot2::geom_point(size=8, alpha=0.5) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
        ggrepel::geom_text_repel(max.overlaps = 100) +
        ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
                       axis.title = element_text(face="bold", size=16), axis.text = element_text(size=14),
                       legend.title = element_text(face="bold", size=16), legend.text = element_text(size=14))
  p3 <- ggplot2::ggplot(tab_pc, aes(x=PC1, y=PC2, shape=pop, label=sample.id)) +
        ggplot2::geom_point(size=8, col="grey20", alpha=0.7) + xlab(paste0("PC1 (",varPC1," %)")) + ylab(paste0("PC2 (",varPC2," %)")) +
        ggrepel::geom_text_repel(max.overlaps = 100) +
        ggplot2::theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                       panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
                       axis.title = element_text(face="bold", size=16), axis.text = element_text(size=14),
                       legend.title = element_text(face="bold", size=16), legend.text = element_text(size=14))
  l <- mget(c("p1","p2","p3"))
  pdf(plotfile, width=10) ; invisible(lapply(l,print)) ; dev.off()
  png(sub(".pdf",".png",plotfile), height=8, width=10, units="in", pointsize=12, bg="white", res=144) ;  print(p2)  ; dev.off()


  return(list(mat = tab_pc,
              var = round(pcaall$varprop[1:3]*100,2)))
}


##  K-means clustering

kclust <- function(vcf_in, popFile, outDir, kmeans.method="euclidean", kmeans.k=3, kmeans.iter=100, plots=TRUE, verbose=TRUE) {

  pca_coor <- plotPCA(vcf.fn=vcf_in, popFile=popFile, outDir=outDir, plots=plots, verbose=verbose)

  # Get best K-means inference (ie. that lowers intraspecific difference the most) over 100 runs
  for(iter in 1:kmeans.iter) {

    #if(verbose) { cat(paste("run",iter,"\n")) }
    km_best <- factoextra::eclust(pca_coor$mat[,c("PC1","PC2")], "kmeans", hc_metric=kmeans.method, stand=TRUE, verbose=FALSE)
    wss_best <- km_best$tot.withinss

    if(kmeans.k==2 | kmeans.k==3) {
      km <- factoextra::eclust(pca_coor$mat[,c("PC1")], "kmeans", hc_metric=kmeans.method, k=kmeans.k, verbose=FALSE)
    } else if(kmeans.k>3) {
      km <- factoextra::eclust(pca_coor$mat[,c("PC1","PC2")], "kmeans", hc_metric=kmeans.method, k=kmeans.k, verbose=FALSE)
    } else {
      print("The number of clusters inferred by K-means ('kmeans.k' option) should be > 2.")
    }
    wss <- km$tot.withinss

    if(iter==1) {
      wss_best_lw <- wss_best
      wss_lw <- wss
      km_best_lw <- km_best
      km_lw <- km

    } else {
      if(wss_best < wss_best_lw) {
        wss_best_lw <- wss_best
        km_best_lw <- km_best
      }
      if(wss < wss_lw) {
        wss_lw <- wss
        km_lw <- km
      }
    }

  }

  return(list(pca_coor = pca_coor,
              km_best_lw = km_best_lw,
              km_lw = km_lw))

}



#########################################
##   5. Identify putative inversions   ##
#########################################


#' Identify putative inverted genomic regions from region-based PCA plots and k-means clustering
#'
#' @param vcf_in The VCF/BCF file used to conduct localPCA analysis (character string)
#' @param tab_regions R data.frame generated by getReg() and containing the identified outlier genomic regions
#' @param popFile Tabulated-text file containing 3 columns with the individual IDs (as in the VCF header), their assigned population (unique alphanumerical name, no space), and the Rcolor to be used for plotting, respectively. The Rcolor column is optional, only black dots will be shown on plots if missing.
#' @param outDir Output directory where to store per-region VCFs and associated PCA figures
#' @param minSize Minimum size of the outlier genomic regions (in Mb) to analyse. (default=2)
#' @param kmeans.method Distance metric used for K-means clustering (default="euclidean")
#' @param kmeans.k Number of cluster groups to be inferred by K-means. The expected number to detect inversions is 3 (or 6 in case of nested inversions), but other values are also possible. (default=3)
#' @param kmeans.iter Number of K-means iterations to run in order to find the best K-mean clustering (default=100)
#' @param plots Whether to generate plots of PCA and heterozygosity boxplots for each region (default=TRUE)
#' @param keep_tmp Whether to keep the per-region VCF files generated (may take up some space if large genome) (default=FALSE)
#' @param vcftools By default, vcftools software is installed in /usr/bin. If this is not the case (no root permission),use this option to specify your user-defined full path to vcftools
#' @param bcftools By default, bcftools software is installed in /usr/bin. If this is not the case (no root permission), use this option to specify your user-defined full path to bcftools
#' @param bgzip By default, bgzip software is installed in /usr/bin. If this is not the case (no root permission), use this option to specify your user-defined full path to bgzip Required if keep_tmp=TRUE
#' @param verbose Whether to display information on the analysis (default=TRUE)

detectInv <- function(vcf_in, tab_regions=regList$tabF, popFile, outDir="localPCA", minSize=2, kmeans.method="euclidean", kmeans.k=3, kmeans.iter=100, keep_tmp=FALSE, plots=TRUE, bcftools="bcftools", vcftools="vcftools", bgzip="bgzip", verbose=TRUE) {

  if(!dir.exists(outDir)) {
    dir.create(path=outDir, showWarnings=FALSE, recursive=TRUE)
    if(verbose) { cat(paste0("\n  ---  Created folder \"", outDir, "\" to store MDS with outliers plots  ---  \n")) }
  }
  if(!endsWith(outDir,"/")) { outDir <- paste0(outDir,"/") }

  ####################################################
  ##   1. Set table & VCF for each outlier region   ##
  ####################################################

  prefix <- basename(sub(".bcf$|.vcf.gz$","",vcf_in))

  tab_vcf <- set_regionID(tab_regions) %>%
               dplyr::filter_at("SizeRegions", dplyr::all_vars(. >= minSize)) %>%
               dplyr::mutate(vcf_reg = NA)

  list_PCAclust <- list()
  list_res <- list()

  for (n in 1:nrow(tab_vcf)) {

    if(verbose) { cat(paste0("-----------------------------\nAnalysing region number: ",n,"\n-----------------------------\n\n")) }

    ## Generate regional VCF file
    chr <- tab_vcf[n,"chr"]
    region <- tab_vcf[n, "regID"]
    start <- tab_vcf[n,"start"]
    end <- tab_vcf[n,"end"]
    vcf_reg <- paste0(outDir, prefix, ".localPCA.", chr, ".", region, ".vcf")
    tab_vcf[n,"vcf_reg"] <- vcf_reg

    extractVCF(vcf.in=vcf_in, chr=chr, start=start, end=end, vcf.out=vcf_reg, bcftools=bcftools, verbose=verbose)

    ## Generate regional PCA k-means plots
    list_PCAclust[[n]] <- kclust(vcf_in=vcf_reg, popFile=pop_file, outDir=outDir, kmeans.method=kmeans.method, kmeans.k=kmeans.k, kmeans.iter=kmeans.iter, plots=plots, verbose=verbose)

    pca_coor <- list_PCAclust[[n]]$pca_coor ; rownames(pca_coor$mat) = NULL
    km_lw <- list_PCAclust[[n]]$km_lw

    hom1 <- pca_coor$mat[which(km_lw$cluster==which(km_lw$centers[,1]==min(km_lw$centers[,1]))),"sample.id"]
    hom2 <- pca_coor$mat[which(km_lw$cluster==which(km_lw$centers[,1]==max(km_lw$centers[,1]))),"sample.id"]
    hz <- pca_coor$mat[which(!pca_coor$mat[,"sample.id"] %in% c(hom1,hom2)),"sample.id"]
    if(verbose) {
      cat(paste0("\nFor file \"", vcf_reg, "\", the three groups are:\n\n  - Group \"hom1\":\n       ", paste0(hom1,collapse=","),
                 "\n\n  - Group \"hom2\":\n       ", paste0(hom2,collapse=","),
                 "\n\n  - Group \"hz\":\n       ", paste0(hz,collapse=","),"\n\n"))
    }

    pca_coor$mat[,"cluster"] <- as.factor(unlist(lapply(pca_coor$mat["sample.id"], function(x) ifelse(x %in% hom1, 1, ifelse(x %in% hz, 2, ifelse(x %in% hom2, 3, NA))))))
    list_PCAclust[[n]] <- append(list_PCAclust[[n]], list(hom1=hom1, hom2=hom2, hz=hz))

    ## Generate separate VCF file for each group (for both outlier/neutral regions)
    df_het <- NULL
    cols_het <- as.data.frame(cbind(group=c("hom1","hz","hom2"), cols=c("darkorange2","seagreen","darkgoldenrod1")))

    for(group in c("hom1","hom2","hz")) {
      het_out <- sub(".vcf$|vcf.gz$", paste0(".", group, ".het"), vcf_reg)
      het_log <- paste0(het_out,".log")
      command_het <- sprintf("%s view -s %s %s | %s --vcf - --out %s --het >> %s 2>&1", bcftools, paste0(get(group),collapse=","), vcf_reg, vcftools, het_out, het_log)
      cat(paste0("Command:\n\n",command_het,"\n\n\nVCFTOOLS LOG output:\n"), file=het_log)
      system(command_het)
      system(sprintf("mv %s %s", paste0(het_out,".het"), het_out))
      df_het <- rbind(df_het, cbind(read.delim(het_out),group=group))
    }

    df_het$O.HET <- (df_het$N_SITES - df_het$O.HOM) / df_het$N_SITES
    df_het$E.HET <- (df_het$N_SITES - df_het$E.HOM) / df_het$N_SITES
    df_het$pop <- pca_coor$mat[match(df_het$INDV,pca_coor$mat$sample.id),"pop"]
    df_het$cols <- cols_het[match(df_het$group,cols_het$group),"cols"]
    df_het$group <- factor(df_het$group, levels=c("hom1","hz","hom2"))

    mat_res <- base::merge(df_het, pca_coor$mat, by.x="INDV", by.y="sample.id")
    mat_res <- mat_res[,c("INDV","pop.x","cols.y","group","cluster","cols.x","PC1","PC2","PC3","N_SITES","O.HET","E.HET")]
    colnames(mat_res) <- c("sample.id","pop","pop.cols","group","cluster","group.cols","PC1","PC2","PC3","het.nsites","het.obs","het.exp")
    #m_het <- reshape2::melt(mat_res, id.vars=c("sample.id","pop","pop.cols","group","group.cols"))
    mat_res$group.cols <- as.factor(mat_res$group.cols)
    mat_res$group.cols <- factor(mat_res$group.cols, levels=mat_res[match(c("hom1","hz","hom2"),mat_res$group),"group.cols"])

    #cmp_h1z = wilcox.test(mat_res[which(mat_res$group=="hom1"),"het.obs"],mat_res[which(mat_res$group=="hz"),"het.obs"])$p.value
    #cmp_h2z = wilcox.test(mat_res[which(mat_res$group=="hom2"),"het.obs"],mat_res[which(mat_res$group=="hz"),"het.obs"])$p.value
    #cmp_h12 = wilcox.test(mat_res[which(mat_res$group=="hom1"),"het.obs"],mat_res[which(mat_res$group=="hom2"),"het.obs"])$p.value
    #ggpubr::compare_means(het.obs ~ group, data=mat_res, method="wilcox.test", label="p.signif", p.adjust.method="bonferroni")
    pval_bonfc_05 <- 0.05 ^ ncol(combn(kmeans.k,2))
    pval_bonfc_01 <- 0.01 ^ ncol(combn(kmeans.k,2))
    pval_bonfc_001 <- 0.001 ^ ncol(combn(kmeans.k,2))
    pval_text <- paste0("Significance after Bonferroni correction for pairwise comparisons between ", kmeans.k, " samples: * < ", pval_bonfc_05, ", ** < ", pval_bonfc_01, ", *** < ", pval_bonfc_001)
    stat.test <- mat_res %>%
                   rstatix::wilcox_test(het.obs ~ group, p.adjust.method = "bonferroni") %>%
                   rstatix::add_xy_position(x="group")
    stat.test$p.adj.bonferroni <- ifelse(stat.test$p < pval_bonfc_001, "***",
                                         ifelse(stat.test$p < pval_bonfc_01, "**",
                                                ifelse(stat.test$p < pval_bonfc_05, "*", "ns"
                                         )))

    # Generate plots
    km_best_lw <- list_PCAclust[[n]]$km_best_lw
    p0 <- km_best_lw$clust_plot +
          ggtitle("Best K-means clustering plot") +
          theme(plot.title = element_text(size=20, face="bold"), axis.title = element_text(size=20, face="bold"),
                legend.title = element_text(size=20, face="bold"), legend.text = element_text(size=18))

    varPC1 <- pca_coor$var[1]
    varPC2 <- pca_coor$var[2]
    p1 <- ggplot(mat_res, aes(x=PC1, y=PC2, col=group.cols, label=sample.id)) +
          ggtitle("PCA based on SNPs present in candidate structural variant") +
          geom_point(size=8, alpha=0.5) +
          xlab(paste0("PC1 (",varPC1," %)")) +
          ylab(paste0("PC2 (",varPC2," %)")) +
          labs(color="Groups") + scale_colour_identity(guide="legend", labels=levels(mat_res$group)) +
          ggrepel::geom_text_repel(max.overlaps = 100) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = 'white'), axis.line=element_line(colour="black"),
                axis.title = element_text(face="bold", size=20), axis.text = element_text(size=18),
                legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=18),
                plot.title = element_text(size=20, face="bold"))

    p2 <- ggplot(mat_res, aes(x=group, y=het.obs, fill=group.cols)) +
          ggtitle("Heterozygosity between the Kmeans-based groups") +
          geom_boxplot() +
          scale_fill_identity("Groups",guide="legend", labels=levels(mat_res$group)) +
          ylab("Proportion of heterozygous sites") +
          scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.2)) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = 'white'), axis.line = element_line(colour="black"),
                axis.text = element_text(size=18), axis.title.y = element_text(size=20, face="bold"), axis.title.x = element_blank(),
                legend.title = element_text(size=18, face="bold"), legend.text = element_text(size=18), legend.position = "right",
                plot.title = element_text(size=20, face="bold")) +
          ggpubr::stat_compare_means(comparisons = list(c("hom1","hz"),c("hom2","hz"),c("hom1","hom2")),
                                     method="wilcox.test",
                                     p.adjust.method="bonferroni",
                                     symnum.args=list(cutpoints=c(0,pval_bonfc_001,pval_bonfc_01,pval_bonfc_05,Inf),
                                                      symbols=c("***","**","*","ns")))

    plotfile <- paste0(sub(".vcf",".boxplot.pdf",vcf_reg))
    l <- mget(c("p0","p1","p2"))
    pdf(plotfile, height=8, width=16) ;  invisible(lapply(l,print))  ; dev.off()
    png(sub(".pdf",".0.bestK.png",plotfile), height=8, width=16, units="in", pointsize=12, bg="white", res=144) ;  print(p0)  ; dev.off()
    png(sub(".pdf",".1.PCA.clustK3.png",plotfile), height=8, width=16, units="in", pointsize=12, bg="white", res=144) ;  print(p1)  ; dev.off()
    png(sub(".pdf",".2.Ho.png",plotfile), height=8, width=16, units="in", pointsize=12, bg="white", res=144) ;  print(p2)  ; dev.off()

    # Summary statistics
    if(verbose) {
      cat("\nObserved heterozygosity: quantile distribution\n")
      print(df_het %>%
              dplyr::group_by(group) %>%
              dplyr::summarize(min = min(O.HET),
                               q1 = quantile(O.HET, 0.25),
                               median = median(O.HET),
                               mean = mean(O.HET),
                               q3 = quantile(O.HET, 0.75),
                               max = max(O.HET))
      )
      cat("\n\n")
    }

    list_res[[n]] <- list(mat = mat_res,
                          pca.var = pca_coor$var[1:3],
                          het.wilcox.test = stat.test,
                          het.bonf.corr = list("*"=pval_bonfc_05, "**"=pval_bonfc_01, "***"=pval_bonfc_001))

  }

  return(list_res)

}


