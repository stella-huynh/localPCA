####################################
##   O. Main local PCA analysis   ##
####################################


#' Local PCA analysis based on Li and Ralph (2019)
#'
#' @param vcf A VCF/BCF file containing populational genome-wide SNP data
#' @param pos_file A 4-columns data.frame reporting the SNP number, the chromosome name, and the starting and end positions of the SNP-based genomic windows on the given chromosome. Column names are c"(ID", "Chr", "posStart", "posEnd")
#' @param outDir Output directory where to store local PCA results and figures
#' @param winsize Size of non-overlapping genomic windows for performing localPCA analysis. (default=1000)
#' @param wintype Units defining the sizes of genomic windows (option="snps","kb"). (default="snps")
#' @param k Number of principal components (PCs) to use for calculating the patterns of relatedness between each genomic window (default=2)
#' @param nMDS Number of MDS dimension to calculate (default=40)
#' @param plots Whether to plot MDS coordinates along chromosomes (default=TRUE)
#' @param verbose Whether to display information on the analysis (default=TRUE)

#' @return Returns a list of 7 elements:
#'         \tabular{ll}{
#'           winfun \tab \pkg{lostruct} win-function to divide the each chromosome into windows (based on bcftools)\cr
#'           pcs_all \tab data.frame containing the patterns of relatedness between individuals based on the first k PCs\cr
#'           rows_na \tab row indices of the data.frame "pcs_all" containing NA\cr
#'           pcs \tab data.frame "pcs_all" but excluding rows containing NA\cr
#'           pcdist \tab dissimilarity matrix between all genomic windows\cr
#'           fit2d \tab list of summarized dissimilarity matrix into nMDS MDS dimensions\cr
#'           data \tab data.frame containing the results of local PCA analysis\cr
#'         }
#'
#' @export loPCA
#' @import lostruct
#' @import ggplot2
#' @importFrom scales unit_format
#' @importFrom utils read.delim
#'


#loPCA_class = setClass(Class="loPCA",
#                       representation(snps="function", pcs="matrix", pcdist="matrix", fit2d="matrix", df="data.frame"),
#                       package="localPCA")

loPCA <- function(vcf, pos_file, outDir="localPCA", winsize=1000, wintype="snp", k=2, nMDS=40, plots=TRUE, verbose=TRUE) {

  ##   !!!!   Requires bcftools to be installed in /usr/bin/   !!!!


  #  ** A-1. Divide genome into windows **

  options(warn=1) #prevent Rscript to exit when vcf_windower() has to ignore too-short windows on chromosome's ends.
  if(verbose) { cat("\n#  ***  A- 1. Divide genome into windows  ***  #\n") }

  snps <- lostruct::vcf_windower(vcf, size=winsize, type=wintype)


  #  ** A-2. Summarize patterns of relatedness in each window **

  if(verbose) { cat("\n#  ***  A- 2. Summarize patterns of relatedness in each window  ***  #\n") }

  pcs_all <- lostruct::eigen_windows(snps, k=k) # k=2 : first two PCs
  pcs_row_na <- unique(which(is.na(pcs_all), arr.ind=T)[,"row"]) #rows windows with NA value
  pcs_row_ok <- unique(which(!is.na(pcs_all), arr.ind=T)[,"row"]) #rows of windows without NA values
  pcs <- pcs_all[pcs_row_ok,] #final data frame with only windows without NA value

  if(length(pcs_row_na)>0) {
    #pcs_na = pcs_all[pcs_row_na,]
    #stop(paste("\nWARNING: There are", length(pcs_row_na), "window(s) with missing data that were not removed. Local PCA analysis will encounter issues to calculate dissimilarity matrices. Use na.rm=TRUE or replace missing values.\n"))
    if(verbose) { cat(paste("\nThere are", length(pcs_row_na), "window(s) with missing data.\n")) }
  }


  #  ** A-3. Measure the dissimilarity in relatedness between each pair of window **

  if(verbose) { cat("\n#  ***  A- 3. Measure the dissimilarity in relatedness between each pair of window  ***  #\n") }

  pcdist <- lostruct::pc_dist(pcs, npc=2)


  #  ** 4a. Summarize dissimilarity matrices using multidimensional scaling (MDS) **

  if(verbose) { cat("\n#  ***  A- 4a. Summarize dissimilarity matrices using multidimensional scaling (MDS)  ***  #\n") }

  fit2d <- stats::cmdscale(pcdist, eig=T, k=nMDS)


  #  ** 4b. Visualize MDS **

  if(plots) {

    if(verbose) { cat("\n#  ***  A- 4b. Generate and Visualize MDS table  ***  #\n") }

    if(!dir.exists(outDir)) {
      dir.create(path=outDir, showWarnings=FALSE, recursive=TRUE)
      if(verbose) { cat(paste0("\n  ---  Created folder \"", outDir, "\" to store MDS plots  ---  \n")) }
    }
    if(!endsWith(outDir,"/")) { outDir <- paste0(outDir,"/") }

    df_win <- read.delim(pos_file) #table with windows start/end positions along the genome
    colnames(df_win) <- c("ID","Chr","posB","posE")
    df <- df_win[pcs_row_ok,] #keep only windows without NA values

    for (i in c(1:nMDS)) {

      MDS <- fit2d$points[,i]
      df <- cbind(df, as.data.frame(MDS))
      plotfile_mds <- paste0(outDir, basename(gsub(".bcf$|.vcf.gz$","",vcf)), ".", winsize, wintype, ".MDS", i, ".pdf")

      if(!file.exists(plotfile_mds)) {

        p <- ggplot(data=df, ggplot2::aes(x=.data$posB, y=.data$MDS, col=.data$Chr)) +
             geom_point(size=5, alpha=0.5) +
             facet_grid(.~Chr, scales="free_x", space="free_x") +
             xlab("Chromosome positions (Mb)") + ylab(paste0("MDS",i," coordinates")) +
             scale_x_continuous(labels = scales::unit_format(unit="", scale=1e-6), breaks = breaks_fun) +
             theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill = 'white'),
                   strip.text = element_text(size=18, face="bold"),
                   axis.line = element_line(colour="black"),
                   axis.text.y = element_text(size=20),
                   axis.text.x = element_text(size=20, angle=75, hjust=1, vjust=1),
                   axis.title.x = element_text(size=30, vjust=-0.4),
                   axis.title.y = element_text(size=30, vjust=1.5),
                   legend.position = "none")

        ggsave(filename=plotfile_mds, plot=p, height=7, width=34, units="in")
        if(verbose) { cat(paste0("\n   --->  plot saved as : ", plotfile_mds, ".\n")) }

      }

      names(df)[names(df) == "MDS"] <- paste0("MDS",i) #rename MDS columns into MDS1, MDS2, ..., MDS40

    }

  } else {

    if(verbose) { cat("\n#  ***  A- 4b. Generate full MDS table  ***  #\n") }

    df_win <- read.delim(pos_file) #table with windows start/end positions along the genome
    df <- cbind(df_win[pcs_row_ok,], as.data.frame(fit2d$points)) #keep only windows without NA values
    colnames(df) <- c("ID","Chr","posB","posE", paste0("MDS",1:nMDS))

  }


  ######################
  #  Function results  #
  ######################

  return(list(winfun=snps, pcs_all=pcs_all, rows_na=pcs_row_na, pcs=pcs, pcdist=pcdist, fit2d=fit2d, data=df))


}

