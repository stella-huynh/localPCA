test_loPCA <- function() {

  #---------- Functions ----------
  `%>%` <- magrittr:: `%>%`
  library(devtools)
  library(ggplot2)

  #---------- Setting ------------
  vcf_in <- "inst/extdata/test_data.vcf.gz"
  pos_file <- "inst/extdata/test_data.snp100.pos.txt"
  pop_file <- "inst/extdata/test_data.popFile.txt"

  #variables
  nMDS=40
  FDR_thr=0.05
  Gtest_thr=0.05
  min_nwin_H1=4
  nwin_gap=5
  nwin_cluster=3
  concat_redundant=TRUE
  clustsize=10
  corr_thr=0.6
  plots=TRUE
  prefix="loPCA_output"
  outDir="localPCA"
  verbose=TRUE
  tab_regions=regList$tabF
  popFile=pop_file
  minSize=2
  kmeans.method="euclidean"
  kmeans.k=3
  kmeans.iter=100
  keep_tmp=FALSE
  plots=TRUE
  bcftools="bcftools"
  vcftools="vcftools"
  bgzip="bgzip"


  #---------- Test 1 - checking main localPCA function ----------
  winList <- loPCA(vcf = vcf_in, pos_file = pos_file, winsize = 100)

  #---------- Test 2 - checking outlier region detection functions ----------
  regList <- getRegion(data = winList$data)

  #---------- Test 3 - checking inversion detection function ----------
  invList <- detectInv(vcf_in="inst/extdata/test_data.vcf.gz", tab_regions=regList$tabF, popFile=pop_file, kmeans.iter=10)


  return(print("Done..."))

}


# https://covid19r.github.io/documentation/unit-tests-for-data-packages.html
