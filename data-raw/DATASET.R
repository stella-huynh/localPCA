## code to prepare `DATASET` dataset goes here
vcf_in <- "inst/extdata/test_data.vcf.gz"
pos_file <- "inst/extdata/test_data.snp100.pos.txt"
pop_file <- "inst/extdata/test_data.popFile.txt"

if(!file.exists("data/winList.rda")) {
  winList <- loPCA(vcf = vcf_in, pos_file = pos_file, winsize = 100)
}
if(!file.exists("data/regList.rda")) {
  regList <- getRegion(data = winList$data)
}
if(!file.exists("data/invList.rda")) {
  invList <- detectInv(vcf_in="inst/extdata/test_data.vcf.gz", tab_regions=regList$tabF, popFile=pop_file, kmeans.iter=10)
}

usethis::use_data(winList, overwrite = TRUE)
usethis::use_data(regList, overwrite = TRUE)
usethis::use_data(invList, overwrite = TRUE)
