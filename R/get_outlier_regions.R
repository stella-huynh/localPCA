#####################################
##   1. Identify outlier windows   ##
#####################################


detect_outlierWin <- function(tab, outDir, nMDS=40, plots=TRUE) {

  ## tab : data.frame with MDS coordinates for all MDS dimensions and all genomic windows


  # B-1. calculate mean-centered MDS

  if(verbose) { cat("\n\n#  ***  B- 1. Calculate mean-centered MDS  ***  #\n\n") }

  tab_cr = tab

  for (i in paste0("MDS",c(1:nMDS))) {
    tab_cr = tab_cr %>% mutate(MDScr = (tab[,i]-mean(tab[,i])) / sd(tab[,i]))
    tab_cr[,i] = NULL
    names(tab_cr)[names(tab_cr) == "MDScr"] = i
  }


  # B-2. assess p-value of each mean-centered MDS

  if(verbose) { cat("\n\n#  ***  B- 2. Assess p-value for each mean-centered MDS  ***  #\n\n") }

  tab_pval = tab_cr

  for (i in paste0("MDS",c(1:nMDS))) {
    tab_pval = tab_pval %>% mutate(pval = pnorm(abs(tab_cr[,i]), mean=0, sd=1, lower.tail=F))
    tab_pval[,i] = NULL
    names(tab_pval)[names(tab_pval) == "pval"] = i
  }


  # B-3. apply FDR (False Discovery Rate) to identify outliers

  if(verbose) { cat("\n\n#  ***  B- 3. Apply FDR (False Discovery Rate) to identify outliers  ***  #\n\n") }

  tab_FDR = tab[,c("ID","Chr","posB","posE")]

  for (i in paste0("MDS",c(1:nMDS))) {
    FDR_H1 = p.fdr(pval=tab_pval[,i], threshold=FDR_thr)$`Reject Vector`
    tab_FDR = cbind(tab_FDR, FDR_H1)
    names(tab_FDR)[names(tab_FDR) == "FDR_H1"] = i
  }


  # plot MDS with outliers

  if(plots) {
    if(verbose) { cat("\n\n   ====   Plot MDS with outliers   ====   \n\n") }

    for (i in c(1:nMDS)) {

      tmp = cbind(tab[,c("ID","Chr","posB","posE")], tab_cr[,paste0("MDScr",i)], tab_FDR[,paste0("FDR",i,"_H1")])

      plotfile = paste0(outDir,basename(gsub(".bcf$|.vcf.gz$","",vcffile)), ".", winsize, wintype,".MDScr",i,".pdf")

      if(!file.exists(plotfile)) {
        p_out = ggplot(tmp, aes(x=posB, y=-log(abs(tmp[,paste0("MDScr",i)]), base=10), col=tmp[,paste0("FDR",i,"_H1")])) +
          geom_point(size=5, alpha=0.5) +
          facet_grid(.~Chr, scales="free_x", space="free_x") +
          xlab("Chromosome positions (Mb)") + ylab(paste0("-log(pvalue) MDS",i)) +
          scale_x_continuous(labels = unit_format(unit="", scale=1e-6), breaks = breaks_fun) +
          theme(panel.grid.major = element_blank(),
                panel.grid.minor = element_blank(),
                panel.background = element_rect(fill = 'white'),
                strip.text = element_text(size=18, face="bold"),
                axis.line = element_line(colour="black"),
                axis.text.y = element_text(size=20),
                axis.text.x = element_text(size=20,angle=75,hjust=1,vjust=1),
                axis.title.x = element_text(size=30,vjust=-0.4),
                axis.title.y = element_text(size=30,vjust=1.5),
                legend.position = "none")
        ggsave(filename=plotfile, plot=p_out, height=7, width=34, units="in")
        if(verbose) { cat(paste0("\n   --->  plot saved as : ", plotfile, ".\n")) }
      }

    }
  }


  return(list(coor.raw = tab,
              coor.cr = tab_cr,
              pval = tab_pval,
              fdr = tab_FDR
              )
         )

}



#######################################
##   2. Clustering outlier windows   ##
#######################################


clusterWin <- function(outWin_IDs, nwin_gap, nwin_cluster) {

  c = 0
  nb_outWin = length(outWin_IDs)
  list_result  =  list()
  new_l = list()
  previous_ID = 1

  #for each significant windows positions
  for (i in outWin_IDs) {
    ID = i
    diff = ID - previous_ID
    if (diff < nwin_gap) {
      if (!(previous_ID %in% new_l) & c>0) {
        new_l[[ (length(new_l) + 1) ]]  =  previous_ID
      }
      new_l[[ (length(new_l) + 1) ]]  =  ID
    } else {
      if (length(new_l) >= nwin_cluster) { # if enough significant windows in the cluster (ex: > 3)
        list_result[[ (length(list_result) + 1) ]]  =  new_l # add a cluster (new_l) to list result
      }
      new_l = list() #if diff > gap a new cluster is created/initialized
    }
    previous_ID = ID
    c = c+1
  }

  #necessary if diff < gap at the last loop
  if ((length(new_l) >= nwin_cluster) & (diff < nwin_gap) & (c == nb_outWin)) {
    list_result[[(length(list_result) + 1)]]  =  new_l # add the last cluster (new_l) to list result
  }

  return(list_result) #return a list of lists, each sub list being a cluster

}



###################################
##   3. Define outlier regions   ##
###################################


#' Define specific genomic regions harbouring population structure deviating from the genome-wide structure.
#'
#' @param winList Output list from \bold{detect_outlierWin}. Must contains 4-levels named: "coor.raw", "coor.cr", "pval" and "fdr"
#' @param nMDS Number of MDS dimension inferred, (integer, default=40)
#' @param Gtest_thr Minimum threshold for G-test significance (float, default=0.05)
#' @param min_nwin Minimum number of genomic windows required to keep a chromosome for further analysis (integer, default=4)
#' @param nwin_gap Mininum number of non-significant genomic windows between two significant genomic windows (integer, default=3)
#' @param nwin_cluster Minimum number of genomic windows within a regions to define it as an "outlier region" potentially carrying strong within-population structure signal
#' @param concat_redundant Concatenate redundant overlapped genomic regions between MDS dimensions (boolean, default=TRUE)
#' @param clustsize Minimum size (in number of genomic windows) of the intersect portion between two overlapping genomic regions. Required if \option{concat_redundant=TRUE}. (integer, default=10)
#' @param corr_thr Minimum threshold for correlation test. Required if \option{concat_redundant=TRUE}. (default=0.6)
#'
#' @return If \option{concat_redundant=TRUE}, returns a list of 3 elements :
#'         \tabular{ll}{
#'           tabF \tab a data.frame of the concatenated outlier regions\cr
#'           tabF_raw \tab a data.frame of the "raw" outlier regions detected\cr
#'           listRegionWin \tab a list of regions with their respective window IDs encompassed
#'         }
#'         If \option{concat_redundant=TRUE}, returns a list of 2 elements (note here that tabF refers to the "raw" regions):
#'         \tabular{ll}{
#'           tabF \tab a data.frame of the "raw" outlier regions detected\cr
#'           listRegionWin \tab a list of regions with their respective window IDs encompassed
#'         }
#'
#' @export
#'


getRegion <- function(winList, nMDS=40, Gtest_thr=0.05, min_nwin_H1=4, nwin_gap=5, nwin_cluster=3, concat_redundant=TRUE, clustsize=10, corr_thr=0.6) {

  df_win = winList$mds[,c("ID","Chr","posB","posE")]
  tab_region = data.frame(ID=integer(), MDS=character(), chr=character(), start=integer(), end=integer(), noutliers=integer(), stringsAsFactors=F)

  count_region = 1
  listRegionH1 = list()

  for (i in 1:nMDS) {

    MDS = paste0("MDS",i)
    tab = as.data.frame(cbind(df_win,
                              coor.raw = winList$coor.raw[,MDS],
                              coor.cr = winList$coor.cr[,MDS],
                              pval = winList$pval[,MDS],
                              fdr = winList$fdr[,MDS]
                              )
                        )

    for (chr in unique(df_win$Chr)) {

      nwin = nrow(tab %>% filter(Chr==chr))
      nwin_H1 = nrow(tab %>% filter(Chr==chr, fdr=="Reject.H0")) #number of significant windows

      ## 1. If enough significant windows are on the chr :
      if (nwin_H1 > min_nwin_H1) {
        nwin_H1_tot = nrow(tab %>% filter(fdr=="Reject.H0")) #number of significant windows for the MDS
        mat_gtest = rbind(tot = c(nwin_H1 = nwin_H1_tot, nwin_H0 = nwin_tot-nwin_H1_tot),
                          chr = c(nwin_H1 = nwin_H1,     nwin_H0 = nwin-nwin_H1        ))
        pval_Gtest = GTest(mat_gtest)$p.value

        ## 2. If windows are significantly gathered on the chromosome :
        if (pval_Gtest < Gtest_thr) {

          ## 3. Clustering of neighbouring significant windows
          tab_p = tab %>% dplyr::select(ID, Chr, posB, posE, coor.raw, fdr) %>%
            filter(Chr==chr) %>%
            filter_at("coor.raw", all_vars(. > 0)) %>%
            filter_at("fdr", all_vars(. == "Reject.H0"))

          if (nrow(tab_p) > 0) {
            clust_pos = clustering(tab_p[,"ID"], nwin_gap, nwin_cluster)

            if (length(clust_pos) > 0) {
              for (clust_p in clust_pos) {
                # nb of outliers windows in the cluster
                nwin_H1_p = length(clust_p)
                # beginning of the cluster
                first_w = tab_p %>% filter(ID==clust_p[1]) %>% dplyr::select(posB)
                first_pos = first_w[1,1]
                # end of the cluster
                last_w = tab_p %>% filter(ID==clust_p[nwin_H1_p]) %>% dplyr::select(posE)
                last_pos = last_w[1,1]
                # unique ID for each region
                ID_region = count_region
                tab_region[nrow(tab_region)+1,] = c(as.list(ID_region),MDS,chr,as.list(c(first_pos,last_pos,nwin_H1_p)))
                count_region = count_region + 1
                listRegionH1[[(length(listRegionH1)+1)]] = clust_p
              }
            }
          }

          ### for negative MDS :
          tab_n = tab %>% dplyr::select(ID, Chr, posB, posE, coor.raw, fdr) %>%
            filter(Chr==chr) %>%
            filter_at("coor.raw", all_vars(. < 0)) %>%
            filter_at("fdr", all_vars(. == "Reject.H0"))

          if (nrow(tab_n) > 0) {
            clust_neg = clusterWin(tab_n[,"ID"], nwin_gap, nwin_cluster)

            if (length(clust_neg) > 0) {
              for (clust_n in clust_neg) {
                # nb of outliers windows in the cluster
                nwin_H1_n = length(clust_n)
                #beginning of the cluster
                first_w = tab_n %>% filter(ID==clust_n[1]) %>% dplyr::select(posB)
                first_pos = first_w[1,1]
                # end of the cluster
                last_w = tab_n %>% filter(ID==clust_n[nwin_H1_n]) %>% dplyr::select(posE)
                last_pos = last_w[1,1]
                # unique ID for each region
                ID_region = count_region
                tab_region[nrow(tab_region) + 1,] = c(as.list(ID_region),MDS,chr,as.list(c(first_pos,last_pos,nwin_H1_n)))
                count_region = count_region + 1
                listRegionH1[[(length(listRegionH1) + 1)]] = clust_n
              }
            }
          }

        }
      }
    }
  }

  tab_region = tab_region %>%
    arrange(by_group=chr,start) %>%
    mutate(SizeRegions=as.numeric(sprintf("%.1f",(end-start)/1000000)))


  # if option "concat_redundant" is TRUE (default), returns a list of 3 elements
  if(concat_redundant) {

    tabF = concatReg(winList=winList, rtab=tab_region, listRegionWin=listRegionH1, clustsize=clustsize, corr_thr=corr_thr)
    return(list(tabF=tabF, tabF_raw=tab_region, listRegionWin=listRegionH1))

  # if option "concat_redundant" is FALSE, returns a list of 2 elements
  } else {

    return(list(tabF=tab_region, listRegionWin=listRegionH1))

  }


}



##################################################
##   4. Concatenate redundant outlier regions   ##
##################################################


concatReg <- function(winList, rtab, listRegionWin, clustsize=10, corr_thr=0.6) {

  options(warn = 2)

  ###  winList : Output list from detect_outlierWin(). Must contains 4-levels named: "coor.raw", "coor.cr", "pval" and "fdr" [winList]
  ###  rtab : table of outlier genomic regions assessed generated by the R-function "get_outlier_regions()" [tab_region]
  ###  listRegWin : list of outlier genomic regions identified by the R-function "get_outlier_regions()", each containing nested list(s) with the indexes of the outlier windows that are comprised in each of them [listRehionH1]
  ###  clustsize : minimum size of intersect portion (in number of windows) between two overlapping regions [nwin_cluster]
  ###  corr_thr : minimum threshold for correlation test [0.6]

  ###  Each genomic region is processed one after the other, by comparing it to the previous region(s). Each time the previous region overlap the current region, it is stored in the R-object "tmp" and will remain until it is concatenated or non-overlapping with the current region. The current region may therefore be compared with more than one previous region if these latter are still overlapping the current region but with no correlated MDS values (in such case they are not concatenated and kept separated). Once the current region does not overlap with a previous region, this latter will be removed from "tmp".
  ###  The previous region(s) is stored in "rtab_i0/rtab_i" and "tmp".
  ###  The current region is stored "rtab_j".

  rtab = rtab %>% mutate(chr = factor(chr, ordered=T, levels = unique(mixedsort(chr)))) %>% arrange(chr,start,end,desc(SizeRegions))
  df_win = winList$coor.raw[,c("ID","Chr","posB","posE")]

  tabF = data.frame(matrix(ncol=ncol(rtab)+1,nrow=0))
  colnames(tabF) = c("nID",colnames(rtab))
  to_remove = list() #list of genomic regions' IDs judged as redundant and therefore to be ignored in final table.

  count=0 #used as IDs for the final concatenated region
  tmp=NULL #subset of rtab containing any region preceding (rtab_i) the current region being assessed (rtab_j), and which have not been merged and are still overlapping it


  for(n in 2:nrow(rtab)) {

    merged_IDs = NULL
    intersect_IDs = NULL

    if(is.null(tmp)) {
      rtab_i0 = rtab[n-1,] #compare to single previous region
    } else {
      rtab_i0 = tmp #compare to all the accumulated overlapping previous regionS
    }
    rtab_j = rtab[n,]

    rtab_i_IDs = rtab_i0[,"ID"]

    for(m in rtab_i_IDs) {

      rtab_i = rtab_i0[which(rtab_i0$ID==m),]
      chr_i = rtab_i$chr
      chr_j = rtab_j$chr

      if(!(rtab_i[,"ID"] %in% to_remove)) {

        #if windows pair are on the same chromosome
        if(chr_i == chr_j) {

          MDS_i   = rtab_i$MDS
          MDS_j   = rtab_j$MDS
          start_i = rtab_i$start
          start_j = rtab_j$start
          end_i   = rtab_i$end
          end_j   = rtab_j$end

          # Nested windows
          # case1 : MDSa1 > MDSa2
          if ((start_i <= start_j) && (end_i >= end_j)) { do_clust=TRUE ; case=1
          # case2 : MDSa2 > MDSa1
          } else if ((start_j <= start_i) && (end_j >= end_i)) { do_clust=TRUE ; case=2
          # Overlapping windows
          # case3 : MDSa1_end overlap MDSa2_start
          } else if ((start_i <= start_j) && (start_j <= end_i) && (end_j >= end_i)) { do_clust=TRUE ; case=3
          # case4 : MDSa1_start overlap MDSa2_end
          } else if ((start_j <= start_i) && (start_i <= end_j) && (end_i >= end_j)) { do_clust=TRUE ; case=4
          #no overlap
          } else { do_clust=FALSE ; case=NA
          }

          #Correlation test (if overlapping windows pair)
          if(do_clust==TRUE) {

            ID_i = rtab_i$ID
            ID_j = rtab_j$ID
            if(is.character(ID_i)){
              outliers_i = unique(sort(unlist(listRegWin[as.numeric(unlist(strsplit(ID_i,",|/")))]))) #extract indexes of the windows comprised in the genomic region rtab_i (for concatenated previous regions)
            } else {
              outliers_i = unlist(listRegWin[ID_i]) #extract indexes of the windows comprised in the genomic region rtab_i (for single previous region)
            }
            outliers_j = unlist(listRegWin[ID_j]) #extract indexes of the windows comprised in the genomic region rtab_j
            firstpos_i = outliers_i[1] #extract position of the first window (ie. beginning of the region)
            firstpos_j = outliers_j[1] #extract position of the first window
            lastpos_i  = rev(outliers_i)[1]
            lastpos_j  = rev(outliers_j)[1]

            tab = as.data.frame(cbind(df_win, winList$coor.raw[,c(MDS_i,MDS_j)]))

            corrmat_i = tab %>% filter((ID>=firstpos_i)&(ID<=lastpos_i)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))  #extract MDS values of those windows for region rtab_i
            corrmat_j = tab %>% filter((ID>=firstpos_j)&(ID<=lastpos_j)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))  #extract MDS values of those windows for region rtab_i
            mcor_i = cor(corrmat_i, method="pearson") #correlation test for MDS values within rtab_i region coordinates
            mcor_j = cor(corrmat_j, method="pearson") #correlation test for MDS values within rtab_j region coordinates
            if(MDS_i == MDS_j) {
              rcor_i = mcor_i[1] ; rcor_j = mcor_j[1]
            } else {
              rcor_i = mcor_i[2] ; rcor_j = mcor_j[2]
            }

            rcorr = mean(c(abs(rcor_i),abs(rcor_j))) #Pearson correlation. A value of 0 implies that there is no linear dependency between the two regions. An absolute value of 1 implies a linear equation describing the relationship between the two regions.

            #if significant correlation ...
            if ( abs(rcor_j) > corr_thr || abs(rcor_i) > corr_thr ) {

              if(case==1) { # --> keep ID_i and remove ID_j
                to_remove = append(to_remove, ID_j)
                if(is.null(tmp)) {
                  tmp = rtab_i
                }

              } else if(case==2) { # ---> keep ID_j and remove ID_i
                to_remove = append(to_remove, ID_i)
                if(is.null(tmp)) {
                  tmp = rtab_j
                } else {
                  tmp = rbind(tmp,rtab_j)
                }

              } else if(case==3) { # --> merge ID_i + ID_j regions. Keep MDS values of the longest region of them two. For that a new region row needs to be created and added to the final result table..
                if(!is.null(tmp)) {
                  if(nrow(tmp)>1) { merged_IDs = c(merged_IDs, c(ID_i,ID_j)) }
                }
                if((end_i-start_i) > (end_j-start_j)) {
                  MDS_ij = MDS_i
                } else {
                  MDS_ij = MDS_j
                }
                new_region = data.frame(ID=paste0(ID_i,",",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_j, noutliers=length(union(outliers_i,outliers_j)), SizeRegions=round((end_j-start_i)/1000000,1))
                if(is.null(tmp) || nrow(tmp)==1) {
                  tmp = new_region
                } else {
                  tmp[which(tmp$ID==m),] = new_region
                }

              } else if(case==4) { # --> same as for case 3.
                if(!is.null(tmp)) {
                  if(nrow(tmp)>1) { merged_IDs = c(merged_IDs, c(ID_i,ID_j)) }
                }
                if((end_i-start_i) > (end_j-start_j)) {
                  MDS_ij = MDS_i
                } else {
                  MDS_ij = MDS_j
                }
                new_region = data.frame(ID=paste0(ID_i,",",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_i, noutliers=length(union(outliers_i,outliers_j)), SizeRegions=round((end_i-start_j)/1000000,1))
                if(is.null(tmp) || nrow(tmp)==1) {
                  tmp = new_region
                } else {
                  tmp[m,] = new_region
                }
              }

              #... else if size overlap > "clustsize" windows ; do corr-test for large intersect only
            } else {
              intersect_ij = intersect( seq(firstpos_i,lastpos_i), seq(firstpos_j,lastpos_j) )

              if(length(intersect_ij) >= clustsize) {

                if(case==1) {
                  corrmat_ij = tab %>% filter((ID>=firstpos_j)&(ID<=lastpos_j)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
                } else if(case==2) {
                  corrmat_ij = tab %>% filter((ID>=firstpos_i)&(ID<=lastpos_i)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
                } else if(case==3) {
                  corrmat_ij = tab %>% filter((ID>=firstpos_j)&(ID<=lastpos_i)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
                } else if(case==4) {
                  corrmat_ij = tab %>% filter((ID>=firstpos_i)&(ID<=lastpos_j)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
                }

                mcor_ij = cor(corrmat_ij, method="pearson")
                if(MDS_i == MDS_j) {
                  rcor_ij = mcor_ij[1]
                } else {
                  rcor_ij = mcor_ij[2]
                }

                if( abs(rcor_ij) > corr_thr ) {
                  if(!is.null(tmp)) {
                    if(nrow(tmp)>1) { intersect_IDs = c(intersect_IDs, c(ID_i,ID_j)) }
                  }

                  if(case==1) { # --> keep both ID_i and ID_j regions separate AND the correlated intersect between them. Keep MDS value of the longest region (here ID_i).
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_j, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=rtab_j[,"SizeRegions"])

                  } else if(case==2) { # --> keep both ID_i and ID_j regions separate AND the correlated intersect between them. Keep MDS value of the longest region (here ID_j).
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_i, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=rtab_i[,"SizeRegions"])

                  } else if(case==3) { # --> keep both ID_i and ID_j regions separate AND the correlated intersect between them. Keep MDS value of the longest region (either ID_i or ID_j).
                    if((end_i-start_i) > (end_j-start_j)) {
                      MDS_ij = MDS_i
                    } else {
                      MDS_ij = MDS_j
                    }
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_i, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=round((end_i-start_j)/1000000,1))

                  } else if(case==4) { # --> keep both ID_i and ID_j regions separate AND the correlated intersect between them. Keep MDS value of the longest region (either ID_i or ID_j).
                    if((end_i-start_i) > (end_j-start_j)) {
                      MDS_ij = MDS_i
                    } else {
                      MDS_ij = MDS_j
                    }
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_j, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=round((end_j-start_i)/1000000,1))
                  }

                  tmp = rbind(tmp, new_region, rtab_j)

                } else {

                  if(is.null(tmp)) {
                    tmp = rtab[(n-1):n,]
                  } else if(!(rtab[n,"ID"] %in% to_remove)) {
                    tmp = rbind(tmp, rtab[n,])
                  }

                }

              } else {

                if(is.null(tmp)) {
                  tmp = rtab[(n-1):n,]
                } else if(!(rtab[n,"ID"] %in% to_remove)) {
                  tmp = rbind(tmp, rtab[n,])
                }

              }
            }

          } else if(do_clust==FALSE) {

            if(is.character(tmp[,"ID"])) {
              tmp_IDs = as.numeric(unlist(strsplit(tmp[,"ID"],",|/")))
            } else {
              tmp_IDs = tmp[,"ID"]
            }

            if(!rtab_j["ID"] %in% tmp_IDs) {

              if(!is.null(tmp)) {

                if(is.character(m)) {
                  m_IDs = as.numeric(unlist(strsplit(m,",|/")))
                } else {
                  m_IDs = m
                }
                m_tmp = tmp[grep(paste0("\\b",m_IDs,"\\b",collapse="|"),tmp[,"ID"]),]

                # if non-overlapping regions are included in other merged regions that are overlapping the current region (rtab_j);
                # remove non-overlapping regions from tmp

                # ... remove non-overlapping regions from tmp
                if(nrow(m_tmp[which(m_tmp$end >= rtab_j$start),]) > 0) {
                  to_remove = append(to_remove, m_tmp[which(m_tmp$end < rtab_j$start),"ID"])

                  # ... otherwise keep longest one if several regions are overlapping each other but not the current region (rtab_j)
                } else {
                  m_tmp = m_tmp[which(m_tmp$end < rtab_j$start),]
                  m_tmp_IDs = m_tmp[,"ID"]
                  if(nrow(m_tmp)>1) { #if more than one region, keep longest one
                    m_tmp = m_tmp[which(m_tmp$SizeRegions==max(m_tmp$SizeRegions)),]
                    if(nrow(m_tmp)>1) { #if more than one region have same max_length, keep the one that has the highest number of outlier windows
                      m_tmp = m_tmp[which(m_tmp$noutliers==max(m_tmp$noutliers)),]
                      if(nrow(m_tmp)>1) { #if more than one region have same max_length and same number of outlier windows, just pick the first one
                        m_tmp = m_tmp[1,]
                      }
                    }
                  }

                  to_remove = append(to_remove, m_tmp_IDs[which(m_tmp_IDs!=m_tmp[,"ID"])])
                  rm(m_tmp_IDs)
                }

                if(nrow(m_tmp)==1) {

                  count = count + 1
                  tabF = rbind(tabF, cbind(nID=count, m_tmp))

                  if(m==tail(rtab_i_IDs,1)) {
                    tmp = NULL
                  } else {
                    tmp = tmp[which(tmp$ID != m),]
                  }

                } else if (nrow(m_tmp)==0) {
                  print(paste0("[",n,",",m,"] !! ISSUE : tmp is unexpectedly empty !!"))
                }
                rm(m_tmp)

              } else {

                count = count + 1
                tabF = rbind(tabF, cbind(nID=count, rtab[n-1,]))

              }

            }

          }

          #if different chromosome
        } else {

          count = count + 1

          if(!is.null(tmp)) {
            tabF = rbind(tabF, cbind(nID=count, rtab_i))
            if(m==tail(rtab_i_IDs,1)) { tmp = NULL } #once all regions in tmp have been processed, reset "tmp" to NULL.

          } else {
            tabF = rbind(tabF, cbind(nID=count, rtab[n-1,]))
          }

        }

      }

    }


    if(!is.null(tmp)) {

      if(m==tail(rtab_i_IDs,1)) {
        tmp = tmp[which(!tmp$ID %in% to_remove),]

        if(nrow(tmp)>0) {

          if(any(duplicated(tmp[,"ID"]))) {
            dup_IDs = as.numeric(tmp[,"ID"][duplicated(tmp[,"ID"])])
            tmp = tmp[!duplicated(tmp[,"ID"]),]
          }

          # remove region that got added to tmp if it has been later merged to another region in tmp
          if(is.character(tmp[,"ID"])) {

            if(any(duplicated(unlist(strsplit(tmp[,"ID"],","))))) {

              dup_IDs = unique(unlist(strsplit(tmp[,"ID"],","))[duplicated(unlist(strsplit(tmp[,"ID"],",")))])
              tmp = tmp[which(!as.character(tmp[,"ID"]) %in% dup_IDs),]

              if(nrow(tmp[grep(paste0("\\b",dup_IDs,"\\b",collapse="|"), tmp[,"ID"]),])>0) {
                cat(paste0("COMMENT : Removed region(s) [",dup_IDs,"] from tmp because it was/were found merged with other region(s) in tmp. \nNew tmp:\n")) ; print(tmp)
              } else {
                cat(paste0("ERROR : check row ",n,", IDs [",dup_IDs,"]."))
              }

            }

          }

          # re-merge merged regions containing common IDs
          if(any(duplicated(merged_IDs))) {
            dup_rows = tmp[grepl(paste0("\\b",merged_IDs[duplicated(merged_IDs)],"\\b",collapse="|"), tmp[,"ID"]),]

            if(nrow(dup_rows[which(dup_rows$start==min(dup_rows$start) & dup_rows$end==max(dup_rows$end)),])>0) {
              new_region = dup_rows[which(dup_rows$start==min(dup_rows$start) & dup_rows$end==max(dup_rows$end)),]

            } else {
              dup_IDs = paste0(unique(unlist(strsplit(dup_rows[,"ID"],","))),collapse=",")
              MDS_ij = unique(dup_rows[which(dup_rows$SizeRegions==min(dup_rows$SizeRegions)),"MDS"])
              noutliers = length(unlist(listRegWin[as.numeric(unique(unlist(strsplit(dup_rows[,"ID"],",|/"))))]))
              new_region = data.frame(ID=dup_IDs, MDS=MDS_ij, chr=chr_j, start=min(dup_rows$start), end=max(dup_rows$end), noutliers=noutliers, SizeRegions=round((max(dup_rows$end)-min(dup_rows$start))/1000000,1))
            }
            tmp = rbind(new_region,tmp[!grepl(paste0("\\b",merged_IDs[duplicated(merged_IDs)],"\\b",collapse="|"), tmp[,"ID"]),])

          }

        } else {
          tmp = NULL
        }

      }

    } else {
      to_remove = list()
    }

  }

  # when reaching the end of the "raw" region table, concatenate the last remaining region.
  if(n==nrow(rtab)) {

    count = count + 1

    if(is.character(tmp[,"ID"])) { tmp_IDs = unique(as.numeric(unlist(strsplit(tmp[,"ID"],",|/")))) } else { tmp_IDs = tmp[,"ID"] }
    if(is.character(tabF[,"ID"])) { tabF_IDs = unique(as.numeric(unlist(strsplit(tabF[,"ID"],",|/")))) } else { tabF_IDs = tabF[,"ID"] }

    if(any(!rtab_j$ID %in% tmp_IDs) && any(!rtab_j$ID %in% tabF_IDs)) {
      tmp = rbind(tmp,rtab[n,])
    }

    tabF = rbind(tabF, cbind(nID=count, tmp))
  }


  return(tabF)

}


