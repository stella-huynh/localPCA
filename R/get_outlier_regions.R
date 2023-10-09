detect_outlierWin <- function(tab, outDir, nMDS=40, plots=TRUE) {

  ## tab : data.frame with MDS coordinates for all MDS dimensions and all genomic windows


  # B-1. calculate mean-centered MDS

  if(verbose) { cat("\n\n#  ***  B- 1. Calculate mean-centered MDS  ***  #\n\n") }

  for (i in c(1:nMDS)) {
    tab = tab %>% mutate(MDScr = (tab[,paste0("MDS",i)]-mean(tab[,paste0("MDS",i)])) / sd(tab[,paste0("MDS",i)]))
    names(tab)[names(tab) == "MDScr"] = paste0("MDScr",i)
  }


  # B-2. assess p-value of each mean-centered MDS

  if(verbose) { cat("\n\n#  ***  B- 2. Assess p-value for each mean-centered MDS  ***  #\n\n") }

  for (i in c(1:nMDS)) {
    tab = tab %>% mutate(pval = pnorm(abs(tab[,paste0("MDScr",i)]), mean=0, sd=1, lower.tail=F))
    names(tab)[names(tab) == "pval"] = paste0("pval",i)
  }


  # B-3. apply FDR (False Discovery Rate) to identify outliers

  if(verbose) { cat("\n\n#  ***  B- 3. Apply FDR (False Discovery Rate) to identify outliers  ***  #\n\n") }

  for (i in c(1:nMDS)) {
    FDR_H1 = p.fdr(pval=tab[,paste0("pval",i)], threshold=FDR_thr)$`Reject Vector`
    tab = cbind(tab, FDR_H1)
    names(tab)[names(tab) == "FDR_H1"] = paste0("FDR",i,"_H1")
  }


  # plot MDS with outliers

  if(verbose) { cat("\n\n   ====   Plot MDS with outliers   ====   \n\n") }

  for (i in c(1:nMDS)) {

    plotfile = paste0(outDir,basename(gsub(".bcf$|.vcf.gz$","",vcffile)), ".", winsize, wintype,".MDScr",i,".pdf")

    if(!file.exists(plotfile)) {
      p_out = ggplot(tab, aes(x=posB, y=-log(abs(tab[,paste0("MDScr",i)]), base=10), col=FDR_H1)) +
        geom_point(size=5, alpha=0.5) +
        facet_grid(.~Chr, scales="free_x", space="free_x") +
        xlab("Chromosome positions") + ylab(paste0("-log(pvalue) MDS",i)) +
        scale_x_continuous(labels = unit_format(unit="Mb", scale=1e-6), breaks = breaks_fun) +
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

  return(tab)

}



clusterWin <- function(list_parameter, gap, lenght_list) {
  c = 0
  size_list_parameter = length(list_parameter)
  list_result  =  list()
  new_l = list()
  previous_ID = 1

  #for each significant windows positions
  for (i in list_parameter) {
    ID = i
    diff = ID-previous_ID
    if ((diff<gap)) {
      if (!(previous_ID %in% new_l) & c>0) {
        new_l[[(length(new_l) + 1)]]  =  previous_ID
      }
      new_l[[(length(new_l) + 1)]]  =  ID
    } else {
      if (length(new_l)>=lenght_list) { #if enough significant windows in the cluster (ex: > 3)
        list_result[[(length(list_result) + 1)]]  =  new_l # add a cluster (new_l) to list result
      }
      new_l = list() #if diff > gap a new cluster is created/initialized
    }
    previous_ID = ID
    c = c+1
  }

  #necessary if diff < gap at the last loop
  if ((length(new_l)>=lenght_list)&(diff<gap)&(c==size_list_parameter)) {
    list_result[[(length(list_result) + 1)]]  =  new_l # add the last cluster (new_l) to list result
  }

  return(list_result) #return a list of lists, each sub list being a cluster

}



getReg <- function(tab, nb_coord=40, Gtest_thr=0.05, min_nwin_H1=4, nwin_gap=5, nwin_cluster=3) {

  tab_region = data.frame(ID=integer(), MDS=character(), chr=character(), start=integer(), end=integer(), noutliers=integer(), stringsAsFactors=F)

  count_region = 1
  listRegionH1 = list()

  for (i in 1:nb_coord) {

    MDS = paste0("MDS",i)
    FDR_MDS = paste0("FDR",i,"_H1")

    for (chr in unique(tab$Chr)) {

      nwin = nrow(tab %>% filter(Chr==chr))
      nwin_H1 = nrow(tab %>% filter(Chr==chr, tab[,FDR_MDS]=="Reject.H0")) #number of significant windows

      ## 1. If enough significant windows are on the chr :
      if (nwin_H1 > min_nwin_H1) {
        nwin_H1_tot = nrow(tab %>% filter(tab[,FDR_MDS]=="Reject.H0")) #number of significant windows for the MDS
        mat_gtest = rbind(tot = c(nwin_H1 = nwin_H1_tot, nwin_H0 = nwin_tot-nwin_H1_tot),
                          chr = c(nwin_H1 = nwin_H1,     nwin_H0 = nwin-nwin_H1        ))
        pval_Gtest = GTest(mat_gtest)$p.value

        ## 2. If windows are significantly gathered on the chromosome :
        if (pval_Gtest < Gtest_thr) {

          ## 3. Clustering of neighbouring significant windows
          tab_p = tab %>% dplyr::select(ID,all_of(FDR_MDS),all_of(MDS),Chr,posB,posE) %>%
            filter(Chr==chr) %>%
            filter_at(vars(all_of(MDS)), all_vars(. < 0)) %>%
            filter_at(vars(all_of(FDR_MDS)), all_vars(. == "Reject.H0"))

          if (length(tab_p[,1]) > 0) {
            clust_pos = clustering(tab_p[,1], nwin_gap, nwin_cluster)

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
          tab_n = tab %>% dplyr::select(ID,all_of(FDR_MDS),all_of(MDS),Chr,posB,posE) %>%
            filter(Chr==chr) %>%
            filter_at(vars(all_of(MDS)), all_vars(. > 0)) %>%
            filter_at(vars(all_of(FDR_MDS)), all_vars(. == "Reject.H0"))

          if (length(tab_n[,1]) > 0) {
            clust_neg = clusterWin(tab_n[,1], nwin_gap, nwin_cluster)

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

  return(list(tab_region=tab_region, listRegionH1=listRegionH1))
}
