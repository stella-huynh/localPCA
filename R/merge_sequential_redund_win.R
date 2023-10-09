options(warn = 2)

merge_sequential_redund_win <- function(tab=tab, rtab=tab_region, listRegionH1=listRegionH1, clustsize=nwin_cluster, corr_thr=0.6, verbose=FALSE) {

  rtab = rtab %>% mutate(chr = factor(chr, ordered=T, levels = unique(mixedsort(chr)))) %>% arrange(chr,start,end,desc(SizeRegions))

  rtab$decision = NA
  tabF = data.frame(matrix(ncol=9,nrow=0))
  colnames(tabF) = c("nID",colnames(rtab))
  list_merge = list()
  listF_merge = list()
  to_remove = list()

  count=0
  tmp=NULL


  for(n in 2:nrow(rtab)) {

    merged_IDs = NULL ; intersect_IDs = NULL
    rtab_j = rtab[n,]

    if(is.null(tmp)) {
      rtab_i0 = rtab[n-1,]
    } else {
      rtab_i0 = tmp
    }

    rtab_i_IDs = rtab_i0[,"ID"]
    m=rtab_i_IDs[1] #to remove once debugging successful

    for(m in rtab_i_IDs) {
      rtab_i = rtab_i0[which(rtab_i0$ID==m),]
      if(isTRUE(verbose)) { cat("tmp\n") ; print(tmp) ; cat("\nrtab_i\n") ; print(rtab_i) ; cat("\nrtab_j\n") ; print(rtab_j) }

      chr_i = rtab_i$chr
      chr_j = rtab_j$chr

      if(!(rtab_i[,"ID"] %in% to_remove)) {

        #if windows pair are on the same chromosome
        if(chr_i == chr_j) {

          MDS_i   = rtab_i$MDS    ;  MDS_j   = rtab_j$MDS
          start_i = rtab_i$start  ;  start_j = rtab_j$start
          end_i   = rtab_i$end    ;  end_j   = rtab_j$end

          #nested windows
          # MDSa1 > MDSa2
          if ((start_i <= start_j) && (end_i >= end_j)) { do_clust=TRUE ; case=1
          # MDSa2 > MDSa1
          } else if ((start_j <= start_i) && (end_j >= end_i)) { do_clust=TRUE ; case=2
          #overlapping windows
          # MDSa1_end overlap MDSa2_start
          } else if ((start_i <= start_j) && (start_j <= end_i) && (end_j >= end_i)) { do_clust=TRUE ; case=3
          # MDSa1_start overlap MDSa2_end
          } else if ((start_j <= start_i) && (start_i <= end_j) && (end_i >= end_j)) { do_clust=TRUE ; case=4
          #no overlap
          } else { do_clust=FALSE ; case=NA
          }

          #Correlation test (if overlapping windows pair)
          if(do_clust==TRUE) {

            ID_i = rtab_i$ID                         ;  ID_j = rtab_j$ID
            if(is.character(ID_i)){
              outliers_i = unique(sort(unlist(listRegionH1[as.numeric(unlist(strsplit(ID_i,",|/")))])))
            } else {
              outliers_i = unlist(listRegionH1[ID_i])
            }
            outliers_j = unlist(listRegionH1[ID_j])
            firstpos_i = outliers_i[1]               ;  firstpos_j = outliers_j[1]
            lastpos_i  = rev(outliers_i)[1]          ;  lastpos_j  = rev(outliers_j)[1]

            corrmat_i = tab %>% filter((ID>=firstpos_i)&(ID<=lastpos_i)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
            corrmat_j = tab %>% filter((ID>=firstpos_j)&(ID<=lastpos_j)) %>% dplyr::select(all_of(MDS_i),all_of(MDS_j))
            mcor_i = cor(corrmat_i, method="pearson")
            mcor_j = cor(corrmat_j, method="pearson")
            if(MDS_i == MDS_j) {
              rcor_i = mcor_i[1] ; rcor_j = mcor_j[1]
            } else {
              rcor_i = mcor_i[2] ; rcor_j = mcor_j[2]
            }

            rcorr = mean(c(abs(rcor_i),abs(rcor_j)))

            #if significant correlation ...
            if ( abs(rcor_j) > corr_thr || abs(rcor_i) > corr_thr ) {

              if(case==1) {
                decision = paste0("Case_1: correlated. Keep ID [",ID_i,"] and remove ID [",ID_j,"].")
                to_remove = append(to_remove, ID_j)
                if(is.null(tmp)) {
                  rtab_i[,"decision"] = sub("^NA / ","",paste(rtab_i[,"decision"],"/",decision))
                  tmp = rtab_i
                } else {
                  tmp[which(tmp$ID==m),"decision"] = sub("^NA / ","",paste(tmp[which(tmp$ID==m),"decision"],"/",decision))
                }

              } else if(case==2) {
                decision = paste0("Case 2: correlated. Keep ID [",ID_j,"] and remove ID [",ID_i,"].")
                to_remove = append(to_remove, ID_i)
                rtab_j[,"decision"] = sub("^NA / ","",paste(rtab_j[,"decision"],"/",decision))
                if(is.null(tmp)) {
                  tmp = rtab_j
                } else {
                  tmp = rbind(tmp,rtab_j)
                }

              } else if(case==3) {
                if(!is.null(tmp)) { if(nrow(tmp)>1) { merged_IDs = c(merged_IDs, c(ID_i,ID_j)) } }
                if((end_i-start_i) > (end_j-start_j)) { MDS_ij = MDS_i } else { MDS_ij = MDS_j }
                decision = paste0("Case 3: correlated. Merge coordinates of ID [",ID_i,"]+[",ID_j,"] and keep longest MDS [",MDS_ij,"].")
                new_region = data.frame(ID=paste0(ID_i,",",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_j, noutliers=length(union(outliers_i,outliers_j)), SizeRegions=round((end_j-start_i)/1000000,1), decision=decision)
                if(is.null(tmp) || nrow(tmp)==1) { tmp = new_region } else { tmp[which(tmp$ID==m),] = new_region }

              } else if(case==4) {
                if(!is.null(tmp)) { if(nrow(tmp)>1) { merged_IDs = c(merged_IDs, c(ID_i,ID_j)) } }
                if((end_i-start_i) > (end_j-start_j)) { MDS_ij = MDS_i } else { MDS_ij = MDS_j }
                decision = paste0("Case 4: correlated. Merge coordinates of ID [",ID_j,"]+[",ID_i,"] and keep longest MDS [",MDS_ij,"].")
                new_region = data.frame(ID=paste0(ID_i,",",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_i, noutliers=length(union(outliers_i,outliers_j)), SizeRegions=round((end_i-start_j)/1000000,1), decision=decision)
                if(is.null(tmp) || nrow(tmp)==1) { tmp = new_region } else { tmp[m,] = new_region }
              }

              rtab[n,"decision"] = sub("^NA / ","",paste(rtab[n,"decision"],"/",decision))
              list_merge = append(list_merge, decision)
              if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nrtab_i\n")) ; print(rtab_i) ; cat("rtab_j\n") ; print(rtab_j) ; cat(paste0("\n",decision,"\n\nto_remove\n")) ; print(to_remove) ; cat("\ntmp\n") ; print(tmp)  }


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
                if(MDS_i == MDS_j) { rcor_ij = mcor_ij[1] } else { rcor_ij = mcor_ij[2] }

                if( abs(rcor_ij) > corr_thr ) {
                  if(!is.null(tmp)) { if(nrow(tmp)>1) { intersect_IDs = c(intersect_IDs, c(ID_i,ID_j)) } }

                  if(case==1) {
                    decision = paste0("Case 1: intersect correlated between [",ID_j,"] & [",ID_i,"]. Keep longest MDS [",MDS_i,"].")
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_j, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=rtab_j[,"SizeRegions"], decision=decision)

                  } else if(case==2) {
                    decision = paste0("Case 2: intersect correlated between [",ID_j,"] & [",ID_i,"]. Keep longest MDS [",MDS_j,"].")
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_i, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=rtab_i[,"SizeRegions"], decision=decision)

                  } else if(case==3) {
                    if((end_i-start_i) > (end_j-start_j)) { MDS_ij = MDS_i } else { MDS_ij = MDS_j }
                    decision = paste0("Case 3: intersect correlated between [",ID_j,"] & [",ID_i,"]. Keep longest MDS [",MDS_ij,"].")
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_j, end=end_i, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=round((end_i-start_j)/1000000,1), decision=decision)

                  } else if(case==4) {
                    if((end_i-start_i) > (end_j-start_j)) { MDS_ij = MDS_i } else { MDS_ij = MDS_j }
                    decision = paste0("Case 4: intersect correlated between [",ID_j,"] & [",ID_i,"]. Keep longest MDS [",MDS_ij,"].")
                    new_region = data.frame(ID=paste0(ID_i,"/",ID_j), MDS=MDS_ij, chr=chr_j, start=start_i, end=end_j, noutliers=length(intersect(outliers_i,outliers_j)), SizeRegions=round((end_j-start_i)/1000000,1), decision=decision)
                  }

                  tmp = rbind(tmp, new_region, rtab_j)
                  if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nrtab_i\n")) ; print(rtab_i) ; cat("rtab_j\n") ; print(rtab_j) ; cat(paste0("\n",decision,"\n")) ; cat("\ntmp\n") ; print(tmp) }

                } else {
                  decision = paste0("Case ",case,": intersect > ",clustsize," windows between [",ID_i,"] and [",ID_j,"]. No correlation.")
                  rtab[n,"decision"] = sub("^NA / ","",paste(rtab[n,"decision"],"/",decision))
                  if(is.null(tmp)) {
                    rtab[n-1,"decision"] = sub("^NA / ","",paste(rtab[n-1,"decision"],"/",decision))
                    tmp = rtab[(n-1):n,]
                  } else if(!(rtab[n,"ID"] %in% to_remove)) {
                    tmp[which(tmp$ID==m),"decision"] = sub("^NA / ","",paste(tmp[which(tmp$ID==m),"decision"],"/",decision))
                    tmp = rbind(tmp, rtab[n,])
                  }
                  if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nrtab_i\n")) ; print(rtab_i) ; cat("rtab_j\n") ; print(rtab_j) ; cat(paste0("\n",decision,"\n\nto_remove\n")) ; print(to_remove) ; cat("\ntmp\n") ; print(tmp) }
                }

              } else {
                decision = paste0("Case ",case,": intersect < ",clustsize," windows between [",ID_i,"] and [",ID_j,"].")
                rtab[n,"decision"] = sub("^NA / ","",paste(rtab[n,"decision"],"/",decision))
                if(is.null(tmp)) {
                  rtab[n-1,"decision"] = sub("^NA / ","",paste(rtab[n-1,"decision"],"/",decision))
                  tmp = rtab[(n-1):n,]
                } else if(!(rtab[n,"ID"] %in% to_remove)) {
                  tmp[which(tmp$ID==m),"decision"] = sub("^NA / ","",paste(tmp[which(tmp$ID==m),"decision"],"/",decision))
                  tmp = rbind(tmp, rtab[n,])
                }
                if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nrtab_i\n")) ; print(rtab_i) ; cat("rtab_j\n") ; print(rtab_j) ; cat(paste0("\n",decision,"\n\nto_remove\n")) ; print(to_remove) ; cat("\ntmp\n") ; print(tmp) }
              }
            }

          } else if(do_clust==FALSE) {

            if(is.character(tmp[,"ID"])) { tmp_IDs = as.numeric(unlist(strsplit(tmp[,"ID"],",|/"))) } else { tmp_IDs = tmp[,"ID"] }

            if(!rtab_j["ID"] %in% tmp_IDs) {

              if(!is.null(tmp)) {

                if(is.character(m)) { m_IDs = as.numeric(unlist(strsplit(m,",|/"))) } else { m_IDs = m }
                m_tmp = tmp[grep(paste0("\\b",m_IDs,"\\b",collapse="|"),tmp[,"ID"]),]

                # if non-overlapping regions are included in other merged regions that are overlapping the current region (rtab_j);
                # remove non-overlapping regions from tmp
                if(nrow(m_tmp[which(m_tmp$end >= rtab_j$start),]) > 0) {
                  to_remove = append(to_remove, m_tmp[which(m_tmp$end < rtab_j$start),"ID"])
                  # otherwise keep longest one if several regions are overlapping each other but not the current region (rtab_j)
                } else {
                  m_tmp = m_tmp[which(m_tmp$end < rtab_j$start),]
                  m_tmp_IDs = m_tmp[,"ID"]
                  if(nrow(m_tmp)>1) { m_tmp = m_tmp[which(m_tmp$SizeRegions==max(m_tmp$SizeRegions)),] }
                  if(nrow(m_tmp)>1) { m_tmp = m_tmp[which(m_tmp$noutliers==max(m_tmp$noutliers)),] }
                  if(nrow(m_tmp)>1) { m_tmp = m_tmp[1,] }
                  to_remove = append(to_remove, m_tmp_IDs[which(m_tmp_IDs!=m_tmp[,"ID"])])
                  rm(m_tmp_IDs)
                }

                if(nrow(m_tmp)==1) {

                  count = count + 1

                  decision = paste0("No more redundant outlier region.")
                  m_tmp[,"decision"] = sub("^NA / ","",paste(m_tmp[,"decision"],"/",decision))
                  tabF = rbind(tabF, cbind(nID=count, m_tmp))
                  listF_merge[[count]] = list_merge
                  list_merge = list()
                  if(m==tail(rtab_i_IDs,1)) { tmp = NULL } else { tmp = tmp[which(tmp$ID != m),] }
                  if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nm_tmp\n")) ; print(m_tmp) ; cat("rtab_j\n") ; print(rtab_j) ; cat("\ntabF\n") ; print(tabF) ; cat("\ntmp\n") ; print(tmp) }
                } else if (nrow(m_tmp)==0) {
                  print(paste0("[",n,",",m,"] !! ISSUE : tmp is empty !!"))
                }
                rm(m_tmp)

              } else {

                count = count + 1

                listF_merge[[count]] = "No redundant outlier region."
                rtab[n-1,"decision"] = "No redundant outlier region."
                tabF = rbind(tabF, cbind(nID=count, rtab[n-1,]))
                if(isTRUE(verbose)) { cat(paste0("n = ",n,"\nrtab_i\n")) ; print(rtab_i) ; cat("rtab_j\n") ; print(rtab_j) ; cat("\nNo redundant outlier region.\n\ntabF\n") ; print(tabF) ; cat("\ntmp\n") ; print(tmp) }
              }

            }

          }

          #if different chromosome
        } else {

          count = count + 1

          if(!is.null(tmp)) {
            tabF = rbind(tabF, cbind(nID=count, rtab_i))
            listF_merge[[count]] = list_merge
            list_merge = list()
            if(m==tail(rtab_i_IDs,1)) { tmp = NULL } #something's missing ???

          } else {
            listF_merge[[count]] = "No redundant outlier region."
            rtab[n-1,"decision"] = "No redundant outlier region."
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
            for(dupid in dup_IDs) {
              decision = tmp[which(tmp$ID==dupid),"decision"]
              tmp[which(tmp$ID==dupid),"decision"] = decision[nchar(decision)==max(nchar(decision))]
            }
            tmp = tmp[!duplicated(tmp[,"ID"]),]
          }

          # remove region that got added to tmp if it has been later merged to another region in tmp
          if(is.character(tmp[,"ID"])) {
            if(any(duplicated(unlist(strsplit(tmp[,"ID"],","))))) {
              dup_IDs = unique(unlist(strsplit(tmp[,"ID"],","))[duplicated(unlist(strsplit(tmp[,"ID"],",")))])
              tmp = tmp[which(!as.character(tmp[,"ID"]) %in% dup_IDs),]
              if(nrow(tmp[grep(paste0("\\b",dup_IDs,"\\b",collapse="|"), tmp[,"ID"]),])>0) {
                if(isTRUE(verbose)) { cat(paste0("COMMENT : Removed region(s) [",dup_IDs,"] from tmp because it was/were found merged with other region(s) in tmp. \nNew tmp:\n")) ; print(tmp) }
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
              noutliers = length(unlist(listRegionH1[as.numeric(unique(unlist(strsplit(dup_rows[,"ID"],",|/"))))]))
              decision = paste0("Case 3/4: correlated. Merge coordinates of duplicated ID [",dup_IDs,"] and keep longest MDS [",MDS_ij,"].")
              new_region = data.frame(ID=dup_IDs, MDS=MDS_ij, chr=chr_j, start=min(dup_rows$start), end=max(dup_rows$end), noutliers=noutliers, SizeRegions=round((max(dup_rows$end)-min(dup_rows$start))/1000000,1), decision=decision)
            }
            tmp = rbind(new_region,tmp[!grepl(paste0("\\b",merged_IDs[duplicated(merged_IDs)],"\\b",collapse="|"), tmp[,"ID"]),])
            if(isTRUE(verbose)) { cat("Duplicated merged regions found !\ndup_rows:\n") ; print(dup_rows) ; cat("\nBest MDS:\n") ; print(MDS_ij) ; cat("\nNew tmp:\n") ; print(tmp) }
          }

        } else {
          tmp = NULL
        }

      }

    } else {
      to_remove = list()
    }

  }

  # when reaching last region, wrap up remaining regions.
  if(n==nrow(rtab)) {

    count = count + 1

    if(!is.null(tmp)) { tmp[,"decision"] = unlist(lapply(tmp[,"decision"],function(x) paste0(x," / End of table."))) }

    if(is.character(tmp[,"ID"])) { tmp_IDs = unique(as.numeric(unlist(strsplit(tmp[,"ID"],",|/")))) } else { tmp_IDs = tmp[,"ID"] }
    if(is.character(tabF[,"ID"])) { tabF_IDs = unique(as.numeric(unlist(strsplit(tabF[,"ID"],",|/")))) } else { tabF_IDs = tabF[,"ID"] }

    if(any(!rtab_j$ID %in% tmp_IDs) && any(!rtab_j$ID %in% tabF_IDs)) {
      rtab[n,"decision"] = "End of table, no redundant region."
      tmp = rbind(tmp,rtab[n,])
    }

    tabF = rbind(tabF, cbind(nID=count, tmp))
    list_merge = append(list_merge, "End of table.")
    listF_merge[[count]] = list_merge
  }


  return(list(rtab=rtab,tabF=tabF,listF_merge=listF_merge))

}


