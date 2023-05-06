#Jianing Wang
#Boston University Dissertation
#Chapter 2

#####################################################
##################### Functions #####################
## Remove the List with the Smallest Catchability ###
#####################################################

remove_smallest_list <- function(dtYfullTable, nDatasets_val){
  if(nDatasets_val == 5){
    dtYfullTable1 <- subset(dtYfullTable, select = - X5)
    ### Recalculate the obs label
    for(i in 1:nrow(dtYfullTable1)){
      dtYfullTable1$obs_label1[i] <- ifelse(sum(dtYfullTable1$X1[i],dtYfullTable1$X2[i],
                                                         dtYfullTable1$X3[i],dtYfullTable1$X4[i])>0,1,0)
    }
    dtYfullTable1 <- dtYfullTable1[, c("X1","X2","X3","X4","target_label","obs_label1","loc_label")]
    colnames(dtYfullTable1) <- c("X1","X2","X3","X4","target_label","obs_label","loc_label")
  }
  if(nDatasets_val == 7){
    dtYfullTable1 <- subset(dtYfullTable_one_simu, select = - X7)
    ### Recalculate the obs label
    for(i in 1:nrow(dtYfullTable1)){
      dtYfullTable1$obs_label1[i] <- ifelse(sum(dtYfullTable1$X1[i],dtYfullTable1$X2[i],
                                                         dtYfullTable1$X3[i],dtYfullTable1$X4[i],
                                                         dtYfullTable1$X5[i],dtYfullTable1$X6[i])>0,1,0)
    }
    dtYfullTable1 <- dtYfullTable1[, c("X1","X2","X3","X4","X5","X6","target_label","obs_label1","loc_label")]
    colnames(dtYfullTable1) <- c("X1","X2","X3","X4","X5","X6","target_label","obs_label","loc_label")
  }
  return(dtYfullTable1)
}