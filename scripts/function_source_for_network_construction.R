##  R functions for cellTag network construction and visualization ##
# Please see tutorials as well.
# updated at 2018, 11, 14
#


convertCellTagMatrix2LinkList <- function(celltag_data){
  # celltag_data should be data frame (N x 3).
  # the columnname of this data frame should be c("CellTagV1", "CellTagV2", "CellTagV3")
  
  
  ### 1.Preprocessing celltag data #### 
  message("Preprocessing data..")
   # pick up cells that have one or more celltag, and remove cells that do not have any celltag.
  Cells_with_tag <- rownames(celltag_data)[!(is.na(celltag_data$CellTagV1) & 
                                               is.na(celltag_data$CellTagV2) & 
                                               is.na(celltag_data$CellTagV3))]
  

  message(paste0(" Cells that have CellTagV1: ", sum(!is.na(celltag_data$CellTagV1))))
  message(paste0(" Cells that have CellTagV2: ", sum(!is.na(celltag_data$CellTagV2))))
  message(paste0(" Cells that have CellTagV3: ", sum(!is.na(celltag_data$CellTagV3))))
  
  
  # remove non-tagged cells
  celltag_data <- celltag_data[Cells_with_tag, ]
  
  # convert NA to "e"
  tags <-  c("CellTagV1", "CellTagV2", "CellTagV3")
  for (i in tags) {
    celltag_data[is.na(celltag_data[ ,i]),i] <- "e"}
  
  ### 2. Constructing LinkList ###
  message("Constructing link list..")
  
  findRoot <- function(cell_id, tag) { # e.g, cell_id = "TGTTCCGGTGAGGCTA-8"; tag = "CellTagV1", "CellTagV2", "CellTagV3"
    tagid <- celltag_data[cell_id,tag]
    tmp <- as.data.frame(t(c(paste0(tag, "_", tagid), cell_id, tag)), stringsAsFactors = F)
    rownames(tmp) <- NULL
    colnames(tmp) <- c("source", "target", "tag")
    return(tmp)
  }
  

  ## first, clonal population that share the same celltag is combined to make subnetwork.
  ## then, subnewtorks will be combined further if they are originated from same mother. 

  all_cell_id <- rownames(celltag_data)
  remaining_cell_id <- all_cell_id
  tags <- c("CellTagV3", "CellTagV2", "CellTagV1")
  linkList <- data.frame()
  
  # 2.1 find connection between "celltag" -> "cells"
  for (tag in tags) {
    remaining_cells <- celltag_data[remaining_cell_id,]
    subcells <- remaining_cells[remaining_cells[,tag] != "e",]
    
    tmp <- foreach(i = rownames(subcells), .combine = rbind) %do% {
      findRoot(i, tag)
    }
    linkList <- rbind(linkList, tmp)
    done_id <- rownames(subcells)
    # remaining_cell_id <- remaining_cell_id[!(remaining_cell_id %in% done_id)]   !!!!algorithm was modified 20180830!!!! in new version, remaining_cell_id will now be updated.
  }
  
  
  # 2.2 hidden link ["CellTagV2" -> "CellTagV3"], or ["CellTagV1" -> "CellTagV3"]
  hiddenlink_D13 <- foreach(i = (unique(celltag_data$CellTagV3)[-1]), .combine = rbind) %do% {
    
    sub_cells <- celltag_data[celltag_data$CellTagV3 == i, ]
    
    prev_tag <- sub_cells$CellTagV2
    prev_tag <- prev_tag[prev_tag != "e"]
    prev_tag <- names(which.max(table(prev_tag)))
    
    if (class(prev_tag) != "NULL") {
      tmp <- as.data.frame(t(c(paste0("CellTagV2", "_", prev_tag),
                               paste0("CellTagV3", "_", i),
                               "CellTagV2")), stringsAsFactors = F)
      rownames(tmp) <- NULL
      colnames(tmp) <- c("source", "target", "tag")
      return(tmp)
    } 
    
    prev_tag <- sub_cells$CellTagV1
    prev_tag <- prev_tag[prev_tag != "e"]
    prev_tag <- names(which.max(table(prev_tag)))
    
    if (class(prev_tag) != "NULL") {
      tmp <- as.data.frame(t(c(paste0("CellTagV1", "_", prev_tag),
                               paste0("CellTagV3", "_", i),
                               "CellTagV1")), stringsAsFactors = F)
      rownames(tmp) <- NULL
      colnames(tmp) <- c("source", "target", "tag")
      return(tmp)
    }   
    
  }
  # 2.3 hidden link ["CellTagV1" -> "CellTagV2"]
  hiddenlink_D3 <- foreach(i = (unique(celltag_data$CellTagV2)[-1]), .combine = rbind) %do% {
    
    sub_cells <- celltag_data[celltag_data$CellTagV2 == i, ]
    
    prev_tag <- sub_cells$CellTagV1
    prev_tag <- prev_tag[prev_tag != "e"]
    prev_tag <- names(which.max(table(prev_tag)))
    
    if (class(prev_tag) != "NULL") {
      tmp <- as.data.frame(t(c(paste0("CellTagV1", "_", prev_tag),
                               paste0("CellTagV2", "_", i),
                               "CellTagV1")), stringsAsFactors = F)
      rownames(tmp) <- NULL
      colnames(tmp) <- c("source", "target", "tag")
      return(tmp)
    }   
    
  }
  rm(remaining_cells, remaining_cell_id, sub_cells, subcells, all_cell_id, done_id, prev_tag, tag, tags)
  
  # 2.4 integrating all links
  
  modifyCellName <- function(linkList){
    # this function change cell name,  like.. "TTCTCCTGTATCACCA-7"  -> "TTCTCCTGTATCACCA-7_D3"
    # in the date processing algorithm v-0.20, cells that have multiple cell tag will show up mutiple times.
    # thus we have make new name to avoid them being overrapped. 
    
    linkList$target_unmodified <- linkList$target
    
    node_cell <- grep("-", linkList$target)
    
    linkList[node_cell, "target"] <- paste0(linkList[node_cell, "target"], 
                                            "_",
                                            stringr::str_split_fixed(linkList[node_cell, "tag"], "g", 2)[,2])
    return(linkList)
  }
  
  # integrate
  linkList <- rbind(linkList, hiddenlink_D3)
  linkList <- rbind(linkList, hiddenlink_D13)
  
  # change cell name
  linkList <- modifyCellName(linkList)
  
  message("finished")
  return(linkList)

}


getNodesfromLinkList <- function(linkList){
  # This function construct Nodes list from linkList.
  # Use "convertCellTagMatrix2LinkList" function before running this function to get linkList.
  
  nodes <- union(linkList$target, linkList$source)
  Nodes <- data.frame(nodes, row.names = nodes, stringsAsFactors = F)
  
  
  #tag
  refferTagid <- function(each_node) {
    cells_or_not <- (sum(c("CellTagV1", "CellTagV2", "CellTagV3") %in% strsplit(each_node, "_")[[1]]) == 0) 
    
    if (cells_or_not) {
      ans <- linkList[linkList$target == each_node, "tag"]
    } else {
      ans <- strsplit(each_node, "_")[[1]][1]
    }
    return(ans)
  }
  
  refferUMname <- function(each_node){
    cells_or_not <-  (sum(c("CellTagV1", "CellTagV2", "CellTagV3") %in% strsplit(each_node, "_")[[1]]) == 0) 
    
    
    if (cells_or_not) {
      ans <- linkList[linkList$target == each_node, "target_unmodified"]
    } else {
      ans <- each_node
    }
    return(ans)
    
  }
  
  
  
  Nodes$tag <- sapply(nodes, refferTagid)
  Nodes$node_name_unmodified <- sapply(nodes, refferUMname)
  
  return(Nodes)
}


addData2Nodes <- function(Nodes, additional_data){
  
  # Nodes: data frame
  # additional_data: data frame
  #
  # the rownames of additional_data should be same format as "node_name_unmodified" in Nedes
  
  return(cbind(Nodes, additional_data[Nodes$node_name_unmodified,]))
}

