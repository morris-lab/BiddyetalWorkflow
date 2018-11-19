


returnDirectlyConnectedNodes <- function(node, linkList){
  tmp_link <- linkList[linkList$source %in% node,]
  tmp_link2 <- linkList[linkList$target %in% node,]
  
  tmp_nodes <- union(tmp_link$target, tmp_link2$source)
  tmp_nodes <- union(tmp_nodes, node)
  return(tmp_nodes)
}



returnAllConnectedNodes <- function(node, linkList){
  for (i in 1:5) {
    node <- returnDirectlyConnectedNodes(node, linkList)
  }
  return(node)
}


drawNetworkGraph <- function(linkList, Nodes, overlay){

  rownames(Nodes) <- 1:nrow(Nodes)
  
  ref <- 1:nrow(Nodes)
  names(ref) <- Nodes$nodes
  linkList$source1 <- ref[as.character(linkList$source)] - 1
  linkList$target1 <- ref[as.character(linkList$target)] - 1
  
  linkList$Value <- 1
  #linkList$Colour <- c("#CD6155", "#566573")[as.numeric(linkList[,3] > 0) + 1]
  
  a <- forceNetwork(Links = linkList, Nodes = Nodes, zoom = T,opacityNoHover = 0.5, 
                    Source = "source1", Target = "target1", arrows = T,
                    NodeID = "nodes", Value ="Value" , #linkColour = linkList$Colour,
                    Group = overlay, opacity = 0.9)
  
  return(a)
  
}


drawSubnet <- function(tag, overlay, linkList, Nodes){
  # e.g. tag; "celltag2.1_698"
  # e.g. color: "cluster" or "tag" or "SuperClone"    
  
  no <- returnAllConnectedNodes(tag, linkList)
  sub_link <- linkList[(linkList$source %in% no) | (linkList$target %in% no),]
  sub_Nodes <- Nodes[Nodes$nodes %in% no ,]
  
  a <- drawNetworkGraph(sub_link, sub_Nodes, overlay)
  
  return(a)
}