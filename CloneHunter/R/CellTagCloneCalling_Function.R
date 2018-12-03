#' Single-cell RNA-seq Whitelisting Function
#'
#' This function conducts whitelist filtering through the single-cell dataset
#' @param celltag.dat raw celltag single-cell data matrix
#' @param tag.cutoff How many tags would you like to be used as a cutoff to say that the cells are tagged?
#' @return Binarized single-cell dataset
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' SingleCellDataWhitelist(sc.mtx, 2)
#' 
SingleCellDataBinarization <- function(celltag.dat, tag.cutoff) {
  # Filter out cells with only 1 CellTag and binarize the data
  # 0 - have less than 2 celltags
  # 1 - have 2 or more celltags
  CellTags <- celltag.dat
  CellTags[CellTags < tag.cutoff] <- 0
  CellTags[CellTags > 0] <- 1
  
  return(CellTags)
}
#' Single-cell RNA-seq Whitelisting Function
#'
#' This function conducts whitelist filtering through the single-cell dataset
#' @param celltag.dat binarized celltag single-cell data matrix
#' @param whitels.cell.tag.file file director to the whitelisted cell tags
#' @return Whitelisted single-cell dataset
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' SingleCellDataWhitelist(sc.mtx, "~/Desktop/My_Favourite_Whitelist.csv", 2)
#' 
SingleCellDataWhitelist <- function(celltag.dat, whitels.cell.tag.file) {
  # Store the cell names
  CellTags <- celltag.dat
  cell.names <- rownames(celltag.dat)
  
  # Process the celltag matrix to format below
  # row - celltag
  # col - cells
  CellTags <- t(CellTags)
  CellTags <- as.data.frame(CellTags)
  celltag.rownames <- row.names(CellTags)
  
  # Filter the matrix using whitelist
  if (endsWith(whitels.cell.tag.file, ".csv")) {
    separator <- ","
  } else {
    if (endsWith(whitels.cell.tag.file, ".txt") | endsWith(whitelisted.cell.tag, ".tsv")) {
      separator <- "\t"
    } else {
      separator <- " "
    }
  }
  whitelist <- read.delim(whitels.cell.tag.file, sep = separator, header = T, stringsAsFactors = F)
  whitelist.names <- whitelist[,1]
  whitelist <- Reduce(intersect, list(whitelist.names, celltag.rownames))
  celltags.whitelisted <- CellTags[whitelist,]
  colnames(celltags.whitelisted) <- cell.names
  
  return(as.data.frame(t(celltags.whitelisted)))
}

#' CellTag Metric Plotting Function
#'
#' This function provides some metric plots for further downstream celltag filtering in the scRNA-seq dataset
#' @param celltag.data The single-cell data with celltags
#' @return A list contain the average occurrence of each celltag across cells and the average frequency of different celltags in each cell
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricPlots(celltags.whitelisted)
#'
MetricPlots <- function(celltag.data) {
  CellTags.per.cell.whitelisted.pf <- rowSums(celltag.data)
  CellTags.per.cell.avg <- mean(CellTags.per.cell.whitelisted.pf)
  CellTags.frequency.whitelisted.pf <- colSums(celltag.data)
  CellTags.freq.avg <- mean(CellTags.frequency.whitelisted.pf)
  plot(CellTags.per.cell.whitelisted.pf)
  plot(CellTags.frequency.whitelisted.pf)
  return(list(CellTags.per.cell.avg, CellTags.freq.avg))
}

#' Metric-Base Filtering Function
#'
#' This function applies further filtering on scRNA-seq data with CellTags based on cutoff values identified from the metric plots
#' @param whitelisted.celltag.data The single-cell data with whitelisted celltag.
#' @param cutoff The cutoff decided from the metric plots
#' @param comparison Would you like to maintain the part less than/greater than the cutoff? Default to less. Choices can be greater or less.
#' @return A new scRNA-seq data frame with CellTags that passes the metric-based filtering
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' MetricBasedFiltering(celltags.whitelisted, 20, "less")
#'
MetricBasedFiltering <- function(whitelisted.celltag.data, cutoff, comparison = "less") {
  # Set up the filtering data frame
  CellTags.per.cell.whitelisted.pf <- as.data.frame(rowSums(whitelisted.celltag.data))
  
  # Set up the filtered celltag dataset object
  if (comparison == "less") {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf < (cutoff + 1))
  } else {
    cell.filter <- subset(CellTags.per.cell.whitelisted.pf, CellTags.per.cell.whitelisted.pf > (cutoff - 1))
  }
  cell.bc.filter <- row.names(cell.filter)
  # Filter celltag dataset
  celltags.whitelisted <- as.data.frame(t(whitelisted.celltag.data))
  celltags.whitelisted.new <- as.data.frame(t(celltags.whitelisted[cell.bc.filter]))
  
  return(celltags.whitelisted.new)
}

#' Jaccard Analysis Function
#'
#' This function conducts Jaccard analysis to calculate the Jaccard similarity between cells
#' @param whitelisted.celltag.data The single-cell data with whitelisted celltag.
#' @param save.mtx Would you like to save the jaccard matrix? Default to TRUE. Save the matrix to the working directory.
#' @param plot.corr Would you like to plot the correlation matrix?
#' @param save.plot Would you like to save the plot? Save the plot to the working directory.
#' @return The jaccard matrix
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' JaccardAnalysis(celltags.whitelisted)
#'
JaccardAnalysis <- function(whitelisted.celltag.data, save.mtx = TRUE, plot.corr = TRUE, save.plot = TRUE) {
  # Calculating the Jaccard matrix
  Jac <- simil(whitelisted.celltag.data, method = "Jaccard")
  Jac <- as.matrix(Jac)
  if (save.mtx) {
    saveRDS(Jac, paste0(getwd(), "/Jaccard_mtx.RDS"))
  }
  
  if (plot.corr) {
    diag(Jac) <- 1
    p1 <- corrplot(Jac, method="color", order="hclust", hclust.method ="ward.D2", cl.lim=c(0,1), tl.cex=0.1)
    if (save.plot) {
      pdf(paste0(getwd(), "/Jaccard_correlation_plot.pdf"), width = 20, height = 20, paper = "special")
      print(p1)
      dev.off()
    } else {
      p1
    }
  }
  return(Jac)
}

#' Clone Calling Function
#'
#' This function conducts clone calling based on the Jaccard results
#' @param Jaccard.Matrix The Jaccard matrix
#' @param output.dir Which directory would you like to save your clone table information?
#' @param output.filename Which CSV file would you like to save the clone table into?
#' @param correlation.cutoff Correlation cutoff for clone membership
#' @return A list of clone table with membership of each cell and clone size counting table.
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CloneCalling(Jac, "~/Desktop/", "My_Favourite_CloneCalling.csv", 0.7)
#'
CloneCalling <- function(Jaccard.Matrix, output.dir, output.filename, correlation.cutoff) {
  # Using the igraph package to facilitate the identification of membership to each clone
  test <- Jaccard.Matrix*lower.tri(Jaccard.Matrix)
  check.corelation <- which(test > correlation.cutoff, arr.ind=TRUE)
  graph.cor <- graph.data.frame(check.corelation, directed = FALSE)
  groups.cor <- split(unique(as.vector(check.corelation)), clusters(graph.cor)$membership)
  conv.groups.cor <- lapply(groups.cor,
                            function(list.cor){
                              rownames(test)[list.cor]}
                            )
  
  # Put clones into tables
  l <- seq(1, length(groups.cor))
  df.conv <- apply(as.matrix(l), 1, 
                   function(x) {
                     data.frame(clone.id = x, 
                                cell.barcode = conv.groups.cor[[x]], 
                                stringsAsFactors = FALSE)
                     }
                   )
  
  df.comb <- rbindlist(df.conv)
  write.csv(df.comb, paste0(output.dir, output.filename), row.names = F, quote = F)
  
  # Calculate the size of each clone
  counts <- table(df.comb$clone.id)
  counts <- as.data.frame(counts)
  colnames(counts) <- c("Clone.ID", "Frequency")
  
  return(list(df.comb, counts))
}



