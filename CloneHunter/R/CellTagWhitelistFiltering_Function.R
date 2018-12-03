#' CellTag Whitelist Filtering Function
#'
#' This function conducts whitelist filtering such that only CellTags with count number over their certain percentile would be considered for clone calling
#' @param count.sorted.table The count table for all CellTags with first column being CellTag and second being counts.
#' @param percentile A fraction cutoff percentile for filtering the CellTags e.g. 0.9 for 90th percentile
#' @param output.dir Which directory would you like to store these files?
#' @param output.count.file Which CSV filename would you like to save to for the filtered count table? Default to NULL. If NULL but save.output = TRUE, automatic names will be generated and used for storage
#' @param save.output Would you like to save your whitelisted counts to file? Default to TRUE
#' @return Whitelisted CellTag count table
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagWhitelistFiltering("data.txt", 0.9, "~/Desktop/", "My_Favourite_Whitelist.csv")
#' 
CellTagWhitelistFiltering <- function(count.sorted.table, percentile, output.dir, 
                                      output.count.file, save.output = TRUE) {
  # Load table and calculate cutoff
  count.cutoff <- quantile(count.sorted.table$Count, probs = percentile)
  count.true.cut <- floor(count.cutoff/10)
  
  # Name of the file
  bnm <- basename(output.count.file)
  bnm.parts <- strsplit(bnm, "[.]")[[1]]
  bnm.only <- paste(bnm.parts[1:(length(bnm.parts) - 1)], collapse = "_")
  
  # Plot
  if (!endsWith(output.dir, "/")) output.dir <- paste0(output.dir, "/")
  quantile.dist.plot.nm <- paste0(output.dir, bnm.only, "_", percentile, "_percentile.pdf")
  pdf(quantile.dist.plot.nm, width = 6, height = 6, paper = "special")
  plot(count.sorted.table$Count, main="CellTag Whitelist",xlab="CellTag",ylab="Reads")
  abline(v=sum(count.sorted.table$Count >= count.true.cut), col="red", lty=2)
  dev.off()
  print(paste0("Abline Threshold: ", sum(count.sorted.table$Count >= count.true.cut)))

  # Subset the ones pass filtering
  whitelist <- subset(count.sorted.table, Count>=count.true.cut)
  
  # Save the file
  if (save.output) {
    write.csv(whitelist, paste0(output.dir, output.count.file), 
              quote = F, row.names = F)
  }
  
  return(whitelist)
}
