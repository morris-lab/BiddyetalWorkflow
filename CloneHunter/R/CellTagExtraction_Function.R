#' CellTag Extraction Function
#'
#' This function extracts CellTags from the raw fastq sequencing file, provides counts of each CellTag and sorts them in desending order.
#' @param fastq.bam.input The input fastq/bam data directory
#' @param short.nt.before.tag A short sequence before the 8nt tag to help more specific identification
#' @param short.nt.after.tag A short sequence after the 8nt tag to help more specific identification
#' @param extraction.output.filename The output file directory, i.e. where would you like to store the 8nt tags with their counts?
#' @param save.fullTag Would you like to save full sequences with the short sequences provided without counts? Defaults to TRUE
#' @param save.onlyTag Would you like to save only the 8nt CellTags without counts? Defaults to TRUE
#' @return A list contains count table of CellTags. If requested to save fullTag counts, i.e. save.fullTag.counts = TRUE, return a list of both 8nt tags and full sequences count. Otherwise, a list of 8nt tags counts. 
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' CellTagExtraction("data.fastq", "v1", "~/Desktop/whitelist.txt")
#' 
CellTagExtraction <- function(fastq.bam.input, celltag.version, extraction.output.filename, save.fullTag = TRUE, save.onlyTag = TRUE) {
  # Set up the output filenames and path
  output.file <- extraction.output.filename
  output.file.dirname <- dirname(output.file)
  output.file.basename <- basename(output.file)
  output.file.bnm.parts <- strsplit(output.file.basename, "[.]")[[1]]
  output.file.bnm <- paste0(output.file.bnm.parts[c(1:(length(output.file.bnm.parts) - 1))], collapse = "_")
  
  if (tolower(celltag.version) == "v1") {
    short.nt.before.tag <- "CCGGT"
  } else if (tolower(celltag.version) == "v2") {
    short.nt.before.tag <- "GTGATG"
  } else {
    short.nt.before.tag <- "TGTACG"
  } 
  short.nt.after.tag <- "GAATTC"
  # Set up the patterns to extract the tags and open the connection to file
  pattern <- paste0(short.nt.before.tag, "[ATCG]{8}", short.nt.after.tag)
  if (endsWith(fastq.bam.input, "fastq")) {
    con <- file(fastq.bam.input, "r")
    
    # Get the sequences containing the tags (with both full tag region and only 8nt tag)
    seq.list <- c()
    filtered.sequences <- c()
    full.tag.seq <- c()
    only.tag.seq <- c()
    print("Reading File......")
    while(TRUE) {
      curr.lines <- readLines(con, 1000000)
      if (length(curr.lines) == 0) break
      else {
        curr.seqs <- curr.lines[seq(2, 1000000, by = 4)]
        seq.list <- c(seq.list, curr.seqs)
        reg.rslt <- regexpr(pattern, curr.seqs, ignore.case = TRUE, perl = TRUE)
        contain.idx <- which(reg.rslt > 0)
        curr.f.seq <- curr.seqs[contain.idx]
        
        filtered.sequences <- c(filtered.sequences, curr.f.seq)
        start.loc <- reg.rslt[contain.idx]
        end.loc <- start.loc + nchar(short.nt.before.tag) + 8 + nchar(short.nt.after.tag) - 1
        curr.full.tag <- substr(curr.f.seq, start = start.loc, stop = end.loc)
        only.tag <- substr(curr.full.tag, start = (nchar(short.nt.before.tag) + 1), stop = (nchar(short.nt.before.tag) + 8))
        full.tag.seq <- c(full.tag.seq, curr.full.tag)
        only.tag.seq <- c(only.tag.seq, only.tag)
      }
    }
    close(con)
  }
  if (endsWith(fastq.bam.input, "bam")) {
    #########################
    # Under construction
    #########################
    # Install Rsamtools
    if (!requireNamespace("BiocManager", quietly = TRUE))
      install.packages("BiocManager")
    if (!requireNamespace("Rsamtools", quietly = TRUE)) {
      BiocManager::install("Rsamtools", version = "3.8")
      require(Rsamtools)
    }
    # Get the bam file
    bamFile <- BamFile(fastq.bam.input)
    yieldSize(bamFile) <- 1000000
    open(bamFile)
    while(TRUE) {
      curr.seqs <- as.character(scanBam(bamFile)[[1]]$seq)
      if (length(curr.seqs) == 0) break
      else {
        seq.list <- c(seq.list, curr.seqs)
        reg.rslt <- regexpr(pattern, curr.seqs, ignore.case = TRUE, perl = TRUE)
        contain.idx <- which(reg.rslt > 0)
        curr.f.seq <- curr.seqs[contain.idx]
        
        filtered.sequences <- c(filtered.sequences, curr.f.seq)
        start.loc <- reg.rslt[contain.idx]
        end.loc <- start.loc + nchar(short.nt.before.tag) + 8 + nchar(short.nt.after.tag) - 1
        curr.full.tag <- substr(curr.f.seq, start = start.loc, stop = end.loc)
        only.tag <- substr(curr.full.tag, start = (nchar(short.nt.before.tag) + 1), stop = (nchar(short.nt.before.tag) + 8))
        full.tag.seq <- c(full.tag.seq, curr.full.tag)
        only.tag.seq <- c(only.tag.seq, only.tag)
      }
    }
  }
  
  # Save them if required
  if (save.fullTag) {
    print("Saving Full CellTag......")
    fullTag.nm <- paste0(output.file.dirname, "/", output.file.bnm, "_fullTag.txt")
    write.table(full.tag.seq, fullTag.nm, sep = "\n", quote = F, row.names = F, col.names = F)
  }
  if (save.onlyTag) {
    print("Saving CellTag......")
    onlyTag.nm <- paste0(output.file.dirname, "/", output.file.bnm, "_onlyTag.txt")
    write.table(only.tag.seq, onlyTag.nm, sep = "\n", quote = F, row.names = F, col.names = F)
  }
  
  rslt <- list(full.tag.seq, only.tag.seq)
  return(rslt)
}

#' CellTag Matrix Generation Function
#'
#' This function uses the extract information from data processed before and generate a Cell Barcode x CellTag matrix
#' @param barcodes.file A .tsv output file from 10x CellRanger pipeline. It contains a list of all cell barcodes identified in the filtered dataset.
#' @param celltag.read.data.file A .tsv output file that contains information about the parsed CellTag read data
#' @return Cell Barcodes x CellTag matrix
#' @keywords single-cell RNA-seq data, CellTagging
#' @export
#' @examples
#' ("data.fastq", "v1", "~/Desktop/whitelist.txt")
#' 
CellTagMatrixCount <- function(barcodes.file, celltag.read.data.file) {
  # Get the basename of this file
  base.fnm.parts <- strsplit(basename(barcodes.file), "[.]")[[1]]
  base.fnm <- base.fnm.parts[c(1:(length(base.fnm.parts) - 1))]
  dir.nm <- dirname(celltag.read.data.file)
  
  #Creates a file name for the summary file.
  summaryOut <- paste0(dir.nm, "/", base.fnm, ".celltag.stats.txt")
  #Creates a file name for the CellTag UMI Count Matrix.
  matrixOut <- paste0(dir.nm, "/", base.fnm, ".matrix.tsv")
  #Creates a file name to save the CellTag UMI Count Matrix as an R object. 
  rdsOut <- paste0(dir.nm, "/", base.fnm, ".celltag.matrix.Rds")
  
  # Read in the cell barcodes identified during alignment
  barcodeList <- fread(barcodes.file, header = FALSE)[[1]]
  #Make sure the list of Cell Barcodes only contains unique barcodes.
  barcodeList <- unique(barcodeList)
  #Split Cell Barcodes into Cell Barcode and suffix number.
  barcodeList <- strsplit(barcodeList, split = "-", fixed = TRUE)
  #Select the Barcode sequence from the list of split barcodes, discarding the trailing
  barcodeList <- sapply(barcodeList, "[[", 1)
  
  #Read in "SAM" File containing Read.ID, Read.Seq, Cell.BC, UMI, Cell.Tag, and Gene if aligned. 
  celltagData <- fread(celltag.read.data.file)
  #With the parsed CellTag reads loaded we can then easily filter the data and generate UMI Counts for each Cell Barcode/Cell Tag combination.
  #-Filters the data.table to just Cell Barcodes, Cell Tags, and UMIs. 
  #-Groups the data.table by Cell Barcode/Cell Tag combination and creates a new column "UMI.Count" which has the number of unique UMI associated with each Cell Barcode/Cell.Tag combination. uniqueN is equivalent to length(unique(UMI))
  #-Removes any rows which do not contain a cell tag.
  celltagCounts <- celltagData[, c("Cell.BC", "UMI", "Cell.Tag")][, .(UMI.Count = uniqueN(UMI)), .(Cell.BC, Cell.Tag)][ !Cell.Tag == "",]           
  #The data is now in a long format and needs to be reshaped. We will cast the long data into a wide format resembling a matrix.
  celltagCountsWide <- dcast(data = celltagCounts, formula = Cell.BC ~ Cell.Tag, value.var = "UMI.Count", fill = 0 )
  
  #Now we have the data we want in the correct format. Next we can add Cells from the barcode list that were not in the celltagData.
  missingCells <- barcodeList[!(barcodeList %in% celltagCountsWide$Cell.BC)]
  #Lets make a data.table with one column Cell.BC which will contain a list of the missing cells. This can then be merged with the UMI Count data table.
  missingCells <- data.table(Cell.BC = missingCells)
  #Bind the missing cells to the data.table containing the Cell Tag UMI Counts.
  alltagCounts <- rbind(celltagCountsWide, missingCells, fill = TRUE)
  #We have added the missing cells, whose values now need to be changed from NA to 0.
  alltagCounts[is.na(alltagCounts)] <- 0
  #Now we can filter out cells which are not in our barcode list.
  alltagCounts <- alltagCounts[Cell.BC %in% barcodeList, ]
  
  #Lets also filter Cell Tags in which no UMIs are counted.
  celltagExpr <- colSums(alltagCounts[, -1])
  tagsRemove <- names(celltagExpr)[celltagExpr == 0]
  alltagCounts[, (tagsRemove):= NULL]
  
  #We now have a final matrix. Next lets generate some stats about the Cell Tags. 
  celltagExpr <- summary(colSums(alltagCounts[, -1]))
  cellsPerTag <- summary(colSums(alltagCounts[, -1] > 0))
  cellExpr <- summary(rowSums(alltagCounts[, -1]))
  
  tagsPerCell <- rowSums(alltagCounts[, -1] > 0)
  tagsPerCellSum <- summary(tagsPerCell)
  
  stats.df <- rbind(celltagExpr, cellsPerTag, cellExpr, tagsPerCellSum)
  rownames(stats.df) <- c("CellTag.UMI.Counts", "Cells.per.CellTag", "Cell.UMI.Counts", "CellTags.per.Cell")
  stats.df <- as.data.frame(stats.df)
  
  sink(summaryOut)
  cat("Number of Cells expressing No Cell Tags\n")
  cat(sum(tagsPerCell == 0))
  cat("\nNumer of Cells expressing 1 or more CellTags\n")
  cat(sum(tagsPerCell >= 1))
  cat("\nNumber of Cells Expressing 5 or more CellTags\n")
  cat(sum(tagsPerCell >= 5))
  cat("\nNumber of Cells Expressing 10 or more CellTags\n")
  cat(sum(tagsPerCell >= 10))
  cat("\n\n")
  
  print(stats.df)
  sink()
  
  fwrite(alltagCounts, file = matrixOut)
  saveRDS(alltagCounts, file = rdsOut)
  
  return(alltagCounts)
}

