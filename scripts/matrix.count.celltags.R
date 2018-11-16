
library(data.table)

commandArguments <- commandArgs(trailingOnly = TRUE)

parseCommands <- function(argumentList){
  
  parsedArgs <- list(summary = argumentList[1], samFile = argumentList[2], prefix = argumentList[3])
  
}

#Creates a list of arguments passed on the command line.
parsedArgs <- parseCommands(commandArguments)

#Creates a file name for the summary file.
summaryOut <- paste0(parsedArgs$prefix, ".celltag.stats.txt")

#Creates a file name for the CellTag UMI Count Matrix.
matrixOut <- paste0(parsedArgs$prefix, ".matrix.tsv")

#Creates a file name to save the CellTag UMI Count Matrix as an R object. 
rdsOut <- paste0(parsedArgs$prefix, ".celltag.matrix.Rds")

#Read in Summary File which contains a list of the Cell Barcodes identified during alignment.
barcodeList <- fread(parsedArgs$summary, header = FALSE)[[1]]

#Make sure the list of Cell Barcodes only contains unique barcodes.
barcodeList <- unique(barcodeList)

#Split Cell Barcodes into Cell Barcode and suffix number.
barcodeList <- strsplit(barcodeList, split = "-", fixed = TRUE)

#Select the Barcode sequence from the list of split barcodes, discarding the trailing
barcodeList <- sapply(barcodeList, "[[", 1)

#Read in "SAM" File containing Read.ID, Read.Seq, Cell.BC, UMI, Cell.Tag, and Gene if aligned. 
celltagData <- fread(parsedArgs$samFile)

#With the parsed CellTag reads loaded we can then easily filter the data and generate UMI Counts for each Cell Barcode/Cell Tag combination.

celltagCounts <- celltagData[, .(Cell.BC, Cell.Tag, UMI)                          #This filters the data.table to just Cell Barcodes, Cell Tags, and UMIs. 
                           ][, .(UMI.Count = uniqueN(UMI)), .(Cell.BC, Cell.Tag)  #This groups the data.table by Cell Barcode/Cell Tag combination and creates a new column "UMI.Count" which has the number of unique UMI associated with each Cell Barcode/Cell.Tag combination. uniqueN is equivalent to length(unique(UMI))
                           ][ !Cell.Tag == "",
                            ]                                                     #This final operation removes any rows which do not contain a cell tag.


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





