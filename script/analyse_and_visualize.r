library(tidyverse)
library(ggembl)
library(patchwork)
library(ggrepel)

if (rstudioapi::isAvailable()) {
  print("We are within Rstudio, so mounted.")
  wdir <- dirname(rstudioapi::getActiveDocumentContext()$path)
}else{
  print("We are on the server side, probably.")
  wdir <- getwd()
}
setwd(wdir)

# PARAMS
plotWindows <- FALSE
experimentName <- "CV001N"
#usedStrain <- "atcc_8492_concatenated"

for (sample in gsub(pattern = "[.]fastq", replacement = "", x = list.files(str_c("../data/experiments/", experimentName, "/raw")))) {
  depthInfo <- read_tsv(str_c("../results/out/", experimentName, "/depth_", sample, ".tab"), col_names = F)
  depthInfo$X1 <- NULL
  colnames(depthInfo) <- c("Position", "coverage")
  
  usedStrain <- str_c(tail(str_split(sample, pattern = "_")[[1]], 4), collapse = "_")
  #sample <- str_replace(string = sample, pattern = usedStrain)
  
  # Join genome nucleotides
  # if (sample == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample2" || sample == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3"){
  #   print("Loading B. vulgatus genome...")
  #   nucs <- read_tsv("../data/genomes_long_format/B_vulgatus_genome_long_format.txt", col_names = F)
  # } else {
  #   nucs <- read_tsv("../data/genomes_long_format/B_uniformis_genome_long_format.txt", col_names = F)
  # }
  nucs <- read_tsv(str_c("../data/genomes_long_format/", usedStrain, "_long.txt"), col_names = F)
  
  
  nucs <- nucs %>% rename(nucleotide = X1) %>% mutate(Position = 1:dim(nucs)[1])
  nucs$nucleotide <- map_chr(nucs$nucleotide, str_to_upper)
  depthInfo <- left_join(depthInfo, nucs, by = 'Position')
  TA <- map2_lgl(depthInfo$nucleotide[1:(length(depthInfo$nucleotide) - 1)], 
                 depthInfo$nucleotide[2:(length(depthInfo$nucleotide))], function(x, y) ifelse(x == "T" && y == "A", TRUE, FALSE))
  depthInfo$TA <- c(TA, T)
  
  # Join insertion info file (basically parsed samtools view containing SAM flag, position, CIGAR flag, and sequence) and parse it further:
  readInsertion <- read_tsv(str_c("../results/out/", experimentName, "/insertionInfoRaw_", sample, ".tab"), col_names = F)
  # Keep only reads that MATCH PERFECTLY (using CIGAR flag).
  # This doesn't generalize to all experiments..
  readInsertion <- readInsertion[str_ends(string = readInsertion$X5, pattern = "M") & str_length(readInsertion$X5) == 3, ]
  # Rename columns
  colnames(readInsertion) <- c("readName", "SAM", "Position", "MAPQ", "CIGAR", "read", "readLength")
  readInsertion$readLength <- as.numeric(readInsertion$readLength)
  # Filter reads based on MAPQ score
  readInsertion <- readInsertion %>% filter(MAPQ == 42)
  
  # Read barcode mapping file and left join it to readInsertion
  # reads_and_their_barcodes <- read_tsv(str_c("results/CV001I10_all_forward_reads_concatenated_140000_reads_per_sample_readname_barcode.tab"), col_names = F)
  reads_and_their_barcodes <- read_tsv(str_c("../results/out/", experimentName, "/readname_barcode_", sample, ".tab"), col_names = F)
  colnames(reads_and_their_barcodes) <- c("readName", "barcode")
  n <- dim(readInsertion)[1]
  readInsertion <- left_join(readInsertion, reads_and_their_barcodes, by = 'readName')
  stopifnot(dim(readInsertion)[1] == n)
  
  # IMPORTANT: Here I remove all reads where barcode == NA
  # I'm not checking for barcode integrity when I filter reads but this file (_sequence_readname_barcode) does have a check for barcode integrity.
  # Thus, we expect to lose some reads here..
  # This can happen when the IR/bridging site could not be identified (probably rare) OR when at least one pjhred score over the barcode was < 30.
  readInsertion <- readInsertion %>% filter(!is.na(barcode))
  
  # The position of reverse-mapped reads (SAM flag == 16) needs to be parsed to get the insertion position. See  some lines below, the following command
  # readInsertion %>% filter(InsertionPosition == 231)
  readInsertion$InsertionPosition <- pmap_dbl(list(readInsertion$SAM, readInsertion$Position, readInsertion$readLength), function(x, y, z){
    if (x == 16){
      return(y + z - 2)
    } else {
      return(y)
    }
  })
  readInsertion$InsertionPosition <- as.integer(readInsertion$InsertionPosition)
  # We've got everything, and things look good
  # observe InsertionPosition == 231 for sample HMW5KAFX2_CV001L3_21s001573-1-1_Voogdt_lane1Sample1_sequence to convince yourself that the parsing works.
  #readInsertion %>% filter(InsertionPosition == 231)
  
  # Remove read info, don't need it anymore
  readInsertion$read <- NULL
  
  # Observe that we have 19 mil reads left (for Sample1)
  dim(readInsertion)
  
  # This function returns the mode/'purity' of barcodes over a given position.
  # Importantly, it does not include reads that did not have a barcode!
  # That means that if we have barcodes in a group like so: c("A", "A", "B", NA, NA),
  # We'd get back A and 0.66 for purity sa opposed to 0.4 for purity!!!
  Mode <- function(x, returnPurity = F) {
    #print(x)
    if (any(is.na(x))){
      x <- x[which(!is.na(x))]
    }
    #print(x)
    if (length(x) == 0){
      return(NA)
    }
    #print(x)
    ux <- unique(x)
    if (returnPurity) {
      return(table(x)[which.max(table(x))] / length(x))
    } else {
      return(ux[which.max(tabulate(match(x, ux)))]) 
    }
  }
  
  # Important: at this point, I have joined readInsertion file with their barcodes. 
  # Thus, here we have the barcode information of all well-mapping mapping reads.
  # We can now also compute the distribution of barcode maps over the entire genome. 
  # So after this, we'll have two complementary tables that will give us information:
  ### 1. The mode/purity of barcodes over insertion sites (computed below)
  ### 2. The distribution of barcodes over insertion sites (computed here)
  
  # IMPORTANT: This is NOT a final map of barcodes -> insertion sites yet. To achieve this,
  # We need to remove barcodes that are 'good' for more than one position (this can happen in multiple integration events!)
  
  readInsertionWithBarcodes <- readInsertion
  
  barcodesOverGenome <- readInsertion %>% filter(!is.na(barcode)) %>% group_by(barcode) %>% count(InsertionPosition) %>% arrange(desc(n)) %>% nest() %>% mutate(goodBarcode = map(data, function(x){
    # Count InsertionPositions and sort so that the most dominant position is in row 1, second-most dominantn in row 2 and so forth.
    modeInsertionPosition <- x$InsertionPosition[1]
    InsertionsInmostDominantPosition <- x$n[1]
    fractionInsertionsInmostDominantPosition <- x$n[1]/sum(x$n)
    if (dim(x)[1] == 1){
      # Only a single insertion site was hit. 
      fractionInsertionsInsecondDominantPosition <- NA
    } else {
      fractionInsertionsInsecondDominantPosition <- x$n[2]/sum(x$n)
    }
    
    if (InsertionsInmostDominantPosition > 10 && fractionInsertionsInmostDominantPosition > (3/4) && (is.na(fractionInsertionsInsecondDominantPosition) || fractionInsertionsInsecondDominantPosition < (1/8))){
      return(list(TRUE, x$InsertionPosition[1], InsertionsInmostDominantPosition, fractionInsertionsInmostDominantPosition, fractionInsertionsInsecondDominantPosition))
    } else {
      return(list(FALSE, NA, InsertionsInmostDominantPosition, fractionInsertionsInmostDominantPosition, fractionInsertionsInsecondDominantPosition))
    }
  })) %>% 
    mutate(InsertionPosition = map_dbl(goodBarcode, function(x) x[[2]])) %>% 
    mutate(InsertionsInmostDominantPosition = map_dbl(goodBarcode, function(x) x[[3]])) %>%
    mutate(fractionInsertionsInmostDominantPosition = map_dbl(goodBarcode, function(x) x[[4]])) %>%
    mutate(fractionInsertionsInsecondDominantPosition = map_dbl(goodBarcode, function(x) x[[5]])) %>%
    mutate(goodBarcode = map_lgl(goodBarcode, function(x) x[[1]]))
  
  
  # Group by insertionPosition and get the number of (perfectly matching) reads that inserted here.
  # Also, compute the mode barcode and the barcode purity at that position.
  readInsertion <- readInsertion %>% group_by(InsertionPosition) %>% summarize(numberInsertions = length(readLength),
                                                                               modeBarCode = Mode(barcode, returnPurity = FALSE),
                                                                               barCodeModePurity = Mode(barcode, returnPurity = TRUE))
  depthInfo <- left_join(depthInfo, readInsertion, by = c('Position' = "InsertionPosition"))
  
  # since we left-join on the depthInfo, we have to replace NAs with 0s (readInsertion table will have no entries where not a single read maps)
  depthInfo$numberInsertions <- map_dbl(depthInfo$numberInsertions, function(x) ifelse(is.na(x), 0, x))
  depthInfo$TAIndex <- NA
  depthInfo$TAIndex[depthInfo$TA] <- 1:sum(depthInfo$TA)
  
  # We're done with parsing the insertion information now.
  
  # Join GFF
  expand <- function(start, end, description, geneName){
    s <- start[1]
    e <- end[1]
    return(data.frame(Position = s:e,
                      info = description[1],
                      geneName = geneName[1]))
  }
  gff <- read_tsv(str_c("../data/prokka/", usedStrain, ".gff"), col_names = F, comment = "#")
  gff <- gff[,c("X3", "X4", "X5", "X7", "X9")]
  colnames(gff) <- c("type", "start", "end", "strandedness", "description")
  # Runs a couple of minutes
  gff$geneName <- map_chr(gff$description, function(x) {
    tmp <- str_split(string = x, pattern = "gene=")[[1]][2]
    tmp <- str_split(string = tmp, pattern = ";")[[1]][1]
  })
  
  gffExpanded <- pmap(list(gff$start, gff$end, gff$description, gff$geneName), expand)
  gffExpanded <- do.call('rbind', gffExpanded)
  
  # When an ORF overlaps, we will retain both positions here. For simplicity, remove a random one
  depthInfo <- left_join(depthInfo, gffExpanded, by = 'Position')
  depthInfo <- depthInfo %>% distinct(Position, .keep_all = T)
  depthInfo$isOrf <- !is.na(depthInfo$info)
  
  print("Mean Coverage: ")
  mean(depthInfo$coverage)
  
  print("Around 10% of the genome is not covered:")
  dim(depthInfo)
  dim(depthInfo %>% filter(coverage == 0))
  
  show <- function(df, mm){
    p <- ggplot()
    p <- p + geom_line(data = df %>% filter(Position > mm[1] & Position < mm[2]), 
                       aes(x = Position, y=log10(coverage+1), group = 1))
    p3 <- ggplot() + geom_point(data = df %>% filter(Position > mm[1] & Position < mm[2]),
                              aes(x = Position, y = barCodeModePurity), stat = 'identity')
    p3 <- p3 + ylim(c(0, 1))
    p <- p + theme_embl()
    p <- p + geom_line(data = df %>% filter(Position > mm[1] & Position < mm[2]) %>%mutate(notCovered = coverage == 0), 
                       aes(x = Position, y = log10(5000), color = notCovered, group = 1))
    p2 <- ggplot() + geom_line(data = df %>% filter(Position > mm[1] & Position < mm[2]) %>% rename(isORF = isOrf), 
                               aes(x = Position, y = log10(10), color = isORF, group = 1), size = 5)
    p2 <- p2 + theme_embl() + theme(axis.ticks = element_blank(),
                                    axis.title = element_blank(),
                                    axis.text = element_blank())
    tmp <- df %>% filter(Position > mm[1] & Position < mm[2]) %>% filter(!is.na(geneName)) %>% group_by(geneName) %>% summarize(meanPos = mean(Position))
    tmp$y <- 10
    p2 <- p2 + geom_text_repel(data = tmp, aes(x = meanPos, y = log10(y), label = geneName), min.segment.length = 0, max.overlaps = Inf) 
    p2 <- p2 + ggtitle(str_c("Coords: ", mm[1], " to ", mm[2]))
    return(p2/p3/p +  plot_layout(heights = c(1, 1, 4)))
  }
  if (plotWindows){
    #tmp <- df %>% filter(Position > mm[1] & Position < mm[2]) %>% filter(!is.na(geneName)) %>% group_by(geneName) %>% summarize(meanPos = mean(Position))
    coord <- 1
    windowSize <- 50000
    dir.create(str_c("plots/", sample), recursive = T)
    pdf(str_c("plots/", sample, "/", "coverage_chunked_with_gene_annotation.pdf"), width = 7.5, height = 4.75)
    while (TRUE){
      #print(str_c('Processing ', coord))
      minMax <- c(coord, coord + windowSize)
      print(show(depthInfo, minMax))
      coord <- coord + windowSize
      if (coord > max(depthInfo$Position)){
        break
      }
    }
    dev.off()
  }

  
  ##################################################
  ### This is the point where we parsed everything #
  ##################################################
  
  # We do have some insertions into non-TA dinucleotide regions. A total of 16k, a total of xx for coverage more than yy.
  # We should probably check what's going on here. Can it be that the transposon sometimes inserts in non-TA sites?
  depthInfo %>% filter(!TA) %>% filter(numberInsertions > 0) 
  
  # Let's analyse stuff!
  
  # Number Of TAs in genome.
  print(str_c("Total number of TAs in genome: ", sum(depthInfo$TA)))
  print(str_c("Out of those, ", 
              depthInfo %>% filter(TA) %>% filter(numberInsertions == 0) %>% dim() %>% magrittr::extract(1), " TAs are never hit with a single read."))
  print(str_c("Out of those, ", 
              depthInfo %>% filter(TA) %>% filter(numberInsertions <= 2) %>% dim() %>% magrittr::extract(1), " TAs are hit less than three times."))
  print(str_c("There are a total of ", depthInfo %>% filter(!is.na(info)) %>% pull(info) %>% unique() %>% length(), " ORFs in the B. uniformis genome. Actually probably a few less (a few of those might be rRNA genes etc.)"))
  print("Let's look at the number of TAs per ORF")
  TAsPerORF <- depthInfo %>% filter(!is.na(info)) %>% group_by(info) %>% summarize(numberTAs = sum(TA))
  hist(TAsPerORF$numberTAs)
  print("Some orfs have a lot of TAs. Probably correlated with size and stuff? whatever...")
  print(str_c("How many ORFs have 0 TAs?: ", TAsPerORF %>% filter(numberTAs == 0) %>% dim() %>% magrittr::extract(1)))
  
  # Define function to define essential genes.
  getEssen <- function(Position, numberInsertions){
    # Trim edge TAs
    if (any(Position != sort(Position))){
      print("Position isnt sorted. Exiting.")
      sadsadd
    }
    TAPosIndex <- 1:length(Position)
    # Trim of (roughly) the leftmost and rightmost 20% of TAs.
    ccL <- quantile(TAPosIndex, 0.2)
    ccR <- quantile(TAPosIndex, 0.8)
    numberInsertions <- numberInsertions[TAPosIndex > ccL & TAPosIndex < ccR]
    if (length(numberInsertions) == 0){
      return(NA)
    }
    if (all(numberInsertions <= 3)){
      return(TRUE)
    } else {
      return(FALSE)
    }
  }
  
  
  
  print("Let's have a look at the average insertion numbers (over all TAs) over all ORFs. This way, we can define essential ORFs")
  # Some of those metrics need to be generated again!!!
  ORFEssentiallityUsingTAs <- depthInfo %>% filter(!is.na(info)) %>% group_by(info) %>% summarize(meanNumberInsertionsOverTAs = mean(numberInsertions), # This is WRONGLY NAMED HURR DURR.
                                                                                                  numberTAs = sum(TA), 
                                                                                                  meanCoverage = mean(coverage), 
                                                                                                  totalInsertions = sum(numberInsertions),
                                                                                                  TARichness = sum(TA)/length(info), 
                                                                                                  geneLength = length(info),
                                                                                                  numberInsertionsByTASites = sum(numberInsertions) / sum(TA), # This is the same as meanNumberInsertionsOverTAs
                                                                                                  uniqueInsertionsByTASites = sum(numberInsertions > 0 & TA)/ sum(TA),
                                                                                                  uniqueInsertions = sum(numberInsertions > 0 & TA))
  
  # for essential gene inferrence, I want to look only at TA sites
  tmp <- depthInfo %>% filter(!is.na(info)) %>% filter(TA) %>% ungroup() %>% group_by(info) %>% summarise(isEssential = getEssen(Position, numberInsertions))
  ORFEssentiallityUsingTAs <- left_join(ORFEssentiallityUsingTAs, tmp, by = 'info')
  
  # tmp %>% filter(isEssential.Sample1 | isEssential.Sample2) %>% dim() %>% magrittr::extract(1)
  # tmp %>% filter(! isEssential.Sample1 & isEssential.Sample2) %>% dim() %>% magrittr::extract(1)
  # tmp %>% filter(isEssential.Sample1 & ! isEssential.Sample2) %>% dim() %>% magrittr::extract(1)
  # tmp %>% filter(isEssential.Sample1 & isEssential.Sample2) %>% dim() %>% magrittr::extract(1)
  
  
  print(str_c("Number of ORFs without any insertions: ",  ORFEssentiallityUsingTAs %>% filter(meanNumberInsertionsOverTAs == 0) %>% dim() %>% magrittr::extract(1)))
  print(str_c("Number of ORFs without any insertions in central TAs: ",  ORFEssentiallityUsingTAs %>% filter(isEssential) %>% dim() %>% magrittr::extract(1)))
  print(str_c("Unique insertions not in ORFs: ", depthInfo %>% filter(!isOrf) %>% pull(TA) %>% sum()))
  
  # Reproduction of Fig S1 of https://mbio.asm.org/content/9/1/e02096-17#DC1
  #hist(ORFEssentiallityUsingTAs$uniqueInsertionsByTASites, 1000)
  
  # Our ad-hoc definition distinguishes the 'two distributions' quite well.
  p <- ggplot(ORFEssentiallityUsingTAs, aes(x = uniqueInsertionsByTASites, group = isEssential, fill = isEssential))
  p <- p + geom_histogram(position = 'identity', alpha = 0.3, bins = 300)
  
  ORFEssentiallityUsingTAs <- ORFEssentiallityUsingTAs %>% arrange(meanNumberInsertionsOverTAs)
  ORFEssentiallityUsingTAs$geneName <- map_chr(ORFEssentiallityUsingTAs$info, function(x) {
    tmp <- str_split(string = x, pattern = "gene=")[[1]][2]
    tmp <- str_split(string = tmp, pattern = ";")[[1]][1]
  })
  ORFEssentiallityUsingTAs$prokkaID <- map_chr(ORFEssentiallityUsingTAs$info, function(x) {
    tmp <- str_split(string = x, pattern = "ID=")[[1]][2]
    tmp <- str_split(string = tmp, pattern = ";")[[1]][1]
  })
  # This is really nice!
  #View(ORFEssentiallityUsingTAs)
  
  # Write out essential genes...
  #write_tsv(x = ORFEssentiallityUsingTAs %>% 
  #            filter(isEssential) %>% 
  #            pull(info) %>% 
  #            map_chr(function(x) str_split(string = x, pattern = "ID=")[[1]][2]) %>% map_chr(function(x) str_split(string = x, pattern = ';')[[1]][1]) %>%
  #            as.data.frame(), path = str_c("B_uniformis_essential_genes_", sample,".tab"), col_names = F)
  
  # And this is also cool to see (mean coverage over ORF and mean insertion over ORFs correlates very well!)
  plot(ORFEssentiallityUsingTAs$meanNumberInsertionsOverTAs, ORFEssentiallityUsingTAs$meanCoverage)
  
  
  # Get QC plots
  if (!dir.exists(str_c("../plots/", experimentName, "/"))){
    dir.create(str_c("../plots/", experimentName, "/"))
  }
  if (sample == "000000000-JLF44_CV001N_21s002101-1-1_Voogdt_lane1Sample3_B_vulgatus_atcc_8482"){
    expectedNumberOfInsertions <- 96
  } else {
    expectedNumberOfInsertions <- 190
  }
  p <- ggplot(depthInfo %>% filter(numberInsertions > 0) %>% arrange(desc(numberInsertions)), aes(x = 1:length(numberInsertions), y = numberInsertions, group = 1))
  p <- p + scale_y_log10() + geom_line() + theme_embl() + geom_vline(xintercept = expectedNumberOfInsertions)
  ggsave(plot = p, filename = str_c("../plots/", experimentName, "/", sample, "__numberInsertions.png"), width = 6, height = 4)
  #plots[[length(plots) + 1]] <- p
  
  # Get also total barcode counts ... 
  tmp <- readInsertionWithBarcodes %>% group_by(barcode) %>% nest() %>% mutate(numBarcodes = map_dbl(data, function(x) dim(x)[1])) %>% arrange(desc(numBarcodes))
  #p <- ggplot(tmp, aes(x = numBarcodes)) + geom_histogram(bins = 100)
  p <- ggplot(tmp, aes(x = numBarcodes)) + geom_histogram(bins = 100) + scale_x_log10()
  ggsave(plot = p, filename = str_c("../plots/", experimentName, "/", sample, "__barCodeCount.png"), width = 6, height = 4)
  
  # ... and also an ordered lineplot
  p <- ggplot(tmp, aes(x = 1:length(numBarcodes), y = numBarcodes, group - 1)) + geom_line() + scale_y_log10()
  p <- p + geom_vline(xintercept = expectedNumberOfInsertions)
  p <- p + scale_x_log10()
  ggsave(plot = p, filename = str_c("../plots/", experimentName, "/", sample, "__barCodeCountAsLinePlot.png"), width = 6, height = 4)
  
  # And finally the last plot: For the [expectedNumberOfInsertions] top-most represented barcodes, count 
  # in how many unique sites you find it!
  
  # OLD
  # tmp <- tmp %>% arrange(desc(numBarcodes)) %>% head(96) %>% mutate(uniqueInsertionEventsCountedByBarcode = map_dbl(data, function(x) {
  #   #print(x)
  #   x <- x %>% count(InsertionPosition) %>% filter(n > 20) %>% dim() %>% magrittr::extract2(1)
  #   return(x)
  # }))
  
  # NEW
  tmp <- tmp %>% arrange(desc(numBarcodes)) %>% head(expectedNumberOfInsertions) %>% mutate(uniqueInsertionEventsCountedByBarcode = map(data, function(x) {
    #print(x)
    x <- x %>% count(InsertionPosition)
    return(x)
  }))
  results <- list()
  for (i in c(1,2,3,4,5, 10, 20, 30, 40, 50, 75, 100, 200, 500)){
    print(i)
    tmpp <- tmp %>% mutate(uniqueInsertionEventsCountedByBarcode = map_dbl(uniqueInsertionEventsCountedByBarcode, function(x){
      x %>% filter(n > i) %>% dim() %>% magrittr::extract2(1)
    })) %>% mutate(threshold = i) %>% select(threshold, uniqueInsertionEventsCountedByBarcode)
    results[[length(results) + 1]] <- tmpp
  }
  results <- do.call('rbind', results)
  
  p <- ggplot(results %>% group_by(threshold) %>% count(uniqueInsertionEventsCountedByBarcode), aes(x = uniqueInsertionEventsCountedByBarcode, y = n)) + geom_bar(stat = 'identity')
  p <- p + facet_wrap(.~threshold,
                      nrow = length(unique(results$threshold)), scales = 'free')
  p <- p + scale_x_continuous(breaks = 0:10)
  p <- p + ylab("Count")
  ggsave(filename = str_c("../plots/", experimentName, "/", sample, "__distinctInsertionEventsPerBarcode.png"), width = 6, height = 12)
  
  # Save files for manual interrogation later, clean up and go on.
  save(sample, depthInfo, ORFEssentiallityUsingTAs, barcodesOverGenome, readInsertionWithBarcodes, file = str_c("../results/out/", experimentName, "/save_file", sample, ".rsave"))
  rm(depthInfo)
  rm(ORFEssentiallityUsingTAs)
  rm(barcodesOverGenome)
  rm(readInsertionWithBarcodes)
  gc()
}


##################################
# HERE COMES SOME CODE FOR LATER #
##################################

# # Compare afonso's p-values to our essentially definition.
# tmp <- map(results, function(x) x[[2]] %>% select(prokkaID, numberInsertionsByTASites, totalInsertions, uniqueInsertions, isEssential))
# tmp[[1]]$type <- "our_agar"
# tmp[[2]]$type <- "our_liquid"
# afonso_agar <- read_csv("/home/nicolai/for_carlos_4/Bacteroides_uniformis_ATCC_8492_alldomains_agar.csv", comment = "#")
# afonso_liquid <- read_csv("/home/nicolai/for_carlos_4/Bacteroides_uniformis_ATCC_8492_alldomains_liquid.csv", comment = "#")
# tmp[[1]] <- left_join(tmp[[1]], afonso_agar %>% select(`Total Tn5 insertions`, `Essentiality p-value`, `Gene ID`, Essentiality), by = c("prokkaID" = "Gene ID"))
# tmp[[1]] <- left_join(tmp[[1]], afonso_liquid %>% select(`Total Tn5 insertions`, `Essentiality p-value`, `Gene ID`, Essentiality), by = c("prokkaID" = "Gene ID"), 
#                       suffix = c("_afonso_agar", "_afonso_liquid"))
# a <- tmp[[1]] %>% group_by(prokkaID) %>% summarize(uniqueInsertions = sum(uniqueInsertions), `Total Tn5 insertions` = sum(`Total Tn5 insertions_afonso_agar`, na.rm = T))
# cor(a$uniqueInsertions, a$`Total Tn5 insertions`)
# a <- tmp[[1]] %>% group_by(prokkaID) %>% summarize(uniqueInsertions = sum(uniqueInsertions), `Total Tn5 insertions` = sum(`Total Tn5 insertions_afonso_liquid`, na.rm = T))
# cor(a$totalInsertions, a$`Total Tn5 insertions`)
# 
# tmp[[2]] <- left_join(tmp[[2]], afonso_agar %>% select(`Total Tn5 insertions`, `Essentiality p-value`, `Gene ID`, Essentiality), by = c("prokkaID" = "Gene ID"))
# tmp[[2]] <- left_join(tmp[[2]], afonso_liquid %>% select(`Total Tn5 insertions`, `Essentiality p-value`, `Gene ID`, Essentiality), by = c("prokkaID" = "Gene ID"), 
#                       suffix = c("_afonso_agar", "_afonso_liquid"))
# a <- tmp[[2]] %>% group_by(prokkaID) %>% summarize(totalInsertions = sum(totalInsertions), `Total Tn5 insertions` = sum(`Total Tn5 insertions_afonso_agar`, na.rm = T))
# cor(a$totalInsertions, a$`Total Tn5 insertions`)
# a <- tmp[[2]] %>% group_by(prokkaID) %>% summarize(totalInsertions = sum(totalInsertions), `Total Tn5 insertions` = sum(`Total Tn5 insertions_afonso_liquid`, na.rm = T))
# cor(a$totalInsertions, a$`Total Tn5 insertions`)
# 
# 
# # Generate Figure 2A/B/C from 
# # https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5099075/
# 
# # A
# tmp <- map(results, function(x) x[[2]] %>% select(prokkaID, numberInsertionsByTASites, meanNumberInsertionsOverTAs, isEssential))
# tmp[[1]]$type <- "our_agar"
# tmp[[2]]$type <- "our_liquid"
# tmp <- full_join(tmp[[1]] %>% select(prokkaID, meanNumberInsertionsOverTAs, type), tmp[[2]] %>% select(prokkaID, meanNumberInsertionsOverTAs, type))
# tmp <- tmp %>% pivot_wider(id_cols = prokkaID, names_from = type, values_from = meanNumberInsertionsOverTAs)
# 
# A <- ggplot(tmp, aes(x = log10(our_agar + 1), y = log10(our_liquid + 1))) + geom_point()
# 
# # B
# tmp <- map(results, function(x) x[[2]] %>% select(prokkaID, numberInsertionsByTASites, totalInsertions, isEssential))
# tmp[[1]]$type <- "our_agar"
# tmp[[2]]$type <- "our_liquid"
# tmp <- full_join(tmp[[1]] %>% select(prokkaID, numberInsertionsByTASites, type), tmp[[2]] %>% select(prokkaID, numberInsertionsByTASites, type))
# #tmp <- tmp %>% pivot_wider(id_cols = prokkaID, names_from = type, values_from = totalInsertions)
# 
# B <- ggplot(tmp, aes(x = numberInsertionsByTASites, fill = type)) + geom_histogram(alpha = 0.3, position = 'identity', bins = 200)
# B <- B + xlab("Mean number of\ninsertion events per ORF")
# B1 <- ggplot(tmp, aes(x = numberInsertionsByTASites, fill = type)) + geom_histogram(alpha = 0.3, position = 'identity', bins = 50) + xlim(c(0, 100))
# B1 <- B1 + xlab("Mean number of\ninsertion events per ORF")
# 
# # C (might be a bit more complicated)
# tmp <- map(results, function(x) x[[1]])
# # This seems to run for quite a bit.
# tmp <- map(tmp, function(x) {
#   x$Position <-  factor(x$Position, levels = x$Position)
#   return(x)
# })
# tmp <- map(tmp, function(x) {
#   x$numberInsertionsRel <- x$numberInsertions/sum(x$numberInsertions)
#   return(x)
# })
# do_rarefy <- function(df, depthToRarefy) {
#   print(sum(df$numberInsertionsRel))
#   stopifnot(near(sum(df$numberInsertionsRel), 1))
#   tmp <- sample(x = df$Position, size = depthToRarefy, replace = T, prob = df$numberInsertionsRel)
#   return(length(unique(tmp)))
# }
# 
# # Have finer curve at start...
# rarefy <- data.frame(subsample_depth = c(seq(from = 100 , to = 0.25 * sum(tmp[[1]]$numberInsertions), length.out = 50),
#                                          seq(from = 0.25 * sum(tmp[[1]]$numberInsertions), to = 0.5 * sum(tmp[[1]]$numberInsertions), length.out = 7)))
# rarefy$subsample_depth <- round(rarefy$subsample_depth)
# rarefy$subsample_depth_percentage_of_all_insertion_events <- (rarefy$subsample_depth/sum(tmp[[1]]$numberInsertions)) * 100
# rarefy$type <- "our_agar"
# tmpp <- data.frame(subsample_depth = c(seq(from = 100, to = 0.25 * sum(tmp[[1]]$numberInsertions), length.out = 50),
#                                        seq(from = 0.25 * sum(tmp[[2]]$numberInsertions), to = 0.5 * sum(tmp[[2]]$numberInsertions), length.out = 7)))
# tmpp$subsample_depth <- round(tmpp$subsample_depth)
# tmpp$subsample_depth_percentage_of_all_insertion_events <- (tmpp$subsample_depth/sum(tmp[[2]]$numberInsertions)) * 100
# tmpp$type <- "our_liquid"
# rarefy <- rbind(rarefy, tmpp)
# names(tmp) <- c("our_agar", "our_liquid")
# 
# rarefy <- rarefy %>% mutate(uniqueInsertionSites = map2_dbl(subsample_depth, type, function(x, y) do_rarefy(df = tmp[[y]], x)))
# 
# C <- ggplot(rarefy, aes(x = subsample_depth_percentage_of_all_insertion_events, y = uniqueInsertionSites, color = type, group = type)) + geom_line()
# 
# C1 <- ggplot(rarefy, aes(x = subsample_depth, y = uniqueInsertionSites, color = type, group = type)) + geom_line()
# 
# A <- A + theme_embl()
# B <- B + theme_embl()
# B1 <- B1 + theme_embl()
# C <- C + theme_embl()
# C1 <- C1 + theme_embl()
# 
# p <- A + (B/B1) + C + C1 + plot_layout(ncol = 4)
# ggsave(filename = "/home/nicolai/for_carlos_4/Fig2abc.pdf", plot = p, width = 15, height = 2.9)

