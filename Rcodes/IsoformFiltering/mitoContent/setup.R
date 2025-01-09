
library(ezRun)
library(tidyverse)
library(parallel)


# samples = ezRead.table("/home/zajac/Genomics/p28443_singlecell_pacbiovsillumina/2024-01-17-sampleInfo.tsv")
# names <- rownames(samples) %>% recode("G18_626_Normal" = "Normal",
#                "G18_626_ccRCC" = "ccRCC_2",
#                "G16_315_ccRCC" = "ccRCC_3",
#                "G18_109_ccRCC" = "ccRCC_4",
#                "G16_667_ccRCC" = "ccRCC_5")
# 
# samples <- ezFrame(Name=names, samples)
# ezWrite.table(samples, file="./2024-05-22-withPDAC-sampleInfo.tsv", head="SampleId")
# #samples <- ezRead.table(file="./2024-05-22-withPDAC-sampleInfo.tsv")
# sampleInfo <- samples[!grepl("PDAC", rownames(samples)), ]
# sampleInfo$nIlluminaReads <- sapply(file.path("/srv/gstore/projects", sampleInfo$CellRangerNoIntron, "metrics_summary.csv"),
#                                     function(summaryFile){
#                                       xx <- data.table::fread(summaryFile)
#                                       xx$`Number of Reads` %>% str_replace_all(",", "") %>% as.integer()
#                                     })
# sampleInfo$percentIlluminaMapped <- sapply(file.path("/srv/gstore/projects", sampleInfo$CellRangerNoIntron, "metrics_summary.csv"),
#                                            function(summaryFile){
#                                              xx <- data.table::fread(summaryFile)
#                                              xx$`Reads Mapped to Genome` %>% str_replace_all("%", "") %>% as.double()
#                                            })
# sampleInfo$nPbSegmentReads <- c(54962298, 58333415, 34180213, 29404003, 35635073)
# sampleInfo$percentPbMapped <- c(99.55,
#                                 99.64,
#                                 99.39,
#                                 99.61,
#                                 99.69)
# ezWrite.table(sampleInfo, file="./2024-05-22-sampleInfo.tsv", head="SampleId")
# 
sampleInfo <- ezRead.table(file="./2024-05-22-sampleInfo.tsv")#[1:2, ]


pbDirs <- setNames(paste0(sampleInfo$PigeonFiltered, "genes_seurat"),
                   sampleInfo$Name)
pbDirs_unfiltered <- setNames(paste0(sampleInfo$PigeonUnfiltered, "genes_seurat"),
                   sampleInfo$Name)
pbDirs_filtered_IntraPriming <- pbDirs %>% str_replace("getMTinfo", "getMTinfo_filteredIntraPriming") %>% setNames(sampleInfo$Name)

illDirs <- setNames(paste0("/srv/gstore/projects/", sampleInfo$CellRangerNoIntron, "/filtered_feature_bc_matrix"),
                    sampleInfo$Name)

isoSeqBam <- setNames(sampleInfo$IsoSeqBam,
                      sampleInfo$Name)

sampleInfo$isoSeqFiltReasons <- dirname(sampleInfo$IsoSeqBam) %>% file.path("../../call-pigeon/execution/scisoseq_classification.filtered_lite_reasons.txt") %>%
  normalizePath()


geneAnno <- ezRead.table("/srv/GT/reference/Homo_sapiens/GENCODE/GRCh38.p13/Annotation/Release_39-2021-12-09/Genes/SMRTLink11.1_annotation/genes_annotation_byGene.txt")
geneAnno$seurat_name <- gsub("_", "-", geneAnno$gene_name)
geneAnno$myType <- geneAnno$type
geneAnno$myType[geneAnno$seqid == "chrM"] <- "mito"
geneAnno$myType[grepl("^RPL|RPS", geneAnno$gene_name)] <- "riboprot"
#table(geneAnno$type, geneAnno$myType)
seuratPal <- scales::hue_pal()(5)
names(seuratPal) <- sampleInfo$Name


# myColorsv1 <-  c("Both"="#440154FF", "Common"="#440154FF", 
#                "Illumina specific"="darkgoldenrod", "PacBio specific"="#ffcf20FF",
#                "Illumina only"="darkgoldenrod", "PacBio only"="#ffcf20FF",
#                "Illumina"="darkgoldenrod", "PacBio"="#ffcf20FF",
#                "pbDetected"="darkgoldenrod",
#                "pbAbsent"="gray99",
#                "illDetected"="darkgoldenrod",
#                "illAbsent"="gray99",
#                seuratPal)

## the RcolorBrewer palette 'Dark2'

# the colors are mainly from: RColorBrewer::brewer.pal(8, "Dark2")

sampleColors <- c(Normal="#7570B3", ccRCC_2="#1B9E77", ccRCC_3="#D95F02",
                  ccRCC_4="#E7298A",
                  ccRCC_5="#66A61E")
#  "#666666"
# pacbio was: "#ffcf20FF"

myColors <-  c("Both"="#440154FF", "Common"="#440154FF",
               "Illumina specific"="#A6761D", "PacBio specific"="#E6AB02",
               "Illumina only"="#A6761D", "PacBio only"="#E6AB02",
               "Illumina"="#A6761D", "PacBio"="#E6AB02",
               "pbDetected"="#A6761D",
               "pbAbsent"="gray99",
               "illDetected"="#A6761D",
               "illAbsent"="gray99",
               sampleColors,
               "polyA"="#D95F02", "TSO"="#E7298A", "TSO+polyA"="#7570B3",
               "unmapped"="#1B9E77", "multi"="#E7298A", "unique"="#7570B3",
               "<= 500"="#E7298A")






#Read in cell annotations
# cellAnnot <- read.delim("/srv/GT/analysis/zajacn/p28443/ccRCC_CellType.tsv")
# cellAnnot$CellBarcode <- reverseComplement(DNAStringSet(str_remove(cellAnnot$CellBarcode, "-1"))) %>% as.character()
# cellAnnot = rbind(cellAnnot, data.frame("CellType" = "non-ccRCC", "CellBarcode" = "CTATCCGGTCACCTTC", "Sample" = "ccRCC_4"))
# cellAnnot = rbind(cellAnnot, data.frame("CellType" = rep("non-ccRCC", ncol(pbData$Pacbio$ccRCC_3)), 
#                                         "CellBarcode" = colnames(pbData$Pacbio$ccRCC_3), 
#                                         "Sample" = rep("ccRCC_3", ncol(pbData$Pacbio$ccRCC_3))))
# cellAnnot = rbind(cellAnnot, data.frame("CellType" = rep("non-ccRCC", ncol(pbData$Pacbio$Normal)), 
#                                         "CellBarcode" = colnames(pbData$Pacbio$Normal), 
#                                         "Sample" = rep("Normal", ncol(pbData$Pacbio$Normal))))
# ezWrite.table(cellAnnot, "2024-05-22-ccRCC_CellType.tsv", row.names = FALSE)
cellAnnot <- ezRead.table("2024-05-22-ccRCC_CellType.tsv", row.names = NULL)


compileMatchedReads_20240603 <- function(sampleInfo, cbStart="AA", fastOverlapOnly=FALSE, nCores=4,
                                         scratchDir="/scratch/hubert/p28443_20240603"){
  
  if (!file.exists(scratchDir)){
    dir.create(scratchDir)
  }
  i <- 1
  sampleInfo$pbGaFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-pbGa.qsd")
  sampleInfo$illGaFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-illGa.qsd")
  sampleInfo$illUnmappedFile <- paste0(scratchDir, "/", sampleInfo$Name, "-", cbStart, "-illUnmapped.qsd")
  for (i in seq_along(sampleInfo$Name)){
    smi <- sampleInfo[i, ]
    message(smi$Name)
    if (!file.exists(smi$pbGaFile)){
      tmpSam <- paste0(scratchDir, "/", smi$Name, "-", cbStart, ".sam")
      illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
      illBam <- sub(".sam$", ".bam", tmpSam)
      if (!file.exists(illBam)){
        ezSystem(paste0("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view -H ", 
                        illDir, "/possorted_genome_bam.bam > ",  tmpSam))
        ezSystem(paste0("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view ", 
                        illDir, "/possorted_genome_bam.bam| grep CB:Z:", cbStart, " >>",  tmpSam))
        ezSystem(paste("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools view -b", tmpSam, ">", illBam))
        ezSystem(paste("/usr/local/ngseq/packages/Tools/samtools/1.11/bin/samtools index", illBam))
        file.remove(tmpSam)
      }
      ## read the illumina alignments
      sbp <- ScanBamParam(what =  c("qname", "seq"), tag = c("GN", "xf", "CB", "UB", "UR", "ts", "pa", "RE", "NH"))
      illGa <- readGAlignments(illBam, param = sbp)
      mcols(illGa)$CB <- sub("-1", "", mcols(illGa)$CB)
      ## replace the undefined corrected barcodes with the raw barcode
      toReplace <- is.na(mcols(illGa)$UB)
      mcols(illGa)$UB[toReplace] <- mcols(illGa)$UR[toReplace]
      mcols(illGa)$UR <- NULL
      
      sbp <- ScanBamParam(what = c("qname", "seq"), tag = c("xf", "CB", "UB", "ts", "pa"), 
                          scanBamFlag(isUnmappedQuery = TRUE))
      illBamList <- scanBam(illBam, param = sbp)[[1]]
      illBamList$tag$CB <- sub("-1", "", illBamList$tag$CB)
      illUnmapped <- cbind(qname=illBamList$qname, data.frame(illBamList$tag))
      illUnmapped$gc <- letterFrequency(illBamList$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(illBamList$seq)
      illUnmapped$qlength <- width(illBamList$seq)
      
      #Load the Pacbio data and keep only same cells in all datasets
      sbp <- ScanBamParam(what = c("qname", "seq"), tag = c( "CB", "XM"), tagFilter=list(rc=1))
      pbGa <- readGAlignments(smi$IsoSeqBam, param = sbp)
      ## reverse complement the barcode to match illumina
      mcols(pbGa)$CB <- mcols(pbGa)$CB %>% DNAStringSet() %>% reverseComplement() %>% as.character()
      
      ## intersect the alignments with the illumina alignments by cell barcode; pb alignments on contain pb real cells
      pbGa <- pbGa[mcols(pbGa)$CB %in% mcols(illGa)$CB]
      illGa <- illGa[mcols(illGa)$CB %in% mcols(pbGa)$CB]
      illUnmapped <- illUnmapped[illUnmapped$CB %in% mcols(pbGa)$CB, ]
      mcols(pbGa)$tagId <- paste(mcols(pbGa)$CB,
                                 as.character(reverseComplement(DNAStringSet(mcols(pbGa)$XM))))
      mcols(illGa)$tagId <- paste(mcols(illGa)$CB, mcols(illGa)$UB)
      illUnmapped$tagId <- paste(illUnmapped$CB, illUnmapped$UB)
      
      pbAnno <- data.table::fread(file.path(smi$PigeonUnfiltered, "scisoseq.annotated.info.csv"), sep="\t")
      pbAnno <- pbAnno[pbAnno$BCrev %in% mcols(illGa)$CB, ]
      
      pbReasons <- data.table::fread(smi$isoSeqFiltReasons, sep=",", skip="filtered_isoform,filter")
      pbAnno$quality <- pbReasons[match(pbAnno$pbid, pbReasons$filtered_isoform), "filter"] %>% unlist() %>% replace_na("Pass")
      stopifnot(pbAnno$id %in% mcols(pbGa)$qname)
      idx <- match(mcols(pbGa)$qname, pbAnno$id)
      mcols(pbGa)$category <- pbAnno$category[idx] %>% replace_na("unannotated")
      mcols(pbGa)$gene <- pbAnno$gene[idx]
      mcols(pbGa)$transcriptLength <- pbAnno$length[idx]
      mcols(pbGa)$quality <- pbAnno$quality[idx] %>% replace_na("Unknown")
      
      commonTagIds <- intersect(mcols(pbGa)$tagId, mcols(illGa)$tagId)
      mcols(illGa)$pbStatus <- ifelse(mcols(illGa)$tagId %in% commonTagIds, "pbDetected", "pbAbsent")
      mcols(pbGa)$illStatus <- ifelse(mcols(pbGa)$tagId %in% commonTagIds, "illDetected", "illAbsent")
      mcols(pbGa)$illUnmapped <- ifelse(mcols(pbGa)$tagId %in% illUnmapped$tagId, "illUnmapped", NA)
      mcols(pbGa)$gc <- letterFrequency(mcols(pbGa)$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(mcols(pbGa)$seq)
      mcols(pbGa)$qlength <- width(mcols(pbGa)$seq)
      mcols(illGa)$gc <- letterFrequency(mcols(illGa)$seq, letters="GC", as.prob = FALSE)[ ,"G|C"] /width(mcols(illGa)$seq)
      mcols(illGa)$qlength <- width(mcols(illGa)$seq)
      
      
      mcols(pbGa)$seq <- NULL
      mcols(illGa)$seq <- NULL
      ## subset for faster analysis:
      ### commonTagIds <- head(commonTagIds, 1000)
      if (!fastOverlapOnly){
        ## turn int GRangesList
        job <- ezJobStart("map")
        pbRgList <- as(pbGa[mcols(pbGa)$tagId %in% commonTagIds], "GRangesList")
        ezWriteElapsed(job, status="Grangeslist")
        # ## regroup by tagid
        pbTagId <- mcols(pbRgList)$tagId
        mcols(pbRgList) <- NULL
        pbRgs <- unlist(pbRgList)
        pbRgs$tagId <- rep(pbTagId, elementNROWS(pbRgList))
        pbRgList <- split(pbRgs, pbRgs$tagId) ## this is a GRangesList
        ## group pbAlignments by tagId and discard the mcols; this is a normal array/list -- insanely slow
        # pbRgList <- as(pbGa, "GRangesList")
        # mcols(pbRgList) <- NULL
        # pbRgList <- tapply(pbRgList[1:10000], mcols(pbGa)$tagId[1:10000], unlist, simplify = FALSE)
        ezWriteElapsed(job, status="Grangeslist regrouped")
        useIlluGa <- mcols(illGa)$tagId %in% commonTagIds
        myIllTagId <- mcols(illGa)$tagId[useIlluGa]
        illCommRgList <-  as(illGa[useIlluGa], "GRangesList") # [mcols(illGa)$tagId %in% names(pbRgList)]
        mcols(illCommRgList) <- NULL
        ezWriteElapsed(job, status=paste("illu Grangeslist subsetted to", length(illCommRgList)))
        
        ### computing the overlap is the slow part
        ### using mapply would be slower than mclapply
        # hasOverlap <- mapply(function(x,y){any(overlapsAny(x,y, minoverlap=5, type="any"))},
        #               illCommRgList, pbRgList[myIllTagId])
        ### also simplifying the illumina ranges with: endoapply(illCommRgList, GenomicRanges::reduce)
        ### is not faster
        hasOverlap <- mclapply(1:length(illCommRgList), function(i){ 
          if (i %% 100000 == 0) ezWriteElapsed(job, status=paste("mapped", i))
          any(overlapsAny(illCommRgList[[i]], pbRgList[[myIllTagId[i]]], minoverlap=5, type="any"))
        }, mc.cores=nCores)
        hasOverlap <- unlist(hasOverlap)
        ezWriteElapsed(job, status="mapping done")
        mcols(illGa)$overlapsPbAlign <- NA
        mcols(illGa)$overlapsPbAlign[useIlluGa] <- ifelse(hasOverlap, "hasOverlapPb", "noOverlapPb")
        
        gc()
        ezWriteElapsed(job, status="cleaned up")
      }
      # fast first version which ignores the cigar information; and maps to the first Pb alignment of the tagId
      job <- ezJobStart("fast map")
      useIlluGa <- mcols(illGa)$tagId %in% commonTagIds
      illCommRgs <- granges(illGa[useIlluGa], use.mcols=FALSE)
      pbCommRgs <- granges(pbGa, use.mcols=FALSE)[match(mcols(illGa)$tagId[useIlluGa], mcols(pbGa)$tagId)]
      ovlStatus <- poverlaps(illCommRgs, pbCommRgs,
                             type="any", minoverlap=10) %>% as.vector()
      mcols(illGa)$overlapsPbAlignRange <- NA
      mcols(illGa)$overlapsPbAlignRange[useIlluGa] <- ifelse(as.vector(ovlStatus), "hasOverlapPb", "noOverlapPb")
      ezWriteElapsed(job, status="done")
      # as.vector(ovlStatus[match(mcols(illCommRgs)$tagId, names(ovlStatus)])
      # illCommRgs <- granges(illGa[match(commonTagIds, mcols(illGa)$tagId)], use.mcols=TRUE)
      # pbCommRgs <- granges(pbGa[match(commonTagIds, mcols(pbGa)$tagId)], use.mcols=TRUE)
      # # 
      # ovlStatus <- poverlaps(illCommRgs, pbCommRgs,
      #                        type="any", minoverlap=10) %>% as.vector()
      # # ### summarize duplicates by tagId
      # ovlStatus <- tapply(ovlStatus, mcols(illCommRgs)$qname, any)
      # mcols(illGa)$overlapsPbAlignRange <- as.vector(ovlStatus[mcols(illGa)$qname])

      isMultiMapper <- mcols(pbGa)$qname %in% mcols(pbGa)$qname[duplicated(mcols(pbGa)$qname)]
      mcols(pbGa)$multiMapping <-  ifelse(isMultiMapper, "Multi mapping", "Uniquely mapping")
      mcols(illGa)$umiCount = ifelse(mcols(illGa)$xf == "25", "Counted", "NotCounted")
      mcols(illGa)$Sample <- smi$Name
      mcols(pbGa)$Sample <- smi$Name
      qs::qsave(illGa, file=smi$illGaFile)
      qs::qsave(pbGa, file=smi$pbGaFile)
      qs::qsave(illUnmapped, file=smi$illUnmappedFile)
      remove(illGa, pbGa, illCommRgs, pbCommRgs)
      gc()
    }
  }
  return(sampleInfo)
}


compileCellStats_20240603 <- function(sampleInfo, 
                                         scratchDir="/scratch/hubert/p28443_20240603"){
  
  if (!file.exists(scratchDir)){
    stop()
  }
  i <- 1
  sampleInfo$cellStatsFile <- paste0(scratchDir, "/", sampleInfo$Name, "-cellStats.qsd")
  for (i in seq_along(sampleInfo$Name)){
    smi <- sampleInfo[i, ]
    message(smi$Name)
    if (!file.exists(smi$cellStatsFile)){
      illDir <- paste0("/srv/gstore/projects/", smi$CellRangerNoIntron)
      sbp <- ScanBamParam(what =  c("qname", "rname", "qwidth"), tag = c("CB"))
      pbBamList <- scanBam(smi$IsoSeqBam, param = sbp)[[1]]
      pbAnno <- data.table::fread(file.path(smi$PigeonUnfiltered, "scisoseq.annotated.info.csv"), sep="\t")
      #pbAnno <- pbAnno[pbAnno$BCrev %in% rownames(cellMito), ]
      pbBamList$category <- pbAnno$category[match(pbBamList$qname, pbAnno$id)]
      isMito <- pbBamList$rname %in% "chrM"
      pbPercentMito <- tapply(isMito, pbBamList$tag$CB, mean) * 100
      #pbPercentMito <- pbPercentMito[!is.na(names(pbPercentMito))]
      pbBc <- as.character(reverseComplement(DNAStringSet(names(pbPercentMito))))
      percentShort <- tapply(pbBamList$qwidth[!isMito] < 500, pbBamList$tag$CB[!isMito], mean) * 100
      catTable <- table(pbBamList$tag$CB[!isMito], pbBamList$category[!isMito], useNA="ifany")
      class(catTable) <- "matrix"
      colnames(catTable) <- colnames(catTable) %>% replace_na("undefined_category")
      

      illBamList <- scanBam(paste0(illDir, "/possorted_genome_bam.bam"), param = sbp)[[1]]
      use <- illBamList$tag$CB %in% paste0(pbBc, "-1")
      xt <- table(illBamList$rname[use], illBamList$tag$CB[use], useNA="always")
      illPercentMito <- xt["chrM", ]/colSums(xt) * 100
      names(illPercentMito) <- sub("-1", "", names(illPercentMito))
      cellStats <- ezFrame(pbPercentMito=as.vector(pbPercentMito), 
                          illPercentMito=as.vector(illPercentMito[pbBc]),
                          percentShort=as.vector(percentShort[names(pbPercentMito)]),
                          catTable[names(pbPercentMito), ],
                          row.names=pbBc
                          )
      qs::qsave(cellStats, smi$cellStatsFile)
    }
  }
  return(sampleInfo)
}



# pbAnnot = sapply(lapply(str_split(pbDirs, "/"), .subset, 1:8), paste, collapse = "/")
# pbAnnot_unfiltered = sapply(lapply(str_split(pbDirs_unfiltered, "/"), .subset, 1:8), paste, collapse = "/")
# pbAnnot_filteredIP = sapply(lapply(str_split(pbDirs_filtered_IntraPriming, "/"), .subset, 1:8), paste, collapse = "/")

# pbAnnot = lapply(pbAnnot, function(x){
#   df = read.csv(paste0(x,"/scisoseq.annotated.info.csv"), sep = "\t")
# })
# 
# pbAnnot_unfiltered = lapply(pbAnnot_unfiltered, function(x){
#   df = read.csv(paste0(x,"/scisoseq.annotated.info.csv"), sep = "\t")
# })
# 
# pbAnnot_filteredIP = lapply(pbAnnot_filteredIP, function(x){
#   df = read.csv(paste0(x,"/scisoseq.annotated.info.csv"), sep = "\t")
# })
# 
# rt = read.delim("../RTswtiching_isoforms.txt", header = F) #preselected from reasons files
# colnames(rt) = c("pbid", "Sample")
# rt$SampleID = case_when(rt$Sample == "Normal" ~ 1, rt$Sample == "ccRCC_2" ~ 2, rt$Sample == "ccRCC_3" ~ 3, rt$Sample == "ccRCC_4" ~ 4, rt$Sample == "ccRCC_5" ~ 5)
