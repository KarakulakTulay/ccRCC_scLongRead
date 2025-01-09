## Natalia Zajac

##Filter isoforms
sample_data = read_delim("/home/zajac/Genomics/p28443_singlecell_pacbiovsillumina/SIB-poster/2024-05-22-sampleInfo.tsv")
sample_data$isoSeqFiltClassification = paste0(str_remove(sample_data$isoSeqFiltReasons, "scisoseq_classification.filtered_lite_reasons.txt"), "scisoseq_classification.filtered_lite_classification.txt")
sample_data$PigeonIsoformFiltering = paste0("/srv/GT/analysis/zajacn/p28443/o31598_customFilteredIsoforms/", str_remove(sapply(str_split(sample_data$PigeonUnfiltered, "/"), .subset, 8), "_getMTinfo_unfiltered"), "_072024IsoformFiltering")

for (i in sample_data$Name){
  smi = sample_data[sample_data$Name == i,]
  df = read.delim(smi$isoSeqFiltClassification)
  fsm = df[df$structural_category == "full-splice_match",]
  others = df[df$structural_category != "full-splice_match",]
  others = others[which(abs(others$dist_to_cage_peak) <= 50 & 
                          others$within_cage_peak == TRUE & 
                          abs(others$diff_to_gene_TSS) <= 50 & 
                          abs(others$diff_to_gene_TTS) <= 50),]
  df = rbind(fsm,others)
  write_delim(df, paste0(smi$PigeonIsoformFiltering, "/scisoseq_classification.filtered_lite_custom_classification.txt"), delim  = "\t")
}

sample_data$isoseqdedupFasta = paste0(sapply(str_split(sample_data$isoSeqFiltClassification, "/call-pigeon"), .subset, 1), "/call-isoseq_dedup/execution/scisoseq.5p--3p.tagged.refined.corrected.sorted.dedup.fasta")
sample_data$isoseqCollapseTranscriptGroup = list.files(paste0(sapply(str_split(sample_data$isoSeqFiltClassification, "execution/"), .subset, 1), "inputs"), full.names = TRUE, recursive=TRUE, pattern = "scisoseq.mapped_transcripts.collapse.group.txt")[c(1,5,3,2,4)]

##Rerun pigeon make seurat
for (i in sample_data$Name){
  smi = sample_data[sample_data$Name == i,]
  cmd = paste("/srv/GT/software/SMRTtools/SMRT_Link_v10/smrtcmds/bin/pigeon make-seurat",  
              "--log-level INFO --log-file pigeon-make-seurat.log --num-threads 7", 
              "--dedup", smi$isoseqdedupFasta, 
              "--group", smi$isoseqCollapseTranscriptGroup, "--out-dir", smi$PigeonIsoformFiltering,
              "--out-prefix scisoseq", 
              paste0(smi$PigeonIsoformFiltering, "/scisoseq_classification.filtered_lite_custom_classification.txt"))
  ezSystem(cmd)
}

