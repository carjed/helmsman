##############################################################################
# Test performance of mutation signature analysis packages
# - SomaticSignatures
# - deconstructSigs
# - signeR
# - maftools
##############################################################################

# pryr package for measuring memory usage of individual objects
# note that this is not the same as the memory usage reported in the manuscript!
require(pryr)
require(yaml)

# put file paths for the vcf, maf, and fasta files in helmsman/test/_config.yaml
keys1 <- yaml.load_file("helmsman/test/_config.yaml")

# due to how the deconstructSigs package parse VCFs, need to use test VCF without RSIDs for each variant
# this file is identical to data/chr22.16050075-16992642.1kg.phase3.vcf.gz, but with all values in the ID column set to "."
chr22.vcf <- keys1$vcf

# maf file for testing maftools package, downloaded from:
# https://portal.gdc.cancer.gov/legacy-archive/files/15ce66c6-0211-4f03-bd41-568d0818a044
lihc.maf <- keys1$maf

# unlike other packages, maftools reads reference genome from file on disk instead of using Bioconductor packages,
# so we need to specify where the fasta file is stored
ref.fasta <- keys1$fasta

#-----------------------------------------------------------------------------
# SomaticSignatures
#-----------------------------------------------------------------------------
require(SomaticSignatures)
require(BSgenome.Hsapiens.1000genomes.hs37d5)
start.time <- Sys.time()

mem_change(testvcf <- readVcfAsVRanges(chr22.vcf))
mem_change(genome(testvcf) <- "hs37d5")
mem_change(seqlevels(testvcf) <- "22")
mem_change(sca_motifs <- mutationContext(testvcf[grepl("1", testvcf$GT),], BSgenome.Hsapiens.1000genomes.hs37d5))
mem_change(sca_matrix <- motifMatrix(sca_motifs, group = "sampleNames", normalize = FALSE))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#-----------------------------------------------------------------------------
# deconstructSigs
#-----------------------------------------------------------------------------
require(deconstructSigs)
start.time <- Sys.time()


mem_change(vcf2 <- vcf.to.sigs.input(chr22.vcf))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#-----------------------------------------------------------------------------
# signeR
#-----------------------------------------------------------------------------
require(signeR)
require(VariantAnnotation)
require(BSgenome.Hsapiens.UCSC.hg19)

# had to rewrite function to handle phased genotypes
genCountMatrixFromVcf2 <- function(bsgenome, vcfobj) {

  vcfobj <- vcfobj[isSNV(vcfobj),]
  contexts <- getSeq(bsgenome, resize(granges(vcfobj), 3, fix="center"))
  alts <- alt(vcfobj)
  refs <- ref(vcfobj)
  
  gtsMat <- geno(vcfobj)$GT
  gtsMat <- structure(sapply(strsplit(gtsMat,"|",fixed=TRUE), unique), dim=dim(gtsMat))
  sample_names <- colnames(vcfobj)
  
  count_matrix <- matrix(0,
                         nrow=length(sample_names), ncol=length(change_triplet))
  
  rownames(count_matrix) <- sample_names
  colnames(count_matrix) <- change_triplet
  
  for(i in 1:nrow(gtsMat)) {
    rb <- refs[[i]]
    cc <- contexts[i]
    if(rb == DNAString("C") || rb == DNAString("T")) {
      cts <- sapply(alts[i], function(ab){paste0(rb,">",ab,":",cc)})
    } else {
      cts <- sapply(alts[i], function(ab){
        paste0(reverseComplement(rb),">",
               reverseComplement(ab),":",
               reverseComplement(cc))})
    }
    for(j in 1:length(sample_names)) {
      gi <- gtsMat[[i,j]]
      if(isGenoIndexRef(gi)) {next}
      if(isBadGeno(gi)){
        warning(paste0("Genotype of Line ",i," sample ",
                       sample_names[j]," is not supported, it will be skipped."))
        next
      }
      ai <- as.integer(getFirstGenoAltIndex(gi))
      ct <- cts[[ai]]
      cx <- match(ct, change_triplet)
      if(is.na(cx)) {
        warning(paste0("ct = ",ct))
        next
      }
      count_matrix[j,cx] = 1 + count_matrix[j,cx]
    }
  }
  return(count_matrix)
}


start.time <- Sys.time()

mem_change(vcfobj <- readVcf(chr22.vcf, "hg19"))
mem_change(seqlevelsStyle(vcfobj) <- "UCSC")
mem_change(seqlevels(vcfobj) <- "chr22")
bsgenome <- BSgenome.Hsapiens.UCSC.hg19
mem_change(mut <- genCountMatrixFromVcf2(bsgenome, vcfobj))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken

#-----------------------------------------------------------------------------
# MutationalPatterns
#-----------------------------------------------------------------------------
# require(MutationalPatterns)
# start.time <- Sys.time()
# 
# mem_change(vcfs <- read_vcfs_as_granges(c(chr22.vcf), genome = bsgenome))
# vcf <- readVcfAsVRanges(chr22.vcf)
# mem_change(mut_mat <- mut_matrix(vcf_list = c(vcf), ref_genome = bsgenome))
# 
# end.time <- Sys.time()
# time.taken <- end.time - start.time
# time.taken

#-----------------------------------------------------------------------------
# maftools
#-----------------------------------------------------------------------------
require(maftools)
start.time <- Sys.time()

mem_change(lihc <- read.maf(maf = lihc.maf, removeSilent = FALSE, useAll = TRUE))
mem_change(lihc.tnm <- trinucleotideMatrix(maf = lihc, ref_genome = ref.fasta, useSyn = TRUE))

end.time <- Sys.time()
time.taken <- end.time - start.time
time.taken
