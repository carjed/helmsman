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

# due to how the deconstructSigs package parse VCFs, need to use test VCF without RSIDs for each variant
chr22.vcf <- "/mnt/norbert/data/1kg/chr22.noid.test.vcf"

# maf file for testing maftools package
lihc.maf <- "/mnt/norbert/data/tcga/15ce66c6-0211-4f03-bd41-568d0818a044/gsc_LIHC_pairs.aggregated.capture.tcga.uuid.automated.somatic.maf"

# unlike other packages, maftools reads reference genome from file on disk instead of using Bioconductor packages
ref.fasta <- "/mnt/norbert/data/ref/human_g1k_v37_min.fasta"

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
