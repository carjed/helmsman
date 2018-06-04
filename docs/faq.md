## FAQ

### What are Helmsman's hardware requirements?
Helmsman has no specific hardware dependencies. We recommend ~2GB of free disk space if storing the demo data.

### What type of sequencing data can Helmsman analyze?
Helmsman's core functionality, the NMF-based or PCA-based mutation signature analysis, can be used to rapidly evaluate mutation patterns present in any collection of sequencing data, and offers an extremely efficient and almost entirely automated (albeit slightly less flexible) alternative to the functions provided by the [SomaticSignatures](http://bioconductor.org/packages/release/bioc/html/SomaticSignatures.html) R package:

>Gehring JS, Fischer B, Lawrence M, Huber W. SomaticSignatures: inferring mutational signatures from single-nucleotide variants. Bioinformatics. 2015;31(22):3673-3675. doi:10.1093/bioinformatics/btv408.

Although we have only tested Helmsman with human data, it should be fully capable of analyzing NGS data from any other organism.

### Helmsman is not reading my input VCF
Helmsman uses the cyvcf2 wrapper for htslib libraries for parsing VCFs. If the input VCF is not formatted properly, htslib will often throw an error/warning or stop completely. If you encounter any issues with reading in a VCF file, first confirm that it is properly formatted. [vcf-validator](https://github.com/EBIvariation/vcf-validator) is a useful utility for this. You can also try indexing the VCF with `tabix -p input.vcf.gz` and running a simple bcftools command, such as `bcftools view input.vcf.gz`. Also check that your VCF contains singleton SNVs, individual genotypes are included, and that the FILTER column contains "PASS" values.

### Why am I getting errors about the reference genome?
To build the input matrix of subtype counts, Helmsman must annotate each site with the surrounding sequence context, using a user-specified reference genome file.

The libraries for parsing fasta files require that the chromosome records in your fasta reference file are formatted identically to the CHROM field of your VCF (e.g., if the fasta header is formatted as `>chrN` but the VCF CHROM field is just `N`). If you encounter any errors, try modifying the the fasta file to either strip or add "chr" to each record.
