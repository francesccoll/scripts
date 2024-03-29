# scripts
variety of different bioinformatic scripts

Script: functionality

* __annotate_repetitive_regions_on_genome.py__: Script to identify repetitive regions in a bacterial reference genome defined by BLASTing the reference genome against itself.
* __create_snippy_consensus.py__: Script to obtain a version of the reference genome with both substitution variants and missing calls initiated from Snippy output files
* __extract_gene_sequence.py__: Script to extract the DNA and protein sequence of a GFF annotated gene (tested on Bakta annotated GFF3 files). The script will output: gene information, gene DNA and protein sequences.
* __extract_gene_sequence_blast.py__: Script to extract the DNA sequence of a gene from an assembly by Blasting the gene sequence against the assembly. The script will output: gene information, gene DNA and protein sequences.
* __get_mapping_stats.py__: Script to obtain short-read mapping statistics from BAM and VCF files for QC purposes
* __get_msa_consensus.py__: Script to obtain a multiple sequence alignment and extract the consensus sequence from this. clustalw2, clustalo and kalign MSA tools supported. The script expects a FASTA file with multiple homologous sequences (DNA or protein), and will output the MSA and consensus sequence.
* __gff_to_table.py__: Python script to parse a GFF file to a CSV table format
* __prepare_vcf_file.py__: This script is used to re-format an input VCF file for downstream analysis. Specifically, it will check the input VCF file is a multi-sample VCF format; it will split multi-allelic sites; make sure GT genotypes are in haploid format, if not, convert diploid to haploid; add variant IDs as CHROM.POS.REF.ALT; and select subset of samples, if chosen.
* __vcf_to_table.py__: Script to convert a multi-sample VCF file into a matrix and saved as a CSV file (samples as rows, variants as columns). If the VCF is annotated, it will output a variant annotation table. If a BED file is specified, it will keep variants within the regions specified only (e.g. genes of interest) in the output table. A file with sample ids (one sample id per line) can be specified to keep a subset of samples in the output table.
