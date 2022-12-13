#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Blast.Applications import NcbitblastnCommandline
from Bio.Blast import NCBIXML

# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------

# This script is used to extract the DNA sequence of a gene from an assembly by Blasting
# the gene sequence against the assembly


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to extract the DNA sequence of a gene from an assembly by Blasting " \
                  "the gene sequence against the assembly.\n" \
                  "The script will output: gene information, gene DNA and protein sequences."
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-f", "--input_fasta", action="store", dest="input_fasta",
        help="input genome in FASTA format (assembly or complete genome)",
        required=True, metavar="GENOME")
    group.add_argument(
        "-g", "--gene_sequence", action="store", dest="gene_sequence",
        help="DNA gene sequence in FASTA format to be extracted",
        required=True, metavar="GENE")
    group.add_argument(
        "-o", "--output_prefix", action="store", dest="output_prefix",
        help="File name prefix to name output files: gene information, gene DNA and protein sequences.",
        required=True, metavar="OUT")
    group.add_argument(
        "-i", "--fasta_id", action="store", dest="fasta_id",
        help="Id to be used in output FASTA files. Recommended to use sample/isolate id for MSA.",
        required=True, metavar="ID")

    group = parser.add_argument_group('Optional Blast arguments')
    group.add_argument(
        "-b", "--blast_tool", action="store", dest="blast_tool",
        help="BLAST tool to be used. Use \'blastn\' when blasting gene DNA sequence against genome, and \'tblastn\' "
             "when blasting a gene protein sequence against genome (default: blastn)",
        required=False, metavar="BTOOL", default="blastn"
    )
    group.add_argument(
        "-x", "--min_perid", action="store", dest="min_perid",
        help="Minimum percentage identity of blast alignment (default: 70)",
        required=False, metavar="PERID", type=float, default=70
    )
    group.add_argument(
        "-y", "--min_perlen", action="store", dest="min_perlen",
        help="Minimum percentage length of blast alignment (default: 80)",
        required=False, metavar="PERLEN", type=float, default=80
    )
    group.add_argument(
        "-e", "--evalue", action="store", dest="evalue",
        help="Minimum blast evalue (default: 0.05)",
        required=False, metavar="PERLEN", type=float, default=0.05
    )

    return parser.parse_args()


# ------------------------------------------------------------------------------------
# Main program
# ------------------------------------------------------------------------------------

def _main():
    # Configure logging
    logging.basicConfig(
        format='%(asctime)s %(levelname)s: %(message)s',
        level=logging.INFO
    )
    # Get arguments
    args = parse_arguments()

    # Making sure input files exist
    input_files = [args.input_fasta, args.gene_sequence]
    for input_file in input_files:
        if not os.path.isfile(input_file):
            logging.error(f'Input file {input_file} not found!')
            sys.exit(-1)

    # Making sure optional Blast parameters are numeric
    options = ["min_perid", "min_perlen", "evalue"]
    values = [args.min_perid, args.min_perlen, args.evalue]
    for option, value in zip(options, values):
        if value is not None:
            if not (type(value) == int or type(value) == float):
                logging.error(f'The chosen value {str(value)} for parameter --{str(option)} is not numeric!')
                sys.exit(-1)

    # Making sure chosen blast tools supported
    blast_tool = args.blast_tool
    blast_tools = ["blastn", "tblastn"]
    if blast_tool not in blast_tools:
        logging.error(f'The chosen --blast_tool {str(blast_tool)} not supported!')
        sys.exit(-1)

    # Reading query and subject sequences and printing information
    logging.info(f"Opening input genome {args.input_fasta}")
    input_records = SeqIO.parse(args.input_fasta, "fasta")
    record_num = 0
    record_nt_length = 0
    for record in input_records:
        record_num += 1
        record_nt_length += len(record.seq)
    logging.info(f"Input genome has {str(record_num)} contigs and {str(record_nt_length)} nucleotide length")

    logging.info(f"Opening input gene {args.gene_sequence}")
    input_records = SeqIO.parse(args.gene_sequence, "fasta")
    record_num = 0
    gene_nt_length = 0
    for record in input_records:
        record_num += 1
        gene_nt_length += len(record.seq)
    # Making sure gene file only contains one sequence
    if record_num > 1:
        logging.error(f'Gene file {args.gene_sequence} contains more than one sequence!')
        sys.exit(-1)
    logging.info(f"Input gene has {str(gene_nt_length)} nucleotide length")

    # Blasting gene against genome
    # NOTE: blastn results are saved as "5 = XML Blast output" to print out alignments
    logging.info(f"Blasting gene {args.gene_sequence} against genome {args.input_fasta}")
    blast_results_xml = args.output_prefix + '.blast_output.xml'
    blastn_cline = ''
    if blast_tool == 'blastn':
        blastn_cline = NcbiblastnCommandline(
            query=args.gene_sequence,
            subject=args.input_fasta,
            out=blast_results_xml,
            outfmt=5,
            task='blastn')
    if blast_tool == 'tblastn':
        blastn_cline = NcbitblastnCommandline(
            query=args.gene_sequence,
            subject=args.input_fasta,
            out=blast_results_xml,
            outfmt=5,
            task='tblastn')
    print(blastn_cline)
    blastn_cline()
    # NOTE: blastn results are saved as a custom "6 = tabular" to keep custom alignment results
    blast_results_tab = args.output_prefix + '.blast_output.tab'
    blastn_cline = ''
    if blast_tool == 'blastn':
        blastn_cline = NcbiblastnCommandline(
            query=args.gene_sequence,
            subject=args.input_fasta,
            out=blast_results_tab,
            outfmt="6 qseqid sseqid sstart send length evalue pident sstrand slen",
            task='blastn')
    if blast_tool == 'tblastn':
        blastn_cline = NcbitblastnCommandline(
            query=args.gene_sequence,
            subject=args.input_fasta,
            out=blast_results_tab,
            outfmt="6 qseqid sseqid sstart send length evalue pident sstrand slen",
            task='tblastn')
    print(blastn_cline)
    blastn_cline()

    # Making sure Blast result output files exist
    if not os.path.isfile(blast_results_xml):
        logging.error(f'Blast result output file {blast_results_xml} not found. Something went wrong!')
        sys.exit(-1)
    if not os.path.isfile(blast_results_tab):
        logging.error(f'Blast result output file {blast_results_tab} not found. Something went wrong!')
        sys.exit(-1)

    # Printing Blast alignments
    logging.info(f"Printing Blast alignments from ({blast_results_xml})")
    blast_record = NCBIXML.read(open(blast_results_xml, "r"))
    for alignment in blast_record.alignments:
        for hsp in alignment.hsps:
            if hsp.expect < 0.01:
                print('\n')
                print('****Alignment****')
                print('sequence:', alignment.title)
                print('length:', alignment.length)
                print('score:', hsp.score)
                print('align_length:', hsp.align_length)
                print('sbjct start:', hsp.sbjct_start)
                print('sbjct end:', hsp.sbjct_end)
                print('sbjct strand:', hsp.strand)
                print('identities:', hsp.identities)
                print('gaps:', hsp.gaps)
                print(hsp.query)
                print(hsp.match)
                print(hsp.sbjct)

    # Parsing Blast alignments to save gene coordinates
    print('\n')
    occurrences = 0  # number of times gene (or gene fragment) found
    out_lines = list()  # list to save gene coordinates in genome, including multiple instances of gene
    logging.info(f"Parsing Blast alignments from ({blast_results_tab}) to save gene coordinates")
    logging.info(f"Minimum percentage identity of blast alignment to keep: {str(args.min_perid)}")
    logging.info(f"Minimum percentage length of blast alignment to keep: {str(args.min_perlen)}")

    for line in open(blast_results_tab, "r"):
        (qseqid, sseqid, sstart, send, length, evalue, pident, sstrand, slen) = line.strip().split('\t')
        per_length = int(length) / int(gene_nt_length) * 100
        if float(per_length) >= float(args.min_perlen) and float(pident) >= float(args.min_perid):
            print('\t --> Blast alignment used to extract gene coordinates')
            print('\t\t --> ' + str(per_length) + ' percentage of query gene covered by alignment')
            print('\t\t --> ' + str(pident) + ' percentage of identical matches')
            print('\t' + line)
            occurrences += 1
            gene_name = qseqid
            gene_contig = sseqid
            contig_length = slen
            gene_from = sstart
            gene_to = send
            align_length = length
            gene_strand = "+" if sstrand == 'plus' else "-"
            if gene_strand == "-":
                tmp = gene_from
                gene_from = gene_to
                gene_to = tmp
            out_items = ["gene", gene_name, gene_contig, contig_length, gene_from, gene_to, gene_strand,
                         str(align_length)]
            print('\t'.join(out_items))
            out_lines.append('\t'.join(out_items))

    # Saving gene information, and DNA and protein sequence of extracted gene sequence
    print('\n')
    output_info = args.output_prefix + ".gene_info.txt"
    logging.info(f"Saving extracted gene into {output_info}")
    output = open(output_info, 'w')
    out_cols = ["sequence", "gene_name", "gene_contig", "contig_length", "gene_from", "gene_to", "gene_strand",
                "gene_length"]
    output.write('\t'.join(out_cols) + '\n')

    fasta_id = args.fasta_id
    output_seq1 = args.output_prefix + ".dna.fa"
    output_seq2 = args.output_prefix + ".protein.fa"
    if occurrences > 0:
        for out_line in out_lines:
            output.write(out_line + '\n')
            [sequence, gene_name, gene_contig, contig_length, gene_from, gene_to, gene_strand, *_] = out_line.split('\t')
            # Opening genome sequence to extract and save gene sequence
            logging.info(f"Opening genome sequence {args.input_fasta} to extract and save gene sequence...")
            input_records = SeqIO.parse(args.input_fasta, "fasta")
            for record in input_records:
                record_id = record.id
                print(record_id)
                if record_id == gene_contig:
                    record_seq = record.seq
                    if gene_strand == "+":
                        gene_seq = record_seq[int(gene_from)-1:int(gene_to)]
                    if gene_strand == "-":
                        gene_seq = record_seq[int(gene_from)-1:int(gene_to)]
                        gene_seq = gene_seq.reverse_complement()
                    print(gene_seq)
                    logging.info(f"Saving extracted gene sequence into {output_seq1}")
                    gene_seq_record = SeqRecord(gene_seq, id=fasta_id)
                    SeqIO.write(gene_seq_record, output_seq1, "fasta")
                    logging.info(f"Saving extracted protein sequence into {output_seq2}")
                    protein_seq = gene_seq.translate()
                    protein_seq_record = SeqRecord(protein_seq, id=fasta_id)
                    SeqIO.write(protein_seq_record, output_seq2, "fasta")
    else:
        # Saving gene information output as "not_found"
        gene_name = ""
        # Extracting gene name
        input_records = SeqIO.parse(args.gene_sequence, "fasta")
        for record in input_records:
            gene_name = record.id
        out_items = ["gene", gene_name, "not_found", "not_found", "not_found", "not_found", "not_found", "not_found"]
        out_line = '\t'.join(out_items)
        output.write(out_line + '\n')

    output.close()


if __name__ == "__main__":
    _main()
