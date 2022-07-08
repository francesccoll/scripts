#!/usr/bin/env python3

import string
import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to obtain a version of the reference genome with both substitution variants and " \
                  "missing calls initiated from Snippy output files"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-r", "--reference", action="store", dest="reference",
        help="FASTA file of reference genome used by Snippy", required=True, metavar="REFERENCE")
    group.add_argument(
        "-c", "--consensus_subs", action="store", dest="consensus_subs",
        help="A version of the reference genome with only substitution variants instantiated "
             "(i.e. Snippy's .consensus.subs.fa output file)",
        required=True, metavar="CONS")
    group.add_argument(
        "-a", "--aligned", action="store", dest="aligned",
        help="A version of the reference genome with zero coverage (-), low coverage (N), and heterozygous or "
             "poor quality genotypes (n) initiated (i.e. Snippy's .aligned.fa output file)",
        required=True, metavar="ALIGNED")
    group.add_argument(
        "-o", "--output", action="store", dest="output",
        help="A version of the reference genome with both substitution variants and missing calls initiated",
        required=True, metavar="OUTPUT")

    group = parser.add_argument_group('optional arguments')
    group.add_argument(
        "-i", "--sample_id", action="store", dest="sample_id",
        help="sample id to be used in header of output fasta file", required=False, metavar="SAMPLE_ID")
    group.add_argument(
        "-s", "--mapping_stats", action="store", dest="mapping_stats",
        help="Output file name to save basic mapping stats: number of consensus and missing SNP calls",
        required=False, metavar="STATS")
    group.add_argument(
        '--version', action='version', version='%(prog)s 1.0')

    return parser.parse_args()


def check_file_exist(file, file_tag):
    if not os.path.isfile(file):
        logging.error(f'{file_tag} {file} not found!')
        sys.exit(-1)
    else:
        logging.info(f'{file_tag} {file} found!')


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

    # Making sure required input files exist
    check_file_exist(args.reference, "Input reference file")
    check_file_exist(args.consensus_subs, "Input Snippy's .consensus.subs.fa file")
    check_file_exist(args.aligned, "Input Snippy's .aligned.fa file")

    # Making sure all three FASTA files have the same length, and save sequences into lists
    # NOTE: by default, the fasta id of .consensus.subs.fa will be used for the output fasta file, unless otherwise
    # specified with the option --output_fasta_id
    fasta_id = ""
    ref_length = 0
    ref_seq = list()
    con_length = 0
    con_seq = list()
    aln_length = 0
    aln_seq = list()
    for record in SeqIO.parse(args.reference, "fasta"):
        ref_length = len(record.seq)
        ref_seq = list(str(record.seq).upper())
    for record in SeqIO.parse(args.consensus_subs, "fasta"):
        con_length = len(record.seq)
        con_seq = list(str(record.seq).upper())
        fasta_id = str(record.id)
    for record in SeqIO.parse(args.aligned, "fasta"):
        aln_length = len(record.seq)
        aln_seq = list(str(record.seq).upper())
    if ref_length != con_length:
        logging.error(f'Length of {args.reference} sequenced ({str(ref_length)} nt) different from '
                      f'{args.consensus_subs} ({str(con_length)})')
        sys.exit(-1)
    if ref_length != aln_length:
        logging.error(f'Length of {args.reference} sequenced ({str(ref_length)} nt) different from '
                      f'{args.aligned} ({str(aln_length)})')
        sys.exit(-1)
    logging.info(f'The length of all three input sequences match ({str(ref_length)})')

    # Replacing characters to make consensus sequence with missing calls
    # Note: three characters (other than ACTG) are expected in .aligned.fa: -, N and n
    # which will become - and N in aln_seq: list(str(record.seq).upper())
    # - must be replaced by N so that all missing calls are encoded with N
    n_calls = ['-', 'N']
    new_seq = ref_seq
    num_snps = 0
    for i, (r, c) in enumerate(zip(ref_seq, con_seq)):
        if r != c:
            new_seq[i] = con_seq[i]
            num_snps += 1
            # print(str(i) + ' ' + new_seq[i] + ' ' + con_seq[i])
    num_n = 0
    for i, (r, a) in enumerate(zip(ref_seq, aln_seq)):
        if r != a:
            n_call = aln_seq[i]
            num_n += 1
            if n_call not in n_calls:
                logging.error(f'Unexpected symbol in {args.aligned}: ({str(n_call)})')
                sys.exit(-1)
            new_seq[i] = "N"
            # print(str(i) + ' ' + new_seq[i] + ' ' + aln_seq[i])
    logging.info(f'Number of consensus SNP calls initiated: {str(num_snps)}')
    logging.info(f'Number of missing SNP calls initiated: {str(num_n)}')

    # Saving consensus sequence with missing calls
    # NOTE: by default, the fasta id of .consensus.subs.fa will be used for the output fasta file, unless otherwise
    # specified with the option --output_fasta_id
    logging.info(f'Saving final consensus sequence with missing calls into {args.output}')
    if args.sample_id is not None:
        fasta_id = args.sample_id
    final_seq = SeqRecord(Seq(''.join(new_seq)), id=fasta_id, description="")
    SeqIO.write(final_seq, args.output, "fasta")

    if args.mapping_stats is not None:
        logging.info(f'Saving basic mapping stats into {args.mapping_stats}')
        output = open(args.mapping_stats, 'w')
        output.write("sample\tnum_consensus_snps\tnumber_N_calls\n")
        output.write(fasta_id + '\t' + str(num_snps) + '\t' + str(num_n) + '\n')
        output.close()


if __name__ == "__main__":
    _main()
