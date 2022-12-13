#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import subprocess
import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio import AlignIO
from Bio.SeqRecord import SeqRecord
from Bio.Align.Applications import ClustalwCommandline
from Bio.Align import MultipleSeqAlignment
from Bio.Align.AlignInfo import SummaryInfo

# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------

# This script is used to is used to obtain a multiple sequence alignment and extract the
# consensus sequence from this. The script expects a FASTA file with multiple homologous
# sequences (DNA or protein) and will output the consensus sequence.

# Developmental notes:
#   Dependencies: clustalw2 and kalign need to be locally installed

# Testing notes:
#   Script tested with CLUSTAL 2.1, Clustal Omega - 1.2.4 and Kalign version 1.04


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "This script is used to obtain a multiple sequence alignment and extract the consensus sequence " \
                  "from this.\n The script expects a FASTA file with multiple homologous sequences (DNA or protein)," \
                  " and will output the MSA and consensus sequence."
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-f", "--input_fasta", action="store", dest="input_fasta",
        help="input FASTA file with multiple homologous (gene) sequences (DNA or protein)",
        required=True, metavar="FASTA")
    group.add_argument(
        "-o", "--output_prefix", action="store", dest="output_prefix",
        help="File name prefix to name output files: multiple sequence alignment and consensus sequence.",
        required=True, metavar="OUT")

    group = parser.add_argument_group('Optional MSA arguments')
    group.add_argument(
        "-m", "--msa_tool", action="store", dest="msa_tool",
        help="Multiple sequence alignment tool. To choose from: clustalw2, clustalo, kalign (default: kalign)\n"
             "NOTE: kalign will return a protein alignment regardless of input sequence type.\n"
             "NOTE: use clustalw2 or clustalo (recommended) to obtain DNA consensus sequences.\n"
             "NOTE: clustalw2 can take a long time for large alignments, recommended to use clustal omega instead.\n",
        required=False, metavar="TOOL", default="kalign"
    )
    group.add_argument(
        "-t", "--threshold", action="store", dest="threshold",
        help="The threshold value required to add a particular character in the output consensus sequence "
             "(default: 0.7)\n",
        required=False, metavar="THRESCHAR", type=float, default=0.7
    )
    group.add_argument(
        "-a", "--ambiguous_char", action="store", dest="ambiguous_char",
        help="The ambiguous character to be added in output consensus sequence when the threshold is not reached "
             "(default: N)\n",
        required=False, metavar="AMBCHAR", default="N"
    )

    return parser.parse_args()


def check_dependency(executable_name):
    """ Returns true if executable exists, else false """
    found = False
    output = subprocess.check_output(['which', executable_name]).strip()
    if output:
        found = True
    return found


def run_command_shell_string(command_line_string):
    """
    This function executes a command line, check for execution errors and but does not return stdout
    This is to be used when the stdout is not needed
    Note: shell=True needs to be set if I/O redirection operators are to be used (e.g. >) in the command line,
    otherwise they will have no special meaning, they are treated as ordinary arguments
    Note: if shell=True is used then the command line must be provided as a string, not a list
    :param command_line_string: it must be a string not a list
    """
    print('\tRunning: ' + command_line_string)
    try:
        subprocess.run(command_line_string,
                       check=True,
                       shell=True,
                       )
    except subprocess.CalledProcessError as err:
        print('ERROR:', err)


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

    # Making sure input file exist
    if not os.path.isfile(args.input_fasta):
        logging.error(f'Input FASTA file {args.input_fasta} not found!')
        sys.exit(-1)

    # Making sure choosen MSA tool within supported options
    msa_tools = ["clustalw2", "kalign", "clustalo"]
    if args.msa_tool not in msa_tools:
        logging.error(f'Chosen MSA {args.msa_tool} not supported!')
        sys.exit(-1)

    # Checking chosen MSA tools is locally installed
    dependency = args.msa_tool
    if check_dependency(dependency):
        logging.info(f'{dependency} is installed!')
    else:
        logging.error(f'{dependency} is NOT installed!')
        sys.exit(-1)

    # Other check for optional parameters
    threshold = float(args.threshold)
    ambiguous_char = args.ambiguous_char
    # Making sure threshold not > 1
    if threshold > 1 or threshold < 0.5:
        logging.error(f'Chosen parameter --threshold {str(threshold)} must be within 0.5 and 1.0!')
        sys.exit(-1)
    # Making sure ambiguous_char contains a single character
    if len(ambiguous_char) > 1:
        logging.error(f'Chosen parameter --ambiguous_char {str(ambiguous_char)} cannot be longer than 1!')
        sys.exit(-1)

    # Reading input sequences and printing information
    logging.info(f"Opening input genome {args.input_fasta}")
    input_records = SeqIO.parse(args.input_fasta, "fasta")
    record_num = 0
    for record in input_records:
        record_num += 1
    logging.info(f"Input FASTA file has {str(record_num)} sequences\n")

    # Performing multiple sequence alignment
    logging.info(f"Performing multiple sequence alignment for the {str(record_num)} sequences in {args.input_fasta}")
    # Clustalw2
    # NOTE: clustalw2 can take a long time for large alignments, recommended to use clustal omega instead
    # NOTE: clustalw2 will output an .aln file with the same prefix as the input FASTA
    if args.msa_tool == "clustalw2":
        logging.info(f"Running Clustalw2...")
        cline = ClustalwCommandline("clustalw2", infile=args.input_fasta)
        stdout, stderr = cline()
        print(stdout)

    # Clustal omega
    # NOTE: by default clustal omega will output  alignment in FASTA format
    if args.msa_tool == "clustalo":
        clustalo_msa = args.output_prefix + ".clustalo.msa.fasta"
        logging.info(f"Running clustal omega...")
        cline = ["clustalo",
                 " -i ", args.input_fasta,
                 " -o ", clustalo_msa]
        run_command_shell_string(''.join(cline))
        logging.info(f"Running clustal omega. DONE")
        check_file_exist(clustalo_msa, " clustal omega output file")

    # kalign omega
    # NOTE: kalign expected to output protein alignment in FASTA format
    if args.msa_tool == "kalign":
        kalign_msa = args.output_prefix + ".kalign.msa.fasta"
        if not os.path.isfile(kalign_msa):
            logging.info(f"Running kalign...")
            cline = ["kalign",
                     " -i ", args.input_fasta,
                     " -o ", kalign_msa]
            run_command_shell_string(''.join(cline))
            logging.info(f"Running kalign. DONE")
            check_file_exist(kalign_msa, " kalign output file")
        else:
            logging.info(f"kalign output MSA {str(kalign_msa)} already found!")

    # Obtaining consensus sequence from MSA
    align_file = args.output_prefix + ".kalign.msa.fasta"
    consensus_file = args.output_prefix + ".kalign.consensus.fasta"
    align_format = "FASTA"

    if args.msa_tool == "clustalo":
        align_file = clustalo_msa
        consensus_file = args.output_prefix + ".clustalo.consensus.fasta"
        align_format = "FASTA"

    if align_format == "FASTA":
        align = AlignIO.read(align_file, "fasta")
        logging.info(f"Printing MSA in {str(align_file)}")
        print(align)
        summary = SummaryInfo(align)
        consensus = summary.dumb_consensus(threshold, ambiguous_char)
        logging.info(f"Printing MSA consensus sequence saved in {str(consensus_file)}")
        print(consensus)
        out_seq_record = SeqRecord(consensus, id="consensus")
        SeqIO.write(out_seq_record, consensus_file, "fasta")


if __name__ == "__main__":
    _main()
