#!/usr/bin/env python3

import argparse
import logging
import os
import sys
from Bio import SeqIO
from Bio import SeqFeature as sf
from BCBio import GFF

# ------------------------------------------------------------------------------------
# Notes
# ------------------------------------------------------------------------------------

# At the moment only GFF format is supported.
# Remove first "source" GFF line, e.g.: "MW2	EMBL/GenBank/SwissProt	source	1	2820462	"


# ------------------------------------------------------------------------------------
# Functions
# ------------------------------------------------------------------------------------

def parse_arguments():
    description = "Python script to parse GFF file to table format"
    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('required arguments')
    group.add_argument(
        "-i", "--input_gff", action="store", dest="input_gff",
        help="input genome annotation in GFF format",
        required=True, metavar="INPUT_GFF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="output genome annotation in table format",
        required=True, metavar="OUTPUT_TABLE"
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
    if not os.path.isfile(args.input_gff):
        logging.error(f'Input genome {args.input_gff} not found!')
        sys.exit(-1)

    # Defining feature types and attributes to extract
    feature_types = ['CDS']
    attributes = ['locus_tag', 'gene', 'product']

    # Creating table header
    header = attributes + ['type', 'start', 'end', 'strand', 'snpEff_locus_id']

    # Loading input genome
    logging.info(f"Opening input annotated genome {args.input_gff}")
    in_handle = open(args.input_gff)
    gff_tab_file = open(args.output_table, 'w')
    gff_tab_file.write('\t'.join(header)+'\n')
    for rec in GFF.parse(in_handle):
        for feature in rec.features:
            if feature.type in feature_types:
                table_line = ''
                print(feature)
                for attribute in attributes:
                    if attribute in feature.qualifiers.keys():
                        table_line = table_line + feature.qualifiers.get(attribute)[0] + '\t'
                    else:
                        table_line = table_line + '-' + '\t'
                table_line = table_line + feature.type + '\t' + str(feature.location.start) + '\t' \
                         + str(feature.location.end) + '\t' + str(feature.location.strand)
                # Creating snpEff locus id: e.g. CDS_EMRSA15_938198_938683
                # Note: feature.location.start extracted is shifted 1bp leftwards, so +1 need to be added
                snpEff_locus_id = feature.type + '_' + rec.id + '_' \
                              + str(feature.location.start+1) + '_' + str(feature.location.end)
                table_line = table_line + '\t' + snpEff_locus_id + '\n'
                print(table_line)
                gff_tab_file.write(table_line)
    in_handle.close()
    gff_tab_file.close()


if __name__ == "__main__":
    _main()
