#!/usr/bin/env python3

import argparse
import logging
import os
import sys
import vcf
import subprocess


# ---------------------------------------------------------------------------------------------------------------------
# Development notes
# ---------------------------------------------------------------------------------------------------------------------
# Tested with PyVCF-0.6.8
# Tested with bcftools 1.9, tabix 1.9

# Tested with multi-sample VCF generated with bcftools-1.9 merge from multiple Snippy zipped (using bgzip/tabix)
# VCF files (snps.vcf):
#   bcftools-1.9 merge -m none -O v *.vcf.gz > input.vcf
# to do:
#   > check for intergenic region annotation

# ---------------------------------------------------------------------------------------------------------------------
# Usage notes
# ---------------------------------------------------------------------------------------------------------------------
# The input VCF must be zipped using bgzip and tabix:
#   bgzip-1.9 -c input.vcf > input.vcf.gz
#   tabix-1.9 -p vcf input.vcf.gz


# ---------------------------------------------------------------------------------------------------------------------
# Functions
# ---------------------------------------------------------------------------------------------------------------------

def parse_arguments():
    description = "Script to convert a multi-sample VCF file into a matrix and saved as a CSV file " \
                  "(samples as rows, variants as columns).\n"

    parser = argparse.ArgumentParser(description=description)

    group = parser.add_argument_group('I/O arguments')
    group.add_argument(
        "-v", "--input_vcf", action="store", dest="input_vcf",
        help="multi-sample VCF file (VCF must be zipped using bgzip and tabix).",
        required=True, metavar="INPUT_VCF"
    )
    group.add_argument(
        "-o", "--output_table", action="store", dest="output_table",
        help="File name of output variant table.",
        required=True, metavar="OUT"
    )
    group.add_argument(
        "-t", "--output_annotation_table", action="store", dest="output_annotation_table",
        help="File name of output variant annotation table.",
        required=True, metavar="OUT_ANN"
    )
    group.add_argument(
        "-b", "--bed_file", action="store", dest="bed_file",
        help="Bed file with regions to keep in output table (e.g. genes of interest).",
        required=False, metavar="BED"
    )
    group.add_argument(
        "-s", "--samples_file", action="store", dest="samples_file",
        help="File with sample ids (one sample id per line) to keep in output table.",
        required=False, metavar="SAMPLES"
    )
    group.add_argument(
        "-p", "--output_prefix", action="store", dest="output_prefix",
        help="Output prefix used to name intermediate temporary files.",
        required=False, metavar="PREFIX", default="tmp"
    )
    group.add_argument(
        "-q", "--frequency", action="store", dest="frequency",
        help="Minimum allele frequency of variants to keep in output variant table.",
        required=False, metavar="FREQ", default=0.0
    )
    group.add_argument(
        "-a", "--annotation", action="store", dest="annotation",
        help="Comma-delimited string indicating what type of variants to keep."
             "Examples: \"frameshift_variant,intergenic_region,missense_variant,stop_gained,synonymous_variant\"",
        required=False, metavar="ANN"
    )
    group.add_argument(
        "-i", "--annotation_impact", action="store", dest="annotation_impact",
        help="Comma-delimited string indicating what type of variants to keep."
             "To choose from: \"LOW,MODIFIER,MODERATE,HIGH\"",
        required=False, metavar="ANN_IMPACT"
    )

    return parser.parse_args()


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


def get_vcf_reader(my_vcf):
    if os.path.splitext(my_vcf)[-1].lower() == '.gz':
        return vcf.Reader(open(my_vcf, 'rb'))
    else:
        return vcf.Reader(open(my_vcf, 'r'))

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

    # Check dependencies exist

    # Check selected samples available on VCF

    # list of variant annotation to keep
    keep_annotations = list()
    if args.annotation is not None:
        keep_annotations = args.annotation.split(',')

    keep_annotations_impact = list()
    if args.annotation_impact is not None:
        keep_annotations_impact = args.annotation_impact.split(',')

    # Making sure input file exists
    if not os.path.exists(args.input_vcf):
        logging.error(f'Input file in --input_vcf {args.input_vcf} not found!')
        sys.exit(-1)
    if not os.path.exists(args.input_vcf + '.tbi'):
        logging.error(f'Input file in --input_vcf {args.input_vcf + ".tbi"} not found!')
        sys.exit(-1)
    if args.bed_file is not None:
        if not os.path.exists(args.bed_file):
            logging.error(f'Input file in --bed_file {args.bed_file} not found!')
            sys.exit(-1)
    if args.samples_file is not None:
        if not os.path.exists(args.samples_file):
            logging.error(f'Input file in --samples_file {args.samples_file} not found!')
            sys.exit(-1)

    # Creating a temporary copy of input VCF to work with
    tmp_vcf_gz = args.output_prefix + ".vcf.gz"
    logging.info(f"Creating a temporary copy of input VCF to work with {tmp_vcf_gz}...")
    run_command_shell_string(
        ' '.join(["cp", args.input_vcf, tmp_vcf_gz])
    )
    run_command_shell_string(
        ' '.join(["cp", args.input_vcf + '.tbi', tmp_vcf_gz + '.tbi'])
    )

    # If a samples file is selected, then keep only samples at those regions
    if args.samples_file is not None:
        logging.info(f"Keeping only samples selected in {args.samples_file}.")
        tmp_vcf2_gz = args.output_prefix + ".vcf2.gz"
        run_command_shell_string(
            ' '.join(["bcftools view", tmp_vcf_gz, "-Oz --force-samples --samples-file", args.samples_file, ">", tmp_vcf2_gz])
        )
        run_command_shell_string(' '.join(["tabix -p vcf", tmp_vcf2_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz, tmp_vcf_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz + '.tbi', tmp_vcf_gz + '.tbi']))

    # If a bed regions file is selected, then keep only variants at those regions
    if args.bed_file is not None:
        logging.info(f"Keeping only variants in selected regions {args.bed_file}.")
        tmp_vcf2_gz = args.output_prefix + ".vcf2.gz"
        run_command_shell_string(
            ' '.join(["bcftools view", tmp_vcf_gz, "-Oz --regions-file", args.bed_file, ">", tmp_vcf2_gz])
        )
        run_command_shell_string(' '.join(["tabix -p vcf", tmp_vcf2_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz, tmp_vcf_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz + '.tbi', tmp_vcf_gz + '.tbi']))

    # If AF selected > 0, then keep only variants with AF > than selected
    if float(args.frequency) > 0:
        logging.info(f"Keeping only variants with MAF greater than {str(args.frequency)}.")
        # first, re-calculating allele frequency after sub-setting VCF
        tmp_vcf2_gz = args.output_prefix + ".vcf2.gz"
        run_command_shell_string(' '.join(["bcftools +fill-tags", tmp_vcf_gz, "-Oz -o", tmp_vcf2_gz, "-- -t AF"]))
        run_command_shell_string(' '.join(["tabix -p vcf", tmp_vcf2_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz, tmp_vcf_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz + '.tbi', tmp_vcf_gz + '.tbi']))
        # filtering by AF
        tmp_vcf2_gz = args.output_prefix + ".vcf2.gz"
        run_command_shell_string(
            ''.join(["bcftools view ", tmp_vcf_gz, " -q ", str(args.frequency), ":nref", " -Oz -o ", tmp_vcf2_gz])
        )
        run_command_shell_string(' '.join(["tabix -p vcf", tmp_vcf2_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz, tmp_vcf_gz]))
        run_command_shell_string(' '.join(["mv", tmp_vcf2_gz + '.tbi', tmp_vcf_gz + '.tbi']))

    # Creating genotype table from VCF file
    genotypes = dict()  # dict to save sample ids, variant ids and their genotypes
    annotation_dict = dict()  # dict to save variant ids and their annotation
    logging.info(f"Reading VCF file {tmp_vcf_gz} to create table with variants at selected regions/samples")
    vcf_reader = get_vcf_reader(tmp_vcf_gz)
    vcf_samples = vcf_reader.samples
    var_ids = list()
    print("\tNumber of samples extracted from VCF: " + str(len(vcf_samples)))
    for record in vcf_reader:
        var_id = record.CHROM + "." + record.REF + "." + str(record.POS) + "." + str(record.ALT[0])
        # SnpEff annotation string:
        snpeff_ann = record.INFO['ANN'][0]
        # Allele | Annotation | Annotation_Impact | Gene_Name | Gene_ID | Feature_Type | Feature_ID |
        # Transcript_BioType | Rank | HGVS.c | HGVS.p | cDNA.pos / cDNA.length | CDS.pos / CDS.length |
        # AA.pos / AA.length | Distance | ERRORS / WARNINGS / INFO'"
        (allele, annotation, annotation_impact, gene_name, gene_id, feature_type, feature_id, transcript_biotype, rank,
         hgvs_c, hgvs_p, cdna_length, ccds_pos, cds_pos, aa_pos, distance, *_) = snpeff_ann.split('|')
        # Saving variant annotation
        ann = 'NA'
        if transcript_biotype == "protein_coding":
            ann = gene_name + '\t' + hgvs_p + '\t' + annotation + '\t' + annotation_impact
        else:
            ann = gene_name + '\t' + hgvs_c + '\t' + annotation + '\t' + annotation_impact

        # Filtering by variant annotation if chosen by the user
        to_keep = True
        if args.annotation is not None:
            to_keep = False
            if annotation in keep_annotations:
                to_keep = True
        if args.annotation_impact is not None:
            to_keep = False
            if annotation_impact in keep_annotations_impact:
                to_keep = True

        # Saving is kept after filtering
        if to_keep is True:
            print(var_id)
            var_ids.append(var_id)
            annotation_dict[var_id] = ann
            # saving genotype calls
            for call in record.samples:
                if call.gt_type is not None:
                    if call.sample not in genotypes:
                        genotypes[call.sample] = dict()
                    genotypes[call.sample][var_id] = str(call.gt_type)

    # Check
    if len(var_ids) == 0:
        logging.error(f"No variants kept. Check chosen filters.")
        sys.exit(-1)

    # Printing output genotype table
    logging.info(f"Printing output table into {args.output_table}")
    output = open(args.output_table, 'w')
    header = "sample"
    for var_id in var_ids:
        header = header + '\t' + var_id
    output.write(header + '\n')

    for sample in vcf_samples:
        newline = sample
        for var_id in var_ids:
            genotype = "0"
            if sample in genotypes:
                genotype = "0"
                if var_id in genotypes[sample]:
                    genotype = genotypes[sample][var_id]
            newline = newline + "\t" + genotype
        output.write(newline + '\n')
    output.close()

    # Printing output genotype table
    logging.info(f"Printing output annotation table into {args.output_annotation_table}")
    output2 = open(args.output_annotation_table, 'w')
    for var_id in var_ids:
        newline = var_id + '\t' + annotation_dict[var_id]
        output2.write(newline + '\n')
    output2.close()


if __name__ == "__main__":
    _main()