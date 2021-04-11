#!/usr/bin/env python3

__version__ = "0.0.1"
__author__ = "Sangjin Lee"


## modules
import os
import sys
import logging
from fastx.gap import seq_gaps
from fastx.head import seq_head
from fastx.sort import seq_sort
from fastx.split import seq_split
from fastx.length import seq_length
from fastx.fastqc import seq_fastqc
from fastx.statistics import seq_statistics
from fastx.parse_arguments import parse_args
from fastx.transform import fasta2fastq, fastq2fasta

def main():
    parser, options = parse_args(program_version=__version__)
    if not os.path.exists(os.path.abspath(options.input)):
        logging.warning("Sequence file is missing")

    if options.input.endswith((".fa", ".fq", ".fasta", ".fastq", ".fa.gz", "fq.gz", ".fofn")):
        if options.sub == "gap": ## return first n lines of sequences
            seq_gaps(options.input, options.output)
        elif options.sub == "head": ## return first n lines of sequences
            seq_head(options.input, options.number, options.output)
        elif options.sub == "stat": ## returns sequence statisics
            seq_statistics(options.input, options.output)
        elif options.sub == "sort": ## sort sequences
            seq_sort(options.input, options.output)
        elif options.sub == "split": ## split and return sequences
            seq_split(options.input, options.directory)
        elif options.sub == "length": ## return sequence id and length
            seq_length(options.input, options.output)
        elif options.sub == "fastqc": ## fasta2fastq
            seq_fastqc(options.input, options.prefix)
        elif options.sub == "fasta2fastq": ## fasta2fastq
            fasta2fastq(options.input, options.output)
        elif options.sub == "fastq2fasta": ## fastq2fasta
            fastq2fasta(options.input, options.output)
        else:
            print("The subcommand does not exist!\n")
            parser.print_help()
            parser.exit()
    else:
        logging.warning("fastx does not support the provided input file")




    # if len(sys.argv) != 1: 
    #     pass
    # else:
    #     options.print_help()
        # logging.info("INPUT: {0}".format(os.path.abspath(options.reads)))
        # logging.info("GENOME: {0}".format(os.path.abspath(options.genome)))
    # logFormatter = logging.Formatter("%(asctime)s [%(levelname)-7.7s]  %(message)s")
    # rootLogger = logging.getLogger()
    # logging.info("****************** Start SVIM, version {0} ******************".format(__version__))
    # logging.info("CMD: python3 {0}".format(" ".join(sys.argv)))
    # logging.info("WORKING DIR: {0}".format(os.path.abspath(options.working_dir)))


if __name__ == "__main__":
    main()
    sys.exit(0)
