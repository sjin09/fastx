#!/usr/bin/env python3

__version__ = "0.0.1"
__author__ = "Sangjin Lee"


## modules
import os
import sys
import logging
from fastx.split import seq_split
from fastx.length import seq_length
from fastx.statistics import seq_statistics
from fastx.parse_arguments import parse_args

def main():
    options = parse_args(program_version=__version__)
    if not os.path.exists(os.path.abspath(options.input)):
        logging.warning("Sequence file is missing")

    if options.input.endswith((".fa", ".fq", ".fasta", ".fastq", ".fa.gz", "fq.gz", ".fofn")):
        if options.sub == "stat": ## sequence statisics
            seq_statistics(options.input, options.output)
        elif options.sub == "split": ## sequence statisics
            seq_split(options.input, options.directory)
        elif options.sub == "length": ## sequence length
            seq_length(options.input, options.output)
    else:
        logging.warning("fastx does not support the input file yet")




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
