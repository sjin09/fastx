#!/usr/bin/env python3

__version__ = "0.0.1"
__author__ = "Sangjin Lee"


# modules
import os
import sys
import logging
from fastx.gap import seq_gaps
from fastx.head import seq_head
from fastx.sort import seq_sort
from fastx.split import seq_split
from fastx.names import seq_names
from fastx.index import seq_index
from fastx.length import seq_length
from fastx.target import get_target
from fastx.fastqc import seq_fastqc
from fastx.tricounts import tricounts
from fastx.parse_args import parse_args
from fastx.statistics import seq_statistics
from fastx.util import fasta2fastq, fastq2fasta, return_blacklist, return_whitelist

# global
INFILE_SUFFIX = (
    ".fa",
    ".fq",
    ".fa.gz",
    ".fq.gz",
    ".fasta",
    ".fastq",
    ".fasta.gz",
    ".fastq.gz",
)


def run_subcommands(parser, options):
    if options.sub == "gap":  # return first n lines of sequences
        seq_gaps(options.input, options.output)
    elif options.sub == "head":  # return first n lines of sequences
        seq_head(options.input, options.number, options.output)
    elif options.sub == "stat":  # returns sequence statisics
        seq_statistics(options.input, options.output)
    elif options.sub == "sort":  # sort sequences
        seq_sort(options.input, options.output)
    elif options.sub == "names":  # return first n lines of sequences
        seq_names(options.input, options.output)
    elif options.sub == "split":  # split and return sequences
        seq_split(options.input, options.directory)
    elif options.sub == "index":  # generate fxi index files
        seq_index(options.input)
    elif options.sub == "length":  # return sequence id and length
        seq_length(options.input, options.output)
    elif options.sub == "target":  # return list of contigs for somatic mutation calling
        get_target(options.input, options.output)
    elif options.sub == "fastqc":  # fastqc
        seq_fastqc(options.input)
    elif options.sub == "blacklist":  # generate fxi index files
        return_blacklist(options.input, options.blacklist, options.output)
    elif options.sub == "whitelist":  # generate fxi index files
        return_whitelist(options.input, options.whitelist, options.output)
    elif options.sub == "fasta2fastq":  # fasta2fastq
        fasta2fastq(options.input, options.output)
    elif options.sub == "fastq2fasta":  # fastq2fasta
        fastq2fasta(options.input, options.output)
    elif options.sub == "tricounts":  # trinucleotide sequence context counts
        tricounts(options.input, options.threshold, options.output)
    else:
        print("The subcommand does not exist!\n")
        parser.print_help()
        sys.exit(0)


def main():
    parser, options = parse_args(program_version=__version__)
    if not os.path.exists(os.path.abspath(options.input)):
        logging.warning("Sequence file is missing")
        sys.exit(0)
    elif not options.input.endswith(INFILE_SUFFIX):
        logging.warning("fastx does not support the provided input file")
        sys.exit(0)
    else:
        run_subcommands(parser, options)


if __name__ == "__main__":
    main()
