# modules
import argparse
import sys


# argparse
def parse_args(program_version, arguments=sys.argv[1:]):

    # main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="fastx reads FASTA, FASTQ, gzipped FASTA or gzipped FASTQ and returns a desired output",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )

    # subcommands
    subparsers = parser.add_subparsers(dest="sub")
    # subcommands: gap
    parser_gap = subparsers.add_parser(
        "gap", help="returns a BED file with gap positions",
    )
    parser_gap.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, or gzipped FASTA \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz",
    )
    parser_gap.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="BED file to return the sequences",
    )

    # subcommands: head
    parser_head = subparsers.add_parser(
        "head", help="returns the first n lines of sequences",
    )
    parser_head.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.,gz .fastq, .fastq.gz",
    )
    parser_head.add_argument(
        "-n",
        "--number",
        type=int,
        default=5,
        required=True,
        help="number of sequences to return",
    )
    parser_head.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return the sequences",
    )

    # subcommands: sort
    parser_sort = subparsers.add_parser("sort", help="returns sorted sequences",)
    parser_sort.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fasta, .fa.gz \
        FASTQ: .fq, .fastq, .fq.gz",
    )
    parser_sort.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return sorted sequences",
    )

    # subcommands: sequence statistics
    parser_statistics = subparsers.add_parser(
        "stat", help="reads returns sequence statistics"
    )
    parser_statistics.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )
    parser_statistics.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .stat file(s)"
    )
    # subcommands: pyfastx index
    parser_index = subparsers.add_parser(
        "index", help="generate .fxi index files"
    )
    parser_index.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )

    # subcommands: names
    parser_names = subparsers.add_parser(
        "names", help="returns sequence id"
    )
    parser_names.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )
    parser_names.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return seq id"
    )
    # subcommands: split
    parser_split = subparsers.add_parser(
        "split", help="splits and return sequences"
    )
    parser_split.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fasta.gz"
    )
    parser_split.add_argument(
        "-d",
        "--dir",
        type=str,
        required=True,
        help="directory to return the split sequences"
    )

    # subcommands: length
    parser_length = subparsers.add_parser(
        "length", help="returns sequence id and length as tab separated value"
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )
    parser_length.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .length file(s)"
    )
    # subcommands: target
    parser_target = subparsers.add_parser(
        "target", help="returns target contigs for somatic mutation calling"
    )
    parser_target.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, or gzipped FASTA \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz"
    )
    parser_target.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return target contigs"
    )
    ## subcommands: blacklist
    parser_blacklist = subparsers.add_parser(
        "blacklist", help="returns filtered.fastq based on blacklist"
    )
    parser_blacklist.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ, gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )
    parser_blacklist.add_argument(
        "--blacklist",
        type=str,
        required=True,
        help="ZMW blacklist"
    )
    parser_blacklist.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return filtered FASTQ file"
    )
    parser_blacklist = subparsers.add_parser(
        "whitelist", help="returns whitelisted fastq"
    )
    parser_blacklist.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ, gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz"
    )
    parser_blacklist.add_argument(
        "--whitelist",
        type=str,
        required=True,
        help="ZMW whitelist"
    )
    parser_blacklist.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return filtered FASTQ file"
    )


    # subcommand: fasta2fastq
    parser_fastq = subparsers.add_parser("fastqc", help="sequence quality control",)
    parser_fastq.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )

    # subcommand: fasta2fastq
    parser_fasta2fastq = subparsers.add_parser(
        "fasta2fastq", help="converts FASTA to FASTQ file",
    )
    parser_fasta2fastq.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA or gzipped FASTA \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz",
    )
    parser_fasta2fastq.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return FASTA file",
    )

    ## subcommands: fastq2fasta
    parser_fastq2fasta = subparsers.add_parser(
        "fastq2fasta", help="converts FASTQ to FASTA file",
    )
    parser_fastq2fasta.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ, gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_fastq2fasta.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return FASTA file",
    )

    ## subcommandS: tricounts
    parser_tricounts = subparsers.add_parser(
        "tricounts",
        help="counts trinucleotide sequence contexts from FASTA and FASTQ file",
    )
    parser_tricounts.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA and FASTQ files \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_tricounts.add_argument(
        "-t",
        "--threshold",
        default=0,
        required=False,
        type=int,
        help="sequence tricounts threshold",
    )
    parser_tricounts.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return trinulceotide context counts",
    )

    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser, parser.parse_args(arguments)
