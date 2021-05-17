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

    parser_head = subparsers.add_parser(
        "gap", help="returns a BED file with gap positions",
    )
    parser_head.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, or gzipped FASTA \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz",
    )
    parser_head.add_argument(
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
    parser_statistics = subparsers.add_parser("sort", help="returns sorted sequences",)
    parser_statistics.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fasta, .fa.gz \
        FASTQ: .fq, .fastq, .fq.gz",
    )
    parser_statistics.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return sorted sequences",
    )

    # subcommands: sequence statistics
    parser_statistics = subparsers.add_parser(
        "stat", help="reads returns sequence statistics",
    )
    parser_statistics.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_statistics.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .stat file(s)",
    )

    # subcommands: split
    parser_statistics = subparsers.add_parser(
        "split", help="splits and return sequences",
    )
    parser_statistics.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fasta.gz",
    )
    parser_statistics.add_argument(
        "-d",
        "--directory",
        type=str,
        required=True,
        help="directory to return the split sequences",
    )

    # subcommands: length
    parser_length = subparsers.add_parser(
        "length", help="returns sequence id and length as tab separated value",
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_length.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .length file(s)",
    )
    # subcommand: fasta2fastq
    parser_length = subparsers.add_parser("fastqc", help="sequence quality control",)
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    # subcommand: fasta2fastq
    parser_length = subparsers.add_parser(
        "fasta2fastq", help="converts FASTA to FASTQ file",
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA or gzipped FASTA \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz",
    )
    parser_length.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return FASTA file",
    )

    parser_length = subparsers.add_parser(
        "fastq2fasta", help="converts FASTQ to FASTA file",
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTQ, gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_length.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return FASTA file",
    )

    parser_length = subparsers.add_parser(
        "tricounts",
        help="counts trinucleotide sequence contexts from FASTA and FASTQ file",
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA and FASTQ files \
        Programs supports the following prefixes: \
        FASTA: .fa, .fa.gz, .fasta, .fasta.gz \
        FASTQ: .fq, .fq.gz, .fastq, .fastq.gz",
    )
    parser_length.add_argument(
        "-t",
        "--threshold",
        default=0,
        required=False,
        type=int,
        help="sequence length threshold",
    )
    parser_length.add_argument(
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
