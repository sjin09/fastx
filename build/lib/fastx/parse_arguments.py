## modules
import sys
import argparse

## argparse
def parse_args(program_version, arguments=sys.argv[1:]):

    ## main_arguments
    parser = argparse.ArgumentParser(
        add_help=True,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        description="fastx reads SAM/BAM, FASTA, FASTQ, gzipped FASTA or gzipped FASTQ and returns a desired output",
    )
    parser.add_argument(
        "-v",
        "--version",
        action="version",
        version="%(prog)s {version}".format(version=program_version),
    )
    
    ## subcommands
    subparsers = parser.add_subparsers(help="fastx subcommands", dest="sub")

    ## subcommands: statistics
    parser_statistics = subparsers.add_parser(
        "stat",
        help="reads FASTA, FASTQ, gzipped FASTA or gzipped FASTA and returns sequence statistics",
    )
    parser_statistics.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fasta, .fa.gz \
        FASTQ: .fq, .fastq, .fq.gz"
    )
    parser_statistics.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .length file(s)"
    )
    
    ## subcommands: length
    parser_length = subparsers.add_parser(
        "length",
        help="reads FASTA, FASTQ, gzipped FASTA or gzipped FASTA and returns sequence id and length as a tsv",
    )
    parser_length.add_argument(
        "-i",
        "--input",
        type=str,
        required=True,
        help="FASTA, FASTQ, gzipped FASTA, or gzipped FASTQ \
        Programs supports the following prefixes: \
        FASTA: .fa, .fasta, .fa.gz \
        FASTQ: .fq, .fastq, .fq.gz"
    )
    parser_length.add_argument(
        "-o",
        "--output",
        required=True,
        type=argparse.FileType("w"),
        help="FILE to return .length file(s)"
    
    )

    if len(arguments) == 0:
        parser.print_help()
        parser.exit()
    else:
        return parser.parse_args(arguments)


