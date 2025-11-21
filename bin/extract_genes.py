#!/usr/bin/env python3
import argparse
from Bio import SeqIO

def parse_args(args=None):
    Description = "Extract genes from FASTA files based on provided coordinates."
    Epilog = """Example usage: python extract_genes.py --fasta sample.fasta --gff sample.gff --interest_genes "gene1,gene2;gene3" --output_prefix "prefix_" """

    parser = argparse.ArgumentParser(description=Description, epilog=Epilog)
    parser.add_argument(
        "-f",
        "--fasta",
        type=str,
        required=True,
        help="Path to input FASTA file (required).",
    )
    parser.add_argument(
        "-o",
        "--output_prefix",
        type=str,
        required=False,
        help="Prefix for output FASTA files with extracted genes.",
    )
    return parser.parse_args(args)




def main(args=None):
    args = parse_args(args)


if __name__ == "__main__":
    main()