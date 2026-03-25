#!/usr/bin/env python3
from genericpath import isfile
import sys, argparse
from pathlib import Path
from Bio.Restriction import CommOnly
from motif_id_lib.input import Sequence, Enzymes

"""
This file contains the main functions that will perform the logic for the program. 
"""

def create_parser() -> argparse.Namespace:
    '''
    Description: 
        Defines each CLI for the program

    Args:
        None

    Returns:

    Raises:
        None
    '''
    motif_parser = argparse.ArgumentParser(description="REcut: Plasmid Sequence Cutting Tool - Generates an annotated plasmid map with restriction enzyme cut sites.",
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    motif_parser.add_argument(
        "-s", "--sequence_filepath",
        type=str, 
        default="inputs/test/pUC19.fasta", 
        help="Path to the input sequence file. Accepted formats: FASTA (fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn) or GenBank (.gb, .gbk)."
    )

    motif_parser.add_argument(
        "-o", "--output",
        type=str,
        default="results/plasmid_map",
        help="Output filepath"
    )

    motif_parser.add_argument(
        "--interface",
        action="store_true",
        help="Enable TUI (terminal user interface) to allow manual selection of Restriction enzymes to search with."
    )


    return motif_parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> None:
    '''
    Description: 
        Validate the command-line arguments

    Args:
        args (argparse.Namespace): Parsed command-line arguments

    Returns:
        bool: True if arguments are valid, False otherwise

    Raises:
        ValueError: If any argument is invalid
    '''
    # Check filetype
    genbank_extensions = ('gb', 'gbk')
    fasta_extensions = ('fasta', 'fas', 'fa', 'fna', 'ffn', 'faa', 'mpfa', 'frn')

    input_extension = args.sequence_filepath.split(".")[-1]
    
    if input_extension not in genbank_extensions and input_extension not in fasta_extensions:
        raise ValueError("Invalid file type (-s, --sequence_filepath). Please provide a FASTA or GenBank file.")

    # Check plasmid filepath exists
    filepath = Path(args.sequence_filepath)
    if not filepath.is_file():
        raise IOError("Filepath (-s, --sequence_filepath) does not exist")


def create_db():
    """
    Description: 
        Create restriction enzyme database file (CSV). Uses "CommOnly" (common type II enzymes) list of REBASE data (Bio.Restriction).
    """
    # Check if database exists:
    filepath = Path("src/database/enzymes.csv")
    if not filepath.is_file():
        sys.stdout.write(f"Creating database {filepath}\n")

        starting_enzyme_list = [
            "EcoRI",
            "HindIII",
            "BamHI",
            "XhoI",
            "NotI",
            "SalI",
            "PstI",
            "KpnI",
            "XbaI",
            "EcoRV",
            "SmaI",
            "NdeI",
            "SacI",
            "SpeI",
            "BglII",
            "ApaI",
            "SphI",
            "MluI",
            "ClaI",
            "HaeIII"
        ]

        with open("src/database/enzymes.csv", "w") as file:
            file.write("Enzyme,motif,cutInfo\n")
            for enzyme in CommOnly:
                if enzyme.__name__ in starting_enzyme_list:
                    file.write(f"{enzyme},{enzyme.site},{enzyme.elucidate()}\n")
    # sys.stdout.write(f"Database {filepath} ready\n")

def main():

    try:
        create_db()

        args = create_parser()
        validate_arguments(args)

        # =======================
        # INPUT.py
        # =======================
        # load and parse sequence 
        user_file = Sequence(args.sequence_filepath)
        user_file.load_sequence()

        # plasmid: list [header, seq]
        plasmid = user_file.sequence


        # load and filter enzymes
        re_list = Enzymes()
        re_list.load_REs()
        re_list.interface()
        re_list.filter_enzymes()
        # Enzymes dict name: [motif, cut]
        enzymes = re_list.filtered
        # =======================

        # Print test output: Plasmid sequence, Enzyme list
        print(f"\nPlasmid:\n{plasmid}\n\nEnzyme list:")
        for enzyme, info in enzymes.items():
            print(f"\t{enzyme}: {info}")


    # Error handling:
    except ValueError as err:
        sys.stderr.write(f"Error: {err}\n")
        return 1
    except IOError as err:
        sys.stderr.write(f"Error: {err}\n")
        return 1

if __name__=="__main__":
    sys.exit(main())
