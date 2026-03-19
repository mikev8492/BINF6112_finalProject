#!/usr/bin/env python3
import sys, argparse
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
    motif_parser = argparse.ArgumentParser(description="")

    motif_parser.add_argument("-s", "--sequence_filepath",
                            type=str, default="inputs/test/pUC19.fasta", 
                            help="Path to the input sequence file. Accepted formats: FASTA (fasta, .fas, .fa, .fna, .ffn, .faa, .mpfa, .frn) or GenBank (.gb, .gbk).")

    return motif_parser.parse_args()


def validate_arguments(args: argparse.Namespace) -> bool:
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
    genbank_extensions = ('.gb', '.gbk')
    fasta_extensions = ('.fasta', '.fas', '.fa', '.fna', '.ffn', '.faa', '.mpfa', '.frn')

    if not args.sequence_filepath.split(".")[-1] in genbank_extensions and not args.sequence_filepath.split(".")[-1] in fasta_extensions:
        raise ValueError("Error: Invalid file type. Please provide a .fasta, .fa, or .txt file.", 'red')
        return False


    return True


def create_db():
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
                # print(f"{enzyme}: {enzyme.site},  {enzyme.elucidate()}")
                file.write(f"{enzyme},{enzyme.site},{enzyme.elucidate()}\n")


def main():

# Could this be used to make it easier for the input.py class to determine the file type?    
    args = create_parser()
    # create_db()


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


    # Print to test output
    print(f"\nPlasmid:\n{plasmid}\n\nEnzyme list:")
    for enzyme, info in enzymes.items():
        print(f"\t{enzyme}: {info}")


if __name__=="__main__":
    main()
