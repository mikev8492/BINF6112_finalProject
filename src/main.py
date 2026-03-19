#!/usr/bin/env python3
import sys, argparse
from Bio.Restriction import CommOnly
from motif_id_lib.input import Sequence, Enzymes

"""
This file contains the main functions that will perform the logic for the program. 
"""


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
    create_db()

    # load and parse sequence 
    user_file = Sequence(sys.argv[1])
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
