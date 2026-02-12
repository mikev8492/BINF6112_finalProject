#!/usr/bin/env python3

"""
This file contains the main functions that will perform the logic for the program. 
"""

def load_REs(file:argv) -> user_enzymes:dict and sequence:str:
   '''
   Read in user sequence file and determine if GenBank or FASTA file. Function has designated workflow dependent on file type. Function will open enzyme database and 
   parse through for selected enzymes. If no enzymes are selected, a default dictionary will be generated instead.

   Parameters
   ----------
   Input file
   User selected enzymes via CLI

   Returns
   -------
   A dictionary that contains selected or default enzymes with motifs 
   Input sequence as a string

   Raises
   ------
   IOError will occur if user does not input file or input file is incorrect format
   '''
  pass

def motif_regx(motifs:dict) -> dict:
  
  pass  

def identify_motifs(sequence:str, regx:dict) -> lst:

  pass
    
def generate_RE_map(sequence:str, regx:dict, identify_motifs:lst) -> file:argv:

  pass
    

def main():


if __name__=="__main__":
    main()
