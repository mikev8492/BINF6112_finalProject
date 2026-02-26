import input 
import sys

class RegX

  def __init__(self, enzyme_dict):
    self.enzyme_dict = enzyme_dict
    self.updated_dict = {}

  def dict_check(self) -> dict:
    '''
    Purpose: Checks the argument is a dictionary before proceeding
    Input: Dictonary
    Output: Updated dictionary
    '''
    if not isinstance(enzyme_dict, dict):
      sys.stderr("Error: Argument is not a dictionary")
    else:
      sys.stdout("Dictionary accepted...")
      
  def motif_regx(self):
    '''
    Purpose: This function generates a regular expression for all selected enzyme motifs and add it to the dictionary.
    Input: Dictionary generated from load_REs()
    Output: Dictionary with regular expression added for each enzyme
    High level steps: 
    - Go through each enzyme entry in the original dictionary and look at the motif
    - Generate a regular expression for the motif and append to the list of values for that enzyme
    - Return updated dictionary
    '''
    
    
