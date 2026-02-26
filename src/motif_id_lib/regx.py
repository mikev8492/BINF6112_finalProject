class RegX
      
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
    
    
