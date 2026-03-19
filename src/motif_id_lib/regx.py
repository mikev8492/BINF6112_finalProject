class RegX
      
  def motif_regx(self):
    '''
    Purpose
    -------
    Generates a regular expression for all selected enzyme motifs and add it to the dictionary.

    Args
    ----
    Dictionary generated from load_REs()

    Returns
    -------
    Dictionary with regular expression added for each enzyme

    High level steps: 
    - Go through each enzyme entry in the original dictionary and look at the motif
    - Generate a regular expression for the motif and append to the list of values for that enzyme
    - Return updated dictionary
    '''
    
    
    
