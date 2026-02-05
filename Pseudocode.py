'''
The following file describes the objectives of each function and the corresponding steps to achieve the correct output.
The link below is a known plasmid that we can use for testing our code:
https://www.ncbi.nlm.nih.gov/nuccore/U13852.1 
'''
def load_REs(file):
  '''
  Purpose: This function is asking the user to select restriction enzymes it wants to locate in the plasmid sequence.
  Input: User input at CLI
  Output: A dictionary that contains the name of the enzyme as the key and the motif, cut position, and type of cut (sticky or blunt ends) in a list from the file or the default
  dictionary if user selects none.
  High level steps: 
    - Display list of enzymes
    - Receive user input of which enzymes to look for within the sequence or use default list if none selected
    - Read in enzymes.csv and extract the selected enzyme and corresponding motif
    - Construct a dictionary with the extracted information
  '''
def motif_regx():
  '''
  Purpose: This function generates a regular expression for all selected enzyme motifs and add it to the dictionary.
  Input: Dictionary generated from load_REs()
  Output: Dictionary with regular expression added for each enzyme
  High level steps: 
  - Go through each enzyme entry in the original dictionary and look at the motif
  - Generate a regular expression for the motif and append to the list of values for that enzyme
  - Return updated dictionary
  '''
def identify_motifs(sequence <=5kb,dictionary):
  '''
  Purpose: This function takes an input sequence less than 5kb long and parses through the sequence to identify motifs within the dictionary.
  Input: A valid sequence and the dictionary generated from motif_regx().
  Output: A list of the names of the enzymes whose motifs were found in the sequence
  High level steps:
    - Check sequence is less than the maximum; if yes continue, if not exit the function with error message
    - Check the sequence for possible motifs found within the dictionary
    - When the program comes across a motif match, it will add the matched start index and the enzyme name to a list called found_REs
    - When finished, it will return the complete list found_REs
  '''
def generate_RE_map():
  '''
  Purpose: This function will generate a visual for the user to see where in the sequence a restriction enzyme is found, the cut point
  within the motif, and which enzyme it is.
  Input: Dictionary, list, sequence
  Output: 1) Table with the name of the enzyme, it's double stranded cut site, index of the cut site as an output file
          2) The full sequence with annotated cut sites*
  High level steps:
    - Generate .csv file and add in information from list and enzyme file 
    - *Unsure specific steps; we want to have the sequence shown with the cut sites and name of enzymes available
  '''

