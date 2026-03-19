# BINF6112_finalProject
**UNCC BINF6112 - Programming II final project.**

## License: 
**GNU General Public License Version 3**

The GNU GPL is a license that ensures code is open-source. GNU GPL allows others to utilize, modify, or distribute code. If other users modify the code, then these users are expected to share their changes to the code under a GNU GPL to maintain the open-source integrity of the code.

We chose to use the GNU GPL to make our code easily accessible for anyone to use, or to further build upon. 

[Project URL](https://github.com/mikev8492/BINF6112_finalProject)

## Group members:
**Michael Villarreal**
mvillar6@charlotte.edu
UNCC ID: 801504997

**Kayla Ball**
kwill423@charlotte.edu
UNCC ID: 801489303

**Bobby Luker**
rluker@charlotte.edu
UNCC ID: 801484356

## Project File Structure:
```
└── 📁BINF6112_finalProject
    └── 📁inputs
        └── 📁test
            ├── pUC19.fasta
            ├── pUC19.gb
    └── 📁results
        ├── plasmid_map.txt
    └── 📁src
        └── 📁database
            ├── enzymes.csv
        └── 📁motif_id_lib
            └── 📁__pycache__
            ├── __init__.py
            ├── id_motifs.py
            ├── input.py
            ├── output.py
            ├── regx.py
        ├── main.py
    ├── dependencies.txt
    ├── environment_mac.yml
    ├── environment.yml
    ├── ideas.md
    ├── LICENSE
    ├── Pseudocode.py
    └── README.md
```
## Development Log:
Date: 3/19/2026
- updated relative filepaths to run main.py from root project directory.
- Added `environment-alternative.yml` to resolve cross compatibility conflicts. 
- Refactored `main.py`: 
    1. Moved CLI header to `input.py` module to allow User interface as optional argument.
    2. Added `ArgParse` functionality to simplify CLI.  

## Testing Instructions:
Date: 3/19/2026

1. Create the environment:

For linux:
```bash
conda env create -f environment.yml
```
For MacOS:
```bash
conda env create -f environment-alternative.yml
```
Conda will automatically create an environment named finalproj with all the specified packages and versions.

**Note:** if environment does not install on MacOS, install dependencies listed in `dependencies.txt` manually and run step 3. using installed python version. 

2. Activate the environment:
```bash
conda activate finalproj
```
3. Run following command to test:
```bash
python src/main.py
```

4. Use arrows to move down checkbox list and click "space" to select enzymes. 

Expected Output:

- Terminal interface displays with application header and checkbox for user to select restriction enzymes to cut with
- prints Plasmid object (list[header, sequence]) and enzyme object (dict{enzyme: [motif, cut]})


## Instructions:
1. Place your plasmid sequence file in the `inputs` folder (GenBank or FASTA format).
2. Run the following script with your sequence filename:
3. Check which restriction enzymes you want to identify within plasmid sequence using CLI
4. Review output located in results folder


## Overview:
REcut is a python3 tool that takes a plasmid sequence file as input and generates an annotated sequence map with Restriction Enzyme cut sites. 

The `src` file contains the following:
1. `database` : Contains a CSV file of restriction enzymes and recognition sequence motifs. 
2. `main.py`: Contains the logic for the program.
    - Loads the user input sequence and restriction enzymes database. 
    - Parses the plasmid sequence with each enzyme motif to identify cut locations.
    - Produces output containing an annotated map of the sequence with the enzyme cut sites in the `results` folder.
3. `motif_id_lib`: Contains modules designated for each function within our program. These modules are:
    - id_motifs.py: Group member responsible- Bobby
    - input.py: Group member responsible- Michael
    - output.py: Group member responsible- Michael
    - regx.py: Group member responsible- Kayla
## Troubleshooting:
1. User inputs incorrect file path and incorrect file type - Check that the file/file path exists, is readable, and is the correct file format (FASTA or GenBank file)
2. No matches found and updated dictionary returns None - Error will occur and state that no matches were identified and exists program
