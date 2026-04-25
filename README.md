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
        ├── plasmid_results.csv
    └── 📁src
        └── 📁database
            ├── enzymes.csv
        └── 📁motif_id_lib
            └── 📁__pycache__
            ├── __init__.py
            ├── csv_output.py
            ├── input.py
            ├── motif_locator.py
            ├── output.py
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

Date 04/17/2026
- Updated `csv_output.py` module:
    - Added `CreateCSV` class to create the CSV output file for the user.
        - `create_csv_output` function writes the enzyme name, motif, cut site, observed count, and the start position locations as a row in the CSV file. 
- Refactored `main.py`:
    - Updated function calls in `main()` to create the CSV output file. 

Date: 04/09/2026
- Updated environment to include `matplotlib`
- Added `annotate_plasmid` function to  `output.py` module:
    - Uses results dictionary and plasmid sequence string to generate a circular map.
    - A circle is drawn with radius of 1.0, and the enzyme cut locations are converted from a bp to an angle on the circle, starting from the top:
        angle = 90° − (bp position / sequence length) × 360°
    - A tick mark is then drawn at the calculated angle. 
    - png graph is generated and saved to `results/plasmid_map.png`

Date: 04/01/2026
- Added `Motifs` class to `motif_locator.py` module:
    - `Motifs` class locates all instances of an enzyme motif in the plasmid sequence. 
    - Class functions added:
        - `array_set`: transforms the plasmid sequence to a numpy array.
        - `motif_search`: searches for all locations of an enzyme motif in the plasmid sequence
        - `get_motif_results`: Reports the results from `motif_search` in a dictionary containing enyzmes as keys and a list containing the corresponding motif sequence, cut site, observed motif counts, and the starting indeces of all observed motif counts as the values.
- Refactored `main.py`:
    - Updated function calls in `main()` to report the resulting dictionary from `Motifs.get_motif_results()`


Date: 03/25/2026
- Refactored `main.py`:
    - Create database function updated to conditional
    - User arguments updated to include: 
        - --interface: optional terminal interface
        - --enzymes: user list of enzymes to map with
        - --output: output filename
    - Added argument validation
    - Refactored main function calls to pass user arguments to input classes
- Refactored `input.py` module:
    - Removed error handling from `Sequence` class (handled with argparse now)
    - Replaced `Enzymes.load_RE()` function, database is read with csv.dictReader object at the filtering function `Enzymes.filter_enzymes()`instead 
    - Updated user interface to be optional
    - Updated default enzyme list to use constant defined in `main.py`


Date: 03/19/2026
- updated relative filepaths to run main.py from root project directory.
- Added `environment-alternative.yml` to resolve cross compatibility conflicts. 
- Refactored `main.py`: 
    1. Moved CLI header to `input.py` module to allow User interface as optional argument.
    2. Added `ArgParse` functionality to simplify CLI.  

## Testing Instructions:
Date: 03/25/2026

### Installation

For linux:
```bash
conda env create -f environment.yml
```
For MacOS:
```bash
conda env create -f environment_mac.yml
```
Conda will automatically create an environment named finalproj with all the specified packages and versions.

**Note:** if environment does not install on MacOS, install dependencies listed in `dependencies.txt` manually and run step 3. using installed python version. 

### Usage

1. Activate the environment:
```bash
conda activate finalproj
```
2. Run following command to test:
```bash
python src/main.py
```
#### Command-Line Arguments:
| Argument                | Description                                  | Default         |
| ----------------------- | -------------------------------------------- | --------------- |
| `-s`. `--sequence_filepath`| Plasmid sequence filepath | inputs/test/pUC19.fasta         |
| `-o`, `--output`        | Output filepath                           | results/plasmid_map|
| `-e`, `--enzymes`       | User list of restriction enzymes             | (list of top 20)|
| `--interface`           | Optional User interface option               | False|

Expected Output:

- Prints Plasmid object (list[header, sequence]) and enzyme object (dict{enzyme: [motif, cut]})

- --interface flag option displays Terminal User Interface



## Overview:
REcut is a python3 tool that takes a plasmid sequence file as input and generates an annotated sequence map with Restriction Enzyme cut sites. 

The `src` file contains the following:
1. `database` : Contains a CSV file of restriction enzymes and recognition sequence motifs. 
2. `main.py`: Contains the logic for the program.
    - Loads the user input sequence and restriction enzymes database. 
    - Parses the plasmid sequence with each enzyme motif to identify cut locations.
    - Produces output containing an annotated map of the sequence with the enzyme cut sites in the `results` folder.
3. `motif_id_lib`: Contains modules designated for each function within our program. These modules are:
    - motif_locator.py: Group member responsible- Bobby/Kayla
    - input.py: Group member responsible- Michael
    - output.py: Group member responsible- Michael
    - csv_output.py: Group member responsible- Bobby/Kayla

## Troubleshooting:
1. User inputs incorrect file path and incorrect file type - Check that the file/file path exists, is readable, and is the correct file format (FASTA or GenBank file)
2. No matches found and updated dictionary returns None - Error will occur and state that no matches were identified and exists program

## References:
Bioinformatics.org. (n.d.). IUPAC codes for nucleotides. https://www.bioinformatics.org/sms/iupac.html
NumPy Developers. (n.d.). NumPy reference: Routines. https://numpy.org/doc/stable/reference/routines.html
NumPy Developers. (n.d.). numpy.lib.stride_tricks.sliding_window_view. https://numpy.org/doc/stable/reference/generated/numpy.lib.stride_tricks.sliding_window_view.html

